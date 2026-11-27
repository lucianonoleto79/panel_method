"""Gera um PDF A4 a partir de um .md de trabalho_cilindro/ (padrao: enunciado.md).

Usa pandoc se estiver instalado; senao, o pacote Python 'markdown' (pip install
markdown). Em ambos os casos o PDF e impresso por um navegador Chromium sem
janela (Edge ou Chrome).

    python trabalho_cilindro/gerar_pdf.py [enunciado.md | gabarito.md]
"""
import html
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

AQUI = os.path.dirname(os.path.abspath(__file__))
NAVEGADORES = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "google-chrome", "chromium", "chromium-browser", "msedge",
]


def para_html(md_path):
    texto = open(md_path, encoding="utf-8").read()
    if shutil.which("pandoc"):
        return subprocess.run(["pandoc", md_path, "-s", "--css", "estilo.css", "--embed-resources",
                               "--resource-path", AQUI], capture_output=True, text=True, check=True,
                              encoding="utf-8").stdout
    import markdown
    meta = {}
    m = re.match(r"^---\n(.*?)\n---\n", texto, re.S)
    if m:
        for linha in m.group(1).splitlines():
            k, _, v = linha.partition(":")
            meta[k.strip()] = v.strip().strip('"')
        texto = texto[m.end():]
    corpo = markdown.markdown(texto, extensions=["tables", "fenced_code", "attr_list", "sane_lists"])
    css = open(os.path.join(AQUI, "estilo.css"), encoding="utf-8").read()
    cab = f'<h1 class="title">{html.escape(meta.get("title", ""))}</h1>'
    if meta.get("subtitle"):
        cab += f'<p class="subtitle">{html.escape(meta["subtitle"])}</p>'
    return (f'<!doctype html><html lang="{meta.get("lang", "pt-BR")}"><head><meta charset="utf-8">'
            f'<title>{html.escape(meta.get("title", ""))}</title><style>{css}</style></head><body>'
            f'<header id="title-block-header">{cab}</header>{corpo}</body></html>')


def main():
    md = os.path.join(AQUI, sys.argv[1] if len(sys.argv) > 1 else "enunciado.md")
    pdf = os.path.splitext(md)[0] + ".pdf"
    nav = next((n for n in NAVEGADORES if os.path.exists(n) or shutil.which(n)), None)
    if nav is None:
        sys.exit("nenhum navegador Chromium encontrado (Edge/Chrome); abra o HTML e imprima em PDF")
    # o executável do Edge no Windows devolve o controle antes de gravar o PDF:
    # espera o arquivo aparecer e parar de crescer antes de apagar o HTML
    tmp = tempfile.mkdtemp()
    h = os.path.join(tmp, "doc.html")
    open(h, "w", encoding="utf-8").write(para_html(md))
    if os.path.exists(pdf):
        os.remove(pdf)
    subprocess.run([nav, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    f"--user-data-dir={os.path.join(tmp, 'perfil')}",
                    f"--print-to-pdf={pdf}", "file:///" + h.replace(os.sep, "/")],
                   capture_output=True, timeout=120)
    tam = -1
    for _ in range(120):
        time.sleep(0.5)
        if os.path.exists(pdf) and os.path.getsize(pdf) == tam and tam > 0:
            break
        tam = os.path.getsize(pdf) if os.path.exists(pdf) else -1
    else:
        sys.exit(f"o PDF não foi gerado; o HTML ficou em {h}")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"{os.path.relpath(pdf)} gerado")


if __name__ == "__main__":
    main()

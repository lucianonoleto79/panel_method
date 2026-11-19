"""Configuracao comum dos exemplos.

Cada exemplo chama ``start("<nome_do_caso>")``: os resultados (figuras e a
saida de texto ``saida.txt``) sao gravados em ``results/<nome_do_caso>/``.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
RESULTS = os.path.join(ROOT, "results")
_case_dir = RESULTS


class _Tee:
    """Escreve simultaneamente no terminal e no arquivo de saida."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


def start(case):
    """Cria results/<case>/ e passa a gravar o stdout em results/<case>/saida.txt."""
    global _case_dir
    _case_dir = os.path.join(RESULTS, case)
    os.makedirs(_case_dir, exist_ok=True)
    log = open(os.path.join(_case_dir, "saida.txt"), "w", encoding="utf-8")
    sys.stdout = _Tee(sys.__stdout__, log)
    return _case_dir


def out(name):
    """Caminho de um arquivo de resultado dentro da pasta do caso atual."""
    return os.path.join(_case_dir, name)


def rel(path):
    """Caminho relativo a raiz do repositorio (para mensagens)."""
    return os.path.relpath(path, ROOT)

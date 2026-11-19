"""Configuracao comum dos exemplos."""
import os
import sys

import matplotlib

matplotlib.use("Agg")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
RESULTS = os.path.join(ROOT, "results")
os.makedirs(RESULTS, exist_ok=True)


def out(name):
    return os.path.join(RESULTS, name)

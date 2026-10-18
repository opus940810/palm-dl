"""Shared configuration for the PalmScan project."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
EMB = ROOT / "data" / "embeddings"
MODELS = ROOT / "models"
for _p in (DATA_PROC, RESULTS, FIGURES, EMB, MODELS):
    _p.mkdir(parents=True, exist_ok=True)

HALF = 15
WINDOW = 2 * HALF + 1
SEED = 42
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
N_FOLDS = 5

# Kyte-Doolittle hydropathy, normalised to [-1, 1]
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": 0.0,
}
# van der Waals volume, normalised later
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0, "X": 130.0,
}
# flexibility (Vihinen), higher = more flexible
FLEX = {
    "A": 0.36, "R": 0.53, "N": 0.46, "D": 0.51, "C": 0.35, "Q": 0.49, "E": 0.50,
    "G": 0.54, "H": 0.32, "I": 0.46, "L": 0.37, "K": 0.47, "M": 0.30, "F": 0.31,
    "P": 0.51, "S": 0.51, "T": 0.44, "W": 0.31, "Y": 0.42, "V": 0.39, "X": 0.44,
}
COMPARTMENTS = ["cytoplasmic", "extracellular", "lumenal", "signal_peptide", "other", "unknown"]

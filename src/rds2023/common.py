"""Small, repository-local plotting and validation utilities."""
from pathlib import Path
import hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = np.array([[145, 14, 38], [255, 70, 111], [246, 178, 147],
                   [109, 173, 209], [49, 54, 145]]) / 256


LABELS = [r"$B_{1h}$", r"$B_{2h}$", r"$B_{3h}$", r"$B_{6h}$", "Cath"]


def comparison(name, actual, expected, *, atol=1e-10, rtol=1e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    if a.shape != b.shape:
        return {"name": name, "pass": False, "shape": list(a.shape), "reference_shape": list(b.shape)}
    finite = np.isfinite(a).all() and np.isfinite(b).all()
    return {"name": name, "pass": bool(finite and np.allclose(a, b, atol=atol, rtol=rtol)),
            "max_abs_error": float(np.max(np.abs(a-b))) if finite else None,
            "atol": atol, "rtol": rtol}


def savefig(fig, destination, name):
    fig.savefig(destination / f"{name}.png", dpi=180)
    fig.savefig(destination / f"{name}.pdf")
    plt.close(fig)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

"""MATLAB and Dantec readers with explicit axis and component conventions."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import h5py
import numpy as np
from scipy.io import loadmat, whosmat

CASES = {"1b": "1block", "2b": "2block", "3b": "3block", "6b": "6block", "ca": "cathedral"}


def premultiplied_spectrum(pair, variance, *, frequency_column=0):
    """Read one explicitly ordered (frequency, PSD) or (PSD, frequency) pair."""
    pair = np.asarray(pair)
    if pair.ndim != 2 or pair.shape[1] != 2 or frequency_column not in (0, 1):
        raise ValueError("Spectrum requires two columns and an explicit frequency column")
    frequency, density = pair[:, frequency_column], pair[:, 1-frequency_column]
    if (len(frequency) < 2 or not np.all(np.isfinite(pair)) or frequency[0] < 0
            or not np.all(np.diff(frequency) > 0) or np.any(density < 0)
            or not np.isfinite(variance) or variance <= 0):
        raise ValueError("Invalid frequency ordering, PSD or reference variance")
    positive = frequency > 0
    return frequency[positive], frequency[positive]*density[positive]/variance


def one_mat(path: Path) -> np.ndarray:
    """Small single-variable MAT only; refuse ambiguous variables."""
    if path.stat().st_size > 100_000_000:
        raise ValueError(f"Use a streaming reader for large MAT file: {path.name}")
    if h5py.is_hdf5(path):
        with h5py.File(path, "r") as f:
            keys = [k for k in f if not k.startswith("#")]
            if len(keys) != 1:
                raise ValueError(f"Expected one variable in {path.name}: {keys}")
            d = f[keys[0]]
            if d.size * d.dtype.itemsize > 100_000_000:
                raise ValueError("Array too large for small MAT reader")
            return d[()].T
    info = whosmat(path)
    if any(np.prod(shape)*8 > 100_000_000 for _, shape, _ in info):
        raise ValueError("Decompressed array too large for small MAT reader")
    values = [v for k, v in loadmat(path).items() if not k.startswith("__")]
    if len(values) != 1:
        raise ValueError(f"Expected one variable in {path.name}")
    return values[0]


def dantec(path: Path) -> tuple[dict, np.ndarray]:
    """Read DSExport.TAB; header numbers must never be treated as data."""
    with path.open(encoding="utf-8", errors="replace") as f:
        header = []
        for line in f:
            header.append(line.rstrip())
            if line.strip() == ">>*DATA*<<":
                break
        else:
            raise ValueError(f"No Dantec DATA marker: {path}")
        columns = next(f).strip().split("\t")
        data = np.loadtxt(f, delimiter="\t", ndmin=2)
    grid = re.search(r"Width=(\d+), Height=(\d+)", "\n".join(header))
    if grid is None:
        raise ValueError(f"No grid size: {path}")
    nx, nz = map(int, grid.groups())
    if data.shape != (nx * nz, len(columns)):
        raise ValueError(f"Dantec grid/column mismatch: {path.name}: {data.shape}")
    return {"nx": nx, "nz": nz, "columns": columns, "header": header}, data


@dataclass
class PlaneStatistics:
    """All fields are (z, horizontal); velocities SI in physical U,V,W order."""
    plane: str
    horizontal: np.ndarray
    z: np.ndarray
    fields: dict[str, np.ndarray]
    source: Path

    def crop(self, start: int, stop: int) -> "PlaneStatistics":
        return PlaneStatistics(self.plane, self.horizontal[:, start:stop], self.z[:, start:stop],
                               {k: a[:, start:stop] for k, a in self.fields.items()}, self.source)


def statistics_path(root: Path, case: str, plane: str) -> Path:
    folder = root / "ORIGINAL" / f"RdS_PIV_{plane}_ORIGINAL" / CASES[case]
    pattern = "*_stats.*.txt" if plane == "XZ" else "*_YZ_stats_NOV.*.txt"
    paths = sorted(folder.glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"Expected one statistics export for {case}/{plane}, got {len(paths)}")
    return paths[0]


def read_statistics(path: Path, plane: str) -> PlaneStatistics:
    if plane not in ("XZ", "YZ"):
        raise ValueError("plane must be XZ or YZ")
    meta, data = dantec(path)
    shape = (meta["nz"], meta["nx"])
    # Zero-based columns of the Dantec export; signs refer to PHYSICAL axes.
    if plane == "XZ":
        mapping = {"U": (2, 1), "V": (4, -1), "W": (3, 1), "sigmaU": (5, 1),
                   "sigmaV": (7, 1), "sigmaW": (6, 1), "uw": (9, 1),
                   "uv": (10, -1), "vw": (11, -1), "skewU": (15, 1), "kurtU": (18, 1)}
        offset = 0.155
    else:
        mapping = {"U": (4, 1), "V": (2, 1), "W": (3, 1), "sigmaU": (7, 1),
                   "sigmaV": (5, 1), "sigmaW": (6, 1), "uw": (11, 1),
                   "uv": (10, 1), "vw": (9, 1), "skewU": (17, 1), "kurtU": (20, 1)}
        offset = 0.175
    fields = {key: sign * data[:, col].reshape(shape) for key, (col, sign) in mapping.items()}
    fields["k"] = 0.5 * sum(fields[k] ** 2 for k in ("sigmaU", "sigmaV", "sigmaW"))
    fields["N"] = data[:, 21].reshape(shape)
    fields["status"] = data[:, -1].reshape(shape)
    # Legacy exports retain masked zeros. Preserve them for reproduction and
    # expose N/status, rather than silently replacing or counting them as valid.
    return PlaneStatistics(plane, data[:, 0].reshape(shape) / 1000,
                           data[:, 1].reshape(shape) / 1000 + offset, fields, path)


def instantaneous_mat(root: Path, case: str, component: str, plane: str = "XZ") -> Path:
    folder = root / "ORIGINAL" / f"RdS_PIV_{plane}_ORIGINAL" / CASES[case]
    middle = "" if case == "1b" and plane == "XZ" else f"_{plane}"
    return folder / f"{component}_contour_PIV{middle}_{case}.mat"


def mat_blocks(path: Path, variable: str, plane: str, batch: int = 100):
    """Yield (time,z,horizontal) without materializing full velocity arrays."""
    if batch < 1:
        raise ValueError("batch must be positive")
    with h5py.File(path, "r") as f:
        d = f[variable]
        if plane == "XZ":
            for start in range(0, d.shape[2], batch):
                yield d[:, :, start:start + batch].transpose(2, 1, 0)
        elif plane == "YZ":
            for start in range(0, d.shape[0], batch):
                yield d[start:start + batch, :, :]
        else:
            raise ValueError("Unknown plane")

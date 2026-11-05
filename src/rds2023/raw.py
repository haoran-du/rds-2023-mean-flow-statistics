"""Full-sample XZ recomputation from archived instantaneous velocity MATs."""
from pathlib import Path
import json
import time
import numpy as np

from .io import CASES, dantec, instantaneous_mat, mat_blocks, one_mat, read_statistics, statistics_path
from .common import comparison
from .statistics import Moments


def text_blocks(root, case, batch=100):
    folder = root / "ORIGINAL/RdS_PIV_XZ_ORIGINAL" / CASES[case]
    paths = sorted(folder.glob(f"{CASES[case]}_.*.txt"))
    if len(paths) != 10000:
        raise ValueError(f"Expected 10000 raw snapshot files, got {len(paths)}")
    for start in range(0, len(paths), batch):
        frames = []
        for path in paths[start:start+batch]:
            meta, a = dantec(path)
            a = a[:, :3].reshape(meta["nz"], meta["nx"], 3)[:, 6:43, :]
            frames.append(np.stack([a[..., 0], -a[..., 2], a[..., 1]], axis=-1))
        yield np.stack(frames)


def validate_xz(root: Path, output: Path, cases=None, batch=100, source_type="mat"):
    root, output = root.resolve(strict=True), output.resolve()
    if root == output or root in output.parents:
        raise ValueError("Output cannot be in source data")
    output.mkdir(parents=True, exist_ok=True)
    reports = []
    for case in cases or CASES:
        start = time.monotonic()
        # Original TXT/profile/flux crop is 7:43 in MATLAB. Archived MATs
        # match 7:43 for 1b but 8:44 for the four other cases: verified against
        # first raw snapshots, not inferred from file names or stale scripts.
        crop_start = 6 if source_type == "txt" or case == "1b" else 7
        source = read_statistics(statistics_path(root, case, "XZ"), "XZ").crop(crop_start, crop_start+37)
        dx = source.horizontal[0, 1]-source.horizontal[0, 0]
        accum = Moments()
        flux = np.zeros((2, 3))
        readers = [mat_blocks(instantaneous_mat(root, case, c), f"{c}_contour", "XZ", batch) for c in "uvw"]
        blocks = (np.stack(arrays, axis=-1) for arrays in zip(*readers, strict=True)) if source_type == "mat" else text_blocks(root, case, batch)
        first = None
        for a in blocks:
            u, v, w = (a[..., j] for j in range(3))
            if first is None:
                first = a[0].copy()
            accum.update(a)
            for i, q in enumerate((w[:, 48, :], v[:, :48, :])):
                flux[i] += np.array([q.sum(), np.maximum(q, 0).sum(), np.minimum(q, 0).sum()]) * dx
        if accum.n != 10000:
            raise ValueError(f"Expected 10000 XZ snapshots, got {accum.n}")
        flux /= accum.n
        cov = accum.covariance(ddof=1)
        std = np.sqrt(np.maximum(np.diagonal(cov, axis1=-2, axis2=-1), 0))
        checks = []
        for j, name in enumerate("UVW"):
            checks.append(comparison(f"{case}/instantaneousMAT/{name}/Dantec", accum.mean[..., j], source.fields[name], atol=2e-6, rtol=2e-6))
            checks.append(comparison(f"{case}/instantaneousMAT/sigma{name}/Dantec", std[..., j], source.fields[f"sigma{name}"], atol=2e-6, rtol=2e-6))
        checks.append(comparison(f"{case}/instantaneousMAT/uw/Dantec", cov[..., 0, 2], source.fields["uw"], atol=2e-6, rtol=2e-6))
        # Verify MAT component mapping/crop against original text, not only
        # against another stored statistical product from the same chain.
        folder = root / "ORIGINAL/RdS_PIV_XZ_ORIGINAL" / CASES[case]
        text_paths = list(folder.glob(f"{CASES[case]}_.*.000000.txt"))
        if len(text_paths) != 1:
            raise ValueError(f"Ambiguous raw snapshot: {case}")
        meta, table = dantec(text_paths[0])
        a = table[:, :3].reshape(meta["nz"], meta["nx"], 3)[:, crop_start:crop_start+37, :]
        physical = np.stack((a[..., 0], -a[..., 2], a[..., 1]), axis=-1)
        checks.append(comparison(f"{case}/firstMAT=firstRawTXT", first, physical))
        if crop_start == 6:
            for i, c in enumerate("wv"):
                for j, direction in enumerate(("total", "p", "n")):
                    path = root / "POSTPROCESSING/RdS_PIV_XZ_MASSFLOWRATE" / f"Q{c}_{direction}_{case}.mat"
                    checks.append(comparison(f"{case}/Q{c}_{direction}/full-recompute", flux[i, j], one_mat(path).item()))
        np.savez_compressed(output / f"{case}_XZ_raw_moments.npz", mean=accum.mean, covariance=cov,
                            n=accum.n, horizontal=source.horizontal, z=source.z, legacy_flux=flux)
        report = {"case": case, "snapshots": accum.n, "checks": checks, "seconds": time.monotonic()-start,
                  "source_type": source_type, "variance_ddof": 1, "masked_zero_policy": "preserved to reproduce legacy exports",
                  "crop_matlab": f"horizontal {crop_start+1}:{crop_start+37}",
                  "flux_validation": "compared" if crop_start == 6 else "not comparable: MAT crop differs from legacy flux crop; use --source txt",
                  "legacy_roof_z_m": float(source.z[48, 0]),
                  "legacy_dx_m": float(dx), "flux_signed": flux.tolist()}
        (output / f"{case}_validation.json").write_text(json.dumps(report, indent=2))
        reports.append(report)
        print(f"{case}: {accum.n} snapshots; {sum(c['pass'] for c in checks)}/{len(checks)} checks pass; {report['seconds']:.1f}s", flush=True)
    summary = {"cases": reports, "passed": sum(c["pass"] for r in reports for c in r["checks"]),
               "failed": [c for r in reports for c in r["checks"] if not c["pass"]]}
    (output / "validation.json").write_text(json.dumps(summary, indent=2))
    return {"passed": summary["passed"], "failed": summary["failed"]}

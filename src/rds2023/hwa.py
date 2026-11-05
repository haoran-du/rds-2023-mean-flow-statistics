"""2023 HWA reproduction, explicit historical corrections, read-only binary IO."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.signal.windows import hamming
from scipy.fft import rfft, rfftfreq
from .io import CASES, one_mat
from .common import comparison


def corrected(data, case, start, stop, channels, *, b6_policy="saved-zero"):
    """Historical corrections applied to a COPY, with B6h version explicit.

    MATLAB ranges are one-based and inclusive. Do not use these corrections
    for a new dataset; they are specific historical interventions.
    saved-zero numerically reconstructs saved 2023 statistics AND spectra;
    recorded-replacement follows the available .m instead. Zero all selected
    channels, including camera triggers, before deriving the analysis window.
    """
    if b6_policy not in ("saved-zero", "recorded-replacement"):
        raise ValueError("Unknown B6h historical preprocessing policy")
    a = np.array(data[start:stop, :][:, channels], copy=True)
    if case == "3b":
        lo, hi = 3782999, 3787500
        left, right = max(start, lo), min(stop, hi)
        if left < right:
            a[left-start:right-start] = 0
    elif case == "6b":
        lo, hi, source = 19375499, 19384500, 19385999
        left, right = max(start, lo), min(stop, hi)
        if left < right:
            a[left-start:right-start] = (0 if b6_policy == "saved-zero" else
                data[source+left-lo:source+right-lo, :][:, channels])
    return a


def welch_stream(data, case, channels, mean, nperseg=32768, fs=15000, *, b6_policy="saved-zero"):
    """MATLAB scalar-window pwelch: symmetric Hamming, half overlap.

    Global mean removed once, NO additional per-segment detrending. The
    explicit convention is verified numerically against saved spectra.
    """
    window = hamming(nperseg, sym=True)
    total = np.zeros((nperseg//2+1, len(channels)))
    count = 0
    for start in range(0, len(data)-nperseg+1, nperseg//2):
        a = corrected(data, case, start, start+nperseg, channels, b6_policy=b6_policy)-mean
        transformed = rfft(a*window[:, None], axis=0)
        total += (transformed*transformed.conj()).real
        count += 1
    if not count:
        raise ValueError("Recording shorter than Welch segment")
    total /= count*fs*np.sum(window**2)
    total[1:-1] *= 2
    return rfftfreq(nperseg, 1/fs), total


def validate_hwa(root: Path, output: Path, cases=None, batch=200000, *, b6_policy="saved-zero"):
    root, output = root.resolve(strict=True), output.resolve()
    if root == output or root in output.parents:
        raise ValueError("Output cannot be in source data")
    output.mkdir(parents=True, exist_ok=True)
    reports = []
    def read_block(data, case, start, stop, channels):
        return corrected(data, case, start, stop, channels, b6_policy=b6_policy)
    for case in cases or CASES:
        started = time.monotonic()
        folder_case = "cathedrale" if case == "ca" else CASES[case]
        folder = root / "ORIGINAL/RdS_HWA_XZ_ORIGINAL" / f"CUBES-25perc_0700_{folder_case}/data"
        files = list(folder.glob("*.f64LEx15"))
        if len(files) != 1 or files[0].stat().st_size % 120:
            raise ValueError(f"Ambiguous or misaligned HWA recording for {case}")
        data = np.memmap(files[0], mode="r", dtype="<f8", shape=(files[0].stat().st_size//120, 15))
        triggers = []
        for i in range(0, len(data), batch):
            # The legacy code uses every sample ==1, not rising edges.
            block = read_block(data, case, i, min(i+batch, len(data)), [14])
            triggers.extend((np.flatnonzero(block[:, 0] == 1)+i).tolist())
        if len(triggers) < 10000:
            raise ValueError(f"Insufficient camera trigger samples: {case}")
        lo, hi = triggers[0]+43, triggers[9999]+44
        if hi > len(data):
            raise ValueError("Trigger-aligned segment exceeds recording")
        sums = np.zeros(10)
        full_sums = np.zeros(5)
        channels = [6, 2, 0, 4, 8]
        for i in range(0, len(data), batch):
            full_sums += read_block(data, case, i, min(i+batch, len(data)), channels).sum(0)
        full_mean = full_sums/len(data)
        for i in range(lo, hi, batch):
            sums += read_block(data, case, i, min(i+batch, hi), list(range(10))).sum(0)
        mean = sums/(hi-lo)
        moments = np.zeros((3, 10))
        full_m2 = np.zeros(5)
        for i in range(lo, hi, batch):
            a = read_block(data, case, i, min(i+batch, hi), list(range(10)))-mean
            square = a*a
            moments += np.stack([square.sum(0), (square*a).sum(0), (square*square).sum(0)])
        for i in range(0, len(data), batch):
            a = read_block(data, case, i, min(i+batch, len(data)), channels)-full_mean
            full_m2 += (a*a).sum(0)
        moments /= hi-lo
        results = {"meanU": mean, "stdU": np.sqrt(moments[0]*(hi-lo)/(hi-lo-1)),
                   "skewU": moments[1]/moments[0]**1.5, "kurtU": moments[2]/moments[0]**2}
        checks = []
        for field, value in results.items():
            path = root / "POSTPROCESSING/RdS_HWA_XZ_STATISTICS" / f"RdS_HWA_XZ_{field}_{case}.mat"
            checks.append(comparison(f"{case}/HWA/{field}", value, one_mat(path).ravel(), atol=1e-9, rtol=1e-9))
        freq, psd = welch_stream(data, case, channels, full_mean, b6_policy=b6_policy)
        folder = root / "POSTPROCESSING/RdS_HWA_XZ_PSD"
        saved = one_mat(folder / f"RdS_HWA_XZ_PSD_{case}.mat")
        checks.append(comparison(f"{case}/HWA/PSD_frequency", freq, saved[:, 0]))
        checks.append(comparison(f"{case}/HWA/PSD", psd, saved[:, 1::2], atol=1e-10, rtol=1e-9))
        checks.append(comparison(f"{case}/HWA/PSD_energy", full_m2/len(data), one_mat(folder / f"RdS_HWA_XZ_PSD_ENERGY_{case}.mat").ravel(), atol=1e-9, rtol=1e-9))
        np.savez_compressed(output / f"{case}_hwa.npz", **results, frequency=freq, psd=psd, variance_full=full_m2/len(data))
        report = {"case": case, "checks": checks, "samples": len(data), "camera_trigger_samples": len(triggers),
                  "statistics_slice_python": [lo, hi], "psd_samples": len(data), "window": "symmetric Hamming 32768; 50% overlap; global demeaning only",
                  "b6_policy": b6_policy,
                  "spike_policy": "3b zeroing; B6h policy explicit; changes in memory only; triggers derived AFTER correction",
                  "seconds": time.monotonic()-started}
        (output / f"{case}_validation.json").write_text(json.dumps(report, indent=2))
        reports.append(report)
        print(f"HWA {case}: {sum(c['pass'] for c in checks)}/{len(checks)} checks pass, {report['seconds']:.1f}s", flush=True)
        del data
    summary = {"b6_policy": b6_policy, "cases": reports, "passed": sum(c["pass"] for r in reports for c in r["checks"]),
               "failed": [c for r in reports for c in r["checks"] if not c["pass"]]}
    (output / "validation.json").write_text(json.dumps(summary, indent=2))
    return {"passed": summary["passed"], "failed": summary["failed"]}

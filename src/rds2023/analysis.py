"""2023 adapters: distinguish recovered exports from raw recalculation."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .common import sha256
from .io import CASES, one_mat, premultiplied_spectrum, read_statistics, statistics_path

from .common import COLORS, LABELS, comparison, savefig


def reproduce(root: Path, output: Path, normalization="paper"):
    if normalization not in ("paper", "legacy"):
        raise ValueError("normalization must be paper or legacy")
    root, output = root.resolve(strict=True), output.resolve()
    if output == root or root in output.parents:
        raise ValueError("Do not write into the original data tree")
    output.mkdir(parents=True, exist_ok=True)
    post = root / "POSTPROCESSING"
    profiles_dir = post / "RdS_PIV_XZ_STATISTICS/RdS_PIV_XZ_STATISTICS_PROFILES"
    checks, provenance, xz, yz, profiles = [], [], {}, {}, {}
    mapping = {"U": "U", "V": "V", "W": "W", "sigmaU": "stdU", "sigmaV": "stdV",
               "sigmaW": "stdW", "uw": "upwp", "uv": "upvp", "vw": "vpwp", "k": "tke"}
    for case in CASES:
        for plane, collection in (("XZ", xz), ("YZ", yz)):
            path = statistics_path(root, case, plane)
            data = read_statistics(path, plane)
            collection[case] = data
            provenance.append({"path": path.relative_to(root).as_posix(), "sha256": sha256(path)})
            cropped = data.crop(7, 44) if plane == "XZ" else data
            folder = (post / "RdS_PIV_XZ_STATISTICS/RdS_PIV_XZ_STATISTICS_CONTOURS" if plane == "XZ"
                      else post / "RdS_PIV_YZ_STATISTICS/RdS_PIV_YZ_STATISTICS_NOV6")
            for field, saved in mapping.items():
                path = folder / f"RdS_PIV_{plane}_{saved}_{case}_STATISTICS_CONTOUR.mat"
                if path.exists():
                    checks.append(comparison(f"{case}/{plane}/{field}/export-to-MAT", cropped.fields[field], one_mat(path).T))
            np.savez_compressed(output / f"{case}_{plane}_statistics_SI.npz",
                                horizontal=data.horizontal, z=data.z, **data.fields)
        # block*_test.m retains MATLAB columns 7:43, whereas the contours
        # use 8:44. Reproduce both intentionally, never silently shift a grid.
        data = xz[case].crop(6, 43)
        profile = {k: v.mean(axis=1) for k, v in data.fields.items()}
        profile["z"] = data.z[:, 0]
        profile["k_from_averaged_std_legacy"] = 0.5 * sum(profile[k]**2 for k in ("sigmaU", "sigmaV", "sigmaW"))
        profiles[case] = profile
        for key, old in (("U", "meanU"), ("V", "meanV"), ("W", "meanW"),
                         ("sigmaU", "sigmaU"), ("sigmaV", "sigmaV"), ("sigmaW", "sigmaW"), ("uw", "upwp")):
            checks.append(comparison(f"{case}/profile/{key}", profile[key], one_mat(profiles_dir / f"PIV_{old}_{case}.mat").ravel()))
        with (output / f"{case}_profiles_SI.csv").open("w", newline="") as f:
            keys = ["z", "U", "V", "W", "sigmaU", "sigmaV", "sigmaW", "uw", "k", "k_from_averaged_std_legacy"]
            writer = csv.writer(f)
            writer.writerow(keys)
            writer.writerows(zip(*(profile[k] for k in keys)))
    plot_profiles(profiles, profiles_dir, output)
    ustar, h2 = ((0.07 * 5.8, 0.116) if normalization == "paper" else (0.064 * 5.8, 0.1165))
    plot_contours(yz, "U", ustar, h2, output, "fig05_YZ_mean", r"$\overline{U}/u_*$", vectors=True)
    plot_contours(yz, "k", ustar**2, h2, output, "fig06_YZ_TKE", r"$\overline{k}/u_*^2$")
    cropped = {case: data.crop(7, 44) for case, data in xz.items()}
    plot_contours(cropped, "V", ustar, h2, output, "fig08_XZ_along_canyon", r"$\overline{V}/u_*$", vectors=True)
    plot_contours(cropped, "uw", -ustar**2, h2, output, "fig09_XZ_shear", r"$-\overline{u'w'}/u_*^2$")
    flux = recover_flux(post, output, checks)
    plot_hwa(post, yz, output, checks)
    report = {"status": "Dantec statistical-export reconstruction; HWA/flux loaded from saved MATLAB outputs",
              "normalization": normalization, "contour_u_star_m_s": ustar, "contour_h2_m": h2,
              "profile_u_star_m_s": 0.406, "profile_h2_m": 0.116,
              "raw_snapshots_recomputed": False, "checks": checks, "inputs": provenance,
              "flux": flux, "issues": [
                  "Profile crop MATLAB 7:43; contour crop 8:44; preserved, not reconciled.",
                  "Legacy contours use u*/Ue=0.064, h2=116.5 mm; paper/Table2 and profile plots use 0.07,116 mm.",
                  "Fig7 legacy TKE squares spatially averaged standard deviations; CSV also supplies mean of local TKE.",
                  "Legacy Fig10 uses Wc=74 mm, roof row49 and 48-row averaging; paper text Wc=73 mm.",
                  "HWA preprocessing edits spikes for 3b/6b; full raw-HWA validation is separate.",
                  "Contours here show measured domains, not manually drawn building outlines/optical masks."]}
    (output / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return {"checks": len(checks), "passed": sum(c["pass"] for c in checks),
            "failed": [c for c in checks if not c["pass"]], "output": str(output)}


def plot_profiles(profiles, reference_dir, output):
    fields = ["U", "W", "V", "sigmaU", "sigmaW", "sigmaV", "uw", "k_from_averaged_std_legacy"]
    labels = [r"$\langle\overline{U}\rangle/u_*$", r"$\langle\overline{W}\rangle/u_*$", r"$\langle\overline{V}\rangle/u_*$",
              r"$\langle\sigma_u\rangle/u_*$", r"$\langle\sigma_w\rangle/u_*$", r"$\langle\sigma_v\rangle/u_*$",
              r"$-\langle\overline{u'w'}\rangle/u_*^2$", r"$k_{\mathrm{legacy}}/u_*^2$"]
    ref = {}
    for key, name in zip(fields[:-1], ["U_mean", "W_mean", "V_mean", "U_std", "W_std", "v_std", "upwp_mean"]):
        ref[key] = one_mat(reference_dir / f"AR07_{name}.mat").ravel()
    ref[fields[-1]] = 0.5 * sum(ref[key]**2 for key in ("sigmaU", "sigmaW", "sigmaV"))
    ref_z = one_mat(reference_dir / "AR07_z_mean.mat").ravel() + 0.003
    fig, axes = plt.subplots(4, 2, figsize=(9, 12), layout="constrained")
    for i, (ax, field, label) in enumerate(zip(axes.flat, fields, labels)):
        norm = (0.406**2 if i >= 6 else 0.406) * (-1 if field == "uw" else 1)
        for color, (case, p), legend in zip(COLORS, profiles.items(), LABELS):
            ax.plot(p[field][:153]/norm, p["z"][:153]/0.116, color=color, label=legend)
        # AllPIVStatistics2.m uses the SAME u* for AR0.7 despite defining 0.38.
        ax.plot(ref[field][:100]/norm, ref_z[:100]/0.05, color="0.4", ls="--", label="AR0.7")
        ax.set(xlabel=label, ylabel=r"$z/h_*$", ylim=(0, 3))
        ax.set_title(f"({chr(97+i)})", loc="left")
        ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8)
    savefig(fig, output, "fig07_XZ_profiles_legacy_definition")


def plot_contours(data, key, norm, h2, output, name, label, vectors=False):
    is_xz = next(iter(data.values())).plane == "XZ"
    fig, axes = plt.subplots(1, 5, figsize=(10, 9), layout="constrained") if is_xz else plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
    lo = min(np.min(p.fields[key]/norm) for p in data.values())
    hi = max(np.max(p.fields[key]/norm) for p in data.values())
    levels = np.linspace(lo, hi, 41)
    for ax, p, legend in zip(axes.flat, data.values(), LABELS):
        x, z = p.horizontal[0]/h2, p.z[:, 0]/h2
        vals = p.fields[key]/norm
        image = ax.contourf(x, z, vals, levels=levels, cmap="RdBu_r", extend="both")
        if vectors:
            hor = "U" if p.plane == "XZ" else "V"
            if p.plane == "XZ":
                # streamplot requires equally spaced coordinates to roundoff.
                ax.streamplot(np.linspace(x[0], x[-1], len(x)), np.linspace(z[0], z[-1], len(z)),
                              p.fields[hor], p.fields["W"], density=.8, color="k", linewidth=.45, arrowsize=.6)
            else:
                ax.quiver(x[::8], z[::8], p.fields[hor][::8, ::8], p.fields["W"][::8, ::8], scale=22, width=.003)
        ax.set_title(legend)
        ax.set(xlabel=(r"$x/h_2$" if p.plane == "XZ" else r"$y/h_2$"), ylabel=r"$z/h_2$")
        ax.set_ylim((0, 3) if p.plane == "XZ" else (1.1, 3.7))
        ax.set_aspect("equal")
    if not is_xz:
        axes.flat[-1].axis("off")
    fig.colorbar(image, ax=list(axes.flat[:5]), label=label, shrink=.8)
    savefig(fig, output, name)


def recover_flux(post, output, checks):
    folder = post / "RdS_PIV_XZ_MASSFLOWRATE"
    result = {}
    fig, axes = plt.subplots(2, 1, figsize=(7, 7), layout="constrained")
    rows = []
    for axis, ax in zip(("w", "v"), axes):
        norm = 0.074 * 5.8 * (48 if axis == "v" else 1)
        values = []
        for case in CASES:
            q = np.array([one_mat(folder / f"Q{axis}_{direction}_{case}.mat").item() for direction in ("total", "p", "n")])
            checks.append(comparison(f"{case}/Q{axis}/net=positive+negative", q[0], q[1]+q[2]))
            j = q * 1.184 / norm
            values.append(j)
            result[f"{case}/{axis}"] = dict(zip(("net", "positive", "negative_signed"), j.tolist()))
            rows.append([case, axis, *q, *j])
        for i, (j, color, label) in enumerate(zip(values, COLORS, LABELS)):
            ax.bar(np.arange(3)+(i-2)*.15, j*np.array([1, 1, -1]), width=.15, color=color, label=label)
        ax.set(xticks=np.arange(3), xticklabels=["Net", "Positive", "Negative magnitude"], ylabel=rf"$j_{{m,{axis.upper()}}}$")
        ax.grid(axis="y", alpha=.2)
    axes[0].legend(ncol=5, fontsize=8)
    savefig(fig, output, "fig10_flux_recovered_MAT")
    with (output / "flux_legacy.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "component", "Qnet", "Qpositive", "Qnegative_signed", "jnet", "jpositive", "jnegative_signed"])
        writer.writerows(rows)
    return result


def plot_hwa(post, yz, output, checks):
    folder = post / "RdS_HWA_XZ_STATISTICS"
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    yh = np.linspace(-.1, .1, 5)
    differences = []
    for ax, field, suffix, scale, label in zip(axes.flat,
            ["U", "sigmaU", "skewU", "kurtU"], ["meanU", "stdU", "skewU", "kurtU"],
            [.406, .406, 1, 1], [r"$\overline{U}/u_*$", r"$\sigma_u/u_*$", r"$S_u$", r"$K_u$"]):
        difference = []
        for (case, p), color, legend in zip(yz.items(), COLORS, LABELS):
            y = p.horizontal[64, 58:115]
            v = p.fields[field][64, 58:115]
            checks.append(comparison(f"{case}/YZ300mm/{field}", v, one_mat(folder / f"RdS_PIV_YZ_{suffix}_300mm_{case}.mat").ravel()))
            if field == "kurtU":
                v = v + 3  # Dantec EXCESS kurtosis vs MATLAB Pearson kurtosis.
            h = one_mat(folder / f"RdS_HWA_XZ_{suffix}_{case}.mat").ravel()[[6, 2, 0, 4, 8]]
            difference.append(np.std(v[::14]-h, ddof=1))
            ax.plot(y/.116, v/scale, color=color, ls="--")
            ax.plot(yh/.116, h/scale, color=color, marker="o", ms=3, label=legend)
        ax.set(xlabel=r"$y/h_2$", ylabel=label)
        ax.grid(alpha=.2)
        differences.append(difference)
    axes[0, 0].legend(fontsize=8)
    savefig(fig, output, "fig02_HWA_saved_vs_PIV_recomputed")
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for i, (color, legend) in enumerate(zip(COLORS, LABELS)):
        ax.bar(np.arange(4)+(i-2)*.15, np.array(differences)[:, i], width=.15, color=color, label=legend)
    ax.set(xticks=np.arange(4), xticklabels=[r"$\overline{U}$", r"$\sigma_u$", r"$S_u$", r"$K_u$"],
           ylabel="Std. of HWA–PIV difference (original units)")
    ax.legend(ncol=5, fontsize=8)
    savefig(fig, output, "fig03_HWA_PIV_difference")
    folder = post / "RdS_HWA_XZ_PSD"
    fig, axes = plt.subplots(3, 2, figsize=(10, 9), layout="constrained")
    for case, color, legend in zip(CASES, COLORS, LABELS):
        p = one_mat(folder / f"RdS_HWA_XZ_PSD_{case}.mat")
        e = one_mat(folder / f"RdS_HWA_XZ_PSD_ENERGY_{case}.mat").ravel()
        for j, ax in enumerate(axes.flat[:5]):
            freq, spectrum = premultiplied_spectrum(p[:, 2*j:2*j+2], e[j])
            ax.semilogx(freq*.2/5.8, spectrum, color=color, label=legend)
    cube = one_mat(folder / "Cube_HWA_XZ_PSD.mat")
    energy = one_mat(folder / "Cube_HWA_XZ_PSD_ENERGY.mat").ravel()
    for j, ax in enumerate(axes.flat[:5]):
        # Cube stores (PSD, frequency), the reverse of the RdS pairs.
        freq, spectrum = premultiplied_spectrum(cube, energy[0], frequency_column=1)
        ax.semilogx(freq*.2/5.8, spectrum, color=".4", ls="--", label="Cubic array")
        ax.set(xlabel=r"$fW_b/U_e$", ylabel=r"$f\Phi_{uu}/\overline{u'^2}$", xlim=(1e-3, 1e2),
               title=rf"$y/h_2={yh[j]/.116:.2f}$")
        ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=7)
    axes.flat[-1].axis("off")
    savefig(fig, output, "fig04_premultiplied_spectra_saved_MAT")

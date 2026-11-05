RdS 2023: mean flow and turbulence statistics
============================================

Independent Python reproduction tools for Du, Savory & Perret (2023),
Building and Environment, https://doi.org/10.1016/j.buildenv.2023.110428.
This repository concerns this paper only. MATLAB and the 2024/2025 code
repositories are not runtime dependencies.

Install and run
---------------

From this repository directory, with Python 3.10 or newer::

    python -m venv .venv
    # Linux/macOS: source .venv/bin/activate
    # Windows PowerShell: .venv\Scripts\Activate.ps1
    python -m pip install -e ".[test]"
    python -m pytest -q
    rds2023 reproduce --data-root "/path/to/DATA" --output outputs/paper
    rds2023 reproduce --data-root "/path/to/DATA" --output outputs/legacy --normalization legacy
    rds2023 validate-xz --data-root "/path/to/DATA" --output outputs/raw_txt --source txt
    rds2023 validate-xz --data-root "/path/to/DATA" --output outputs/raw_mat --source mat
    rds2023 validate-hwa --data-root "/path/to/DATA" --output outputs/hwa --b6-policy saved-zero

Data are external and read-only. The same DATA layout used by the original
study is required; no drive letter or personal directory is hardcoded.
``python -m rds2023.cli`` is equivalent to the ``rds2023`` entry point.

Scope and verified baseline
---------------------------

* Numerical figures corresponding to Figs. 2–10, with simplified layouts,
  SI arrays, profile CSVs and explicit normalization: 165/165 checks.
* Five cases and all 50,000 original XZ TXT snapshots: 70/70 checks.
* Instantaneous XZ MATs with their actual case-specific crops: 46/46 checks.
* Five raw HWA recordings with the saved-result processing: 35/35 checks.

These are the validated pre-split data results, retained in
``docs/validation_summary.json``. Splitting the repository does not by itself
constitute another full-data rerun. Automated CI uses synthetic tests, not
the unpublished experimental data. See ``docs/SPLIT_VERIFICATION.rst`` for
checks performed after extraction.

B6h defaults to zeroing its historical anomalous interval, including trigger
channels, before identifying camera triggers. This reconstructs saved
statistics and spectra. ``--b6-policy recorded-replacement`` instead follows
the available .m version and retains its documented mismatch. Neither mode
changes original data. See ``docs/AUDIT_FINDINGS.rst``.

The uncertainty table and manual geometry/optical-mask overlays are not
fully reproduced. Historical crop, normalization and TKE-averaging
differences are explicit; do not transfer these settings unexamined to a
new experiment. The future idealised-campaign interface remains in the
original research workspace, outside this published-paper repository.

Contents and provenance
-----------------------

``src/rds2023`` contains only the statistics, HWA, raw-frame and figure
workflows. ``tests`` is data-independent. ``archive/matlab`` is a local-only
byte-exact source archive; its public hash index is in ``docs``. Raw arrays,
MATLAB originals and publisher PDFs are excluded from GitHub.

See ``THIRD_PARTY.rst`` for redistribution status and ``CITATION.cff`` for
this paper. Commands return status 1 when a validation fails; generated
figures are not by themselves evidence that all checks passed.

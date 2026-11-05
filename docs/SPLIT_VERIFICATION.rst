Post-split verification
=======================

2026-11-05. Independent Python 3.12 environment: 21 tests passed.
The command-line entry point is installed and runs without the other paper
packages or the original combined ``rds`` package.

Figure reconstruction: 165/165. Raw B6h HWA: 7/7. The all-case XZ and
HWA baseline checks remain inherited evidence, not newly claimed reruns.

The existing scientific checks and tolerances are unchanged. Source
function/class bodies were also compared by AST; see PROVENANCE.rst.
Numerical arrays, figures and raw data remain local. The accompanying
``split_validation.json`` contains check results only.

Each repository has its own Linux/Windows/macOS x Python 3.10/3.12 workflow.
Local success is not a claim of remote CI completion; the Actions tab shows
the actual results for the published commit.

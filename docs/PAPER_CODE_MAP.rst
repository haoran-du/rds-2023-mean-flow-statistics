Paper-to-code map
=================

Paths below are DATA-relative and retained in the local archive.

2023: mean flow and turbulence statistics
---------------------------------------

* Geometry/source: local paper_RdS.pdf; code does not recreate CAD geometry.
* Figs. 2–3: ``PREPROCESSING/RdS_HWA_XZ/RdS_HWA_XZ_PREPROCESSING.m``;
  ``PREPROCESSING/RdS_PIV_YZ/RdS_PIV_YZ_STATISTICS_CONTOURS_PREPROCESSING.m``;
  ``POSTPROCESSING/RdS_HWA_XZ_STATISTICS/RdS_HWA_XZ_STATISTICS_POSTPROCESSING.m``.
* Fig. 4: same HWA preprocessing; ``POSTPROCESSING/RdS_HWA_XZ_PSD/``.
* Figs. 5–6: YZ statistical exports -> preprocessing ->
  ``POSTPROCESSING/RdS_PIV_YZ_STATISTICS/RdS_PIV_YZ_STATISTICS_NOV6/``.
  OCT20 is retained as an earlier variant, not automatically substituted.
* Fig. 7: ``POSTPROCESSING/RdS_PIV_XZ_STATISTICS/RdS_PIV_XZ_STATISTICS_PROFILES/``;
  block*_test.m/cathedral_test.m -> saved PIV_*.mat -> AllPIVStatistics2.m.
  AR07_v_std.mat has lower-case v: the original Windows reference used V.
* Figs. 8–9: ``POSTPROCESSING/RdS_PIV_XZ_STATISTICS/RdS_PIV_XZ_STATISTICS_CONTOURS/``.
* Fig. 10: flux section of
  ``PREPROCESSING/RdS_PIV_XZ/RdS_PIV_XZ_STATISTICAL_ERROR_PREPROCESSING.m``
  -> Qw/Qv files -> ``POSTPROCESSING/RdS_PIV_XZ_MASSFLOWRATE/``.
* Table 3: statistical-error sections of that same preprocessing file;
  uncertainty methodology is explicitly not approved for automatic reuse.

Python: see this repository's ``src/rds2023`` package.


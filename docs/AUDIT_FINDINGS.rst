2023 paper: conventions and audit findings
=================================================

1. Coordinate/component conventions
----------------------------------

XZ Dantec U,V,W columns map to physical U,W,-V; the vertical calibration
offset is +155 mm. YZ columns map to physical V,W,U, with +175 mm offset.
Python arrays use (time,z,horizontal[,component]) and physical U,V,W order.
MATLAB v7.3/HDF5 storage dimensions are reversed; explicit transposes are
tested. The YZ coordinate named x_contour in old files is physically y.

2. One-column discrepancy, verified from original snapshots
---------------------------------------------------------

Published-profile and flux scripts retain MATLAB horizontal columns 7:43.
Contour scripts retain 8:44. The 1b instantaneous MAT matches 7:43;
2b/3b/6b/ca instantaneous MATs match 8:44, despite old per-case source
scripts showing 7:43. All original TXT snapshots therefore provide the
authoritative input for legacy flux reconstruction. The two MAT crops
must not be averaged or relabelled as the same coordinates.

3. Normalization and statistical definitions
-------------------------------------------

* Profile scripts and paper Table 2 use u*/Ue=0.07; contour scripts use
  0.064. h2 is 116 mm in the paper/profiles, 116.5 mm in several scripts.
  Both plotting modes are provided. This is not a decision that either
  archived plotting script necessarily generated the final publisher PDF.
* The Dantec XZ exported variances/covariances match the N-1 convention
  when checked against all original snapshots. Generic estimators expose
  ddof rather than hiding it.
* Fig. 7's saved plotting logic computes TKE from squared spatially
  averaged standard deviations. The mean of local TKE is different.
  Both values are exported; new-campaign analysis uses the latter.
* The flux script uses Wc=74 mm, roof row 49, dx-weighted sums and a
  48-row mean for along-canyon flux, while the paper reports Wc=73 mm
  and a continuous integral over h2. Exact discrete legacy values are
  preserved; this stencil is not copied to the new experiment.
* Positive/negative fluxes are separated on instantaneous velocities.
  The negative contribution is signed internally; plots show its magnitude.
* Dantec kurtosis is excess kurtosis in the HWA comparison; +3 is needed
  to compare with the MATLAB Pearson kurtosis output.
* The saved cubic-array reference spectrum uses (PSD, frequency), whereas
  RdS files use (frequency, PSD). Both are read with explicit column order
  and validated frequency monotonicity before plotting Fig. 4.
* Statistical exports retain zero-filled masked positions. Reproduction
  preserves them and exposes N/status. New data require an explicit mask.

4. B6h HWA version difference: numerically resolved
--------------------------------------------------

The available preprocessing source replaces samples 1291.7–1292.3 s with
1292.4–1293.0 s. Exact one-based inclusive indexing has been translated.
That implementation reproduces 1b/2b/3b/ca but initially did not reconstruct
the saved B6h products. A follow-up check on 2026-11-05 resolved this:
zeroing the SAME 1291.7–1292.3 s interval in ALL columns, then finding the
camera triggers, reconstructs all four statistics, frequency, PSD and
full-record variance. Maximum errors are about 1.6e-11 for the four
statistics, 1.5e-14 for PSD and 1.4e-12 for full-record variance.

Trigger handling matters: zeroing removes four camera-trigger samples.
The zero-based statistics interval changes from [1786539,22059882) to
[1786539,22067992). Zeroing velocities but keeping the old trigger interval
does not reproduce the statistics. This was missing from the initial
diagnosis; it must not be described as a non-reproducible B6h dataset.

``saved-zero`` is now the explicit default for saved-result reconstruction;
``recorded-replacement`` retains the current archived source behavior.
The latter's old failed checks remain separately recorded. No interval or
tolerance was adjusted to fit results. Numerical equivalence establishes
the matching operation, not the exact unavailable historical source version.
Raw data and saved products are unchanged.

7. Uncertainty analysis needs separate review
--------------------------------------------

The XZ statistical-error script calls normalized MATLAB autocorr, fits an
exponential plus offset, integrates the fit and then divides by velocity
variance. That extra normalization requires dimensional/methodological
review before reuse. The main 2023 mean/statistics/flux validation does
not depend on adopting this uncertainty model.

Implementation references checked
---------------------------------

MATLAB scalar-window pwelch uses Hamming windows and default half overlap:
https://www.mathworks.com/help/signal/ref/pwelch.html

Hamming's default symmetric window:
https://www.mathworks.com/help/signal/ref/hamming.html

Default skewness/kurtosis bias conventions:
https://www.mathworks.com/help/stats/skewness.html
https://www.mathworks.com/help/stats/kurtosis.html

Python uses explicit conventions, verified numerically against saved
spectra; SciPy's default Welch window/detrending is not used implicitly.

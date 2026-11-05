Expected inputs
===============

Pass the original DATA directory to ``--data-root``. It must contain:

* ``ORIGINAL/RdS_PIV_XZ_ORIGINAL/<case>/`` and
  ``ORIGINAL/RdS_PIV_YZ_ORIGINAL/<case>/``: Dantec statistical exports.
* ``POSTPROCESSING/RdS_PIV_XZ_STATISTICS/`` and
  ``POSTPROCESSING/RdS_PIV_YZ_STATISTICS/``: saved statistical profiles and contours.
* ``POSTPROCESSING/RdS_PIV_XZ_MASSFLOWRATE/``: sectional fluxes.
* ``POSTPROCESSING/RdS_HWA_XZ_STATISTICS/`` and
  ``POSTPROCESSING/RdS_HWA_XZ_PSD/``: hot-wire statistics and spectra.

The XZ profile directory also contains the AR0.7 reference profiles; the
spectra directory contains the cubic-array reference spectra.
Raw-data processing additionally requires the original instantaneous PIV
fields or hot-wire recordings, as selected by the command.

Case folders are ``1block``, ``2block``, ``3block``, ``6block`` and ``cathedral``;
their command-line identifiers are ``1b``, ``2b``, ``3b``, ``6b`` and ``ca``.
Inputs are read-only. Output paths must be outside the input data tree.

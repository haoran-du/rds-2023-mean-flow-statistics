Effect of morphology and an upstream tall building on the mean turbulence statistics of a street canyon flow
============================================================================================================

Authors: Haoran Du, Eric Savory, and Laurent Perret

Affiliations:

* Haoran Du and Eric Savory: Department of Mechanical and Materials Engineering,
  University of Western Ontario, London, Ontario, Canada.
* Laurent Perret: Nantes Universite, Ecole Centrale de Nantes, CNRS, LHEEA,
  UMR 6598, Nantes, France.

Corresponding author: Haoran Du, hdu59@uwo.ca

Published in Building and Environment (2023).
Article: https://doi.org/10.1016/j.buildenv.2023.110428

Code Repository
---------------

This repository contains Python routines for the mean-flow and turbulence
statistics of the Rue des Suisses (RdS) street canyon. The workflow reads
the PIV and hot-wire measurements and saved statistical outputs to produce
velocity profiles, contour plots, spectra and sectional air fluxes.

This is a code-only release. Experimental data, processed arrays and
generated figures are supplied or produced separately.

Repository layout
-----------------

* ``src/rds2023/analysis.py``: statistical profiles, contours, spectra and flux plots.
* ``src/rds2023/statistics.py``: velocity moments, stresses and sectional fluxes.
* ``src/rds2023/raw.py``: processing of instantaneous XZ PIV fields.
* ``src/rds2023/hwa.py``: hot-wire processing and spectral analysis.
* ``src/rds2023/io.py``: measurement-file readers and array conventions.
* ``data/README.rst``: expected inputs and directory layout.
* ``tests/``: automated tests for the Python routines.
* ``docs/``: implementation and processing notes.

Setup
-----

Use Python 3.10 or later. From the repository directory::

    python -m venv .venv

Activate the environment on Linux/macOS::

    source .venv/bin/activate

Or in Windows PowerShell::

    .\.venv\Scripts\Activate.ps1

Then install the package::

    python -m pip install -e .

Run order
---------

1. Obtain the external measurement and statistical files listed in ``data/README.rst``.
   Keep their original directory layout.
2. Reconstruct the numerical plots corresponding to Figs. 2–10::

       rds2023 reproduce --data-root "/path/to/DATA" --output outputs/paper

3. To recompute XZ statistics or hot-wire statistics from the original recordings::

       rds2023 validate-xz --data-root "/path/to/DATA" --output outputs/xz --source txt
       rds2023 validate-hwa --data-root "/path/to/DATA" --output outputs/hwa

   Add ``--cases 6b`` to either command to process B6h only.

All generated files are written to the selected output directory. ``outputs/``
is ignored by Git. ``python -m rds2023.cli`` is an alternative to ``rds2023``.

Processing scope
----------------

The default plots use the paper's reference scales. ``--normalization legacy``
selects the earlier plotting-script convention. The hot-wire workflow defaults
to the processing associated with the saved statistics and spectra.

The numerical plots use simplified layouts. The uncertainty table and manual
building outlines and optical masks are not included in this workflow.

Data and citation
-----------------

Input files are external to this repository; their expected layout is given
in ``data/README.rst``. Please cite the article linked above when using this
code. Citation metadata are also provided in ``CITATION.cff``.

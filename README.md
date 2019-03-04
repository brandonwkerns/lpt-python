# lpt-python
Python version of LPT. Geared for real time use at NCEP.

Python module dependencies:
- numpy
- scipy.signal, scipy.ndimage
- NetCDF4.Dataset

Code organization:
- Functions for reading data are in readdata.py.
- Functions for input/output are in io.py.
- Supporting functions for calculations are in helpers.py.
- The following driver scripts are included:
++ lpt_real_time_driver_tmpa.py
++ lpt_real_time_driver_cmorph.py
- Test scripts are under tests/ directory.

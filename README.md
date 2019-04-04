# lpt-python
Python version of LPT. Geared for real time use at NCEP.

*This version DOES NOT yet include splitting up LPT system groups in to track branches.*  
*This version DOES NOT yet include MJO identification.*


## Python module dependencies (see below for full environment I used):
- numpy
- scipy.signal, scipy.ndimage
- NetCDF4.Dataset


## Code organization:
- Real time data download scripts, crontab setup, and log files are in crontab/.
- The main code directory is lpt/.
  * Functions for reading data are in lpt/eaddata.py.
  * Functions for LP object and LPTs input/output are in lpt/io.py.
  * Supporting functions for calculations are in lpt/helpers.py.
  * Plotting functions are in lpt/plotting.py
  * The following driver scripts are included in lpt/:
    + lpt_real_time_driver_tmpa.py
    + lpt_real_time_driver_cmorph.py
    + lpt_real_time_driver.py  (The "master" real time driver function, called by the above two.)
- Some test scripts are under tests/ directory. Mainly used for development purposes.


## Setting up LPT on a new system:
1) Clone this repository to your system, or download the zip file format.
2) First, set the data download directories and crontab directories in the crontab/ directory.
   Set the same directory as dataset['raw_data_parent_dir'] in the driver scripts
3) Determine where the digital data output and images should go.
   Set those as output['img_dir'] and output['data_dir'] in the driver scripts.
4) (optional, but it helps to start with a 60+ day record for the LPTs tracking step)
   To start with an up to date set of plots and digital data (back to January 2019), download the following image and data directories and untar them (tar -zxvf) on your system:
     https://orca.atmos.washington.edu/~bkerns/realtime_mjo_tracking/lpt/data.tgz (128 MB)
     https://orca.atmos.washington.edu/~bkerns/realtime_mjo_tracking/lpt/images.tgz (699 MB).
5) Recommend doing a manual test run of the data download and the LPT scripts.
6) Turn on the crontab. Either add it to your current crontab, or set it as a new crontab
    (crontab lpt.cron()


-------------------------------------------------------------------------------
## In my implementation:
- The repository is in /home/orca/bkerns/lib/lpt/lpt-python
- TMPA data are downloaded to /home/orca/data/satellite/trmm_global_rainfall
- CMORPH data are downloaded to /home/orca/data/satellite/cmorph
- Digital data from LPT code is saved under


I used Anaconda Python 3.6.2 with the following environment:

```
 $ conda list                                                                                                          [12:34:44]
   packages in environment at /home/disk/atmos/bkerns/anaconda3/envs/meteo:  
backports                 1.0                      py36_1    conda-forge  
backports.functools_lru_cache 1.5                      py36_0    conda-forge  
basemap                   1.1.0                    py36_4    conda-forge  
basemap-data-hires        1.1.0                         0    conda-forge  
blas                      1.0                         mkl    
ca-certificates           2018.1.18                     0    conda-forge  
certifi                   2018.1.18                py36_0    conda-forge
cmocean                   2.0                        py_0    conda-forge
curl                      7.59.0                        0    conda-forge
cycler                    0.10.0                   py36_0    conda-forge
dbus                      1.13.0               h3a4f0e9_0    conda-forge
decorator                 4.2.1                    py36_0    conda-forge
expat                     2.2.5                         0    conda-forge
fontconfig                2.13.1               h65d0f4c_0    conda-forge
freetype                  2.9.1                h6debe1e_4    conda-forge
geos                      3.6.2                         1    conda-forge
gettext                   0.19.8.1                      0    conda-forge
glib                      2.55.0                        0    conda-forge
gst-plugins-base          1.12.5               hde13a9d_0    conda-forge
gstreamer                 1.12.5               h61a6719_0    conda-forge
h5py                      2.7.1                    py36_3    conda-forge
hdf4                      4.2.13                        0    conda-forge
hdf5                      1.10.1                        2    conda-forge
icu                       58.2                          0    conda-forge
intel-openmp              2018.0.0                      8  
ipython                   6.2.1                    py36_1    conda-forge
ipython_genutils          0.2.0                    py36_0    conda-forge
jedi                      0.11.1                   py36_0    conda-forge
jpeg                      9c                   h470a237_1    conda-forge
kiwisolver                1.0.1                    py36_1    conda-forge
krb5                      1.14.2                        0    conda-forge
libedit                   3.1.20170329                  0    conda-forge
libffi                    3.2.1                         3    conda-forge
libgcc-ng                 7.2.0                hdf63c60_3  
libgfortran               3.0.0                         1  
libgfortran-ng            7.2.0                hdf63c60_3  
libiconv                  1.15                          0    conda-forge
libnetcdf                 4.5.0                         3    conda-forge
libpng                    1.6.34                        0    conda-forge
libssh2                   1.8.0                         2    conda-forge
libstdcxx-ng              7.2.0                hdf63c60_3    conda-forge
libtiff                   4.0.9                he6b73bb_2    conda-forge
libuuid                   2.32.1               h470a237_2    conda-forge
libxcb                    1.13                          0    conda-forge
libxml2                   2.9.8                         0    conda-forge
matplotlib                2.2.3            py36h8e2386c_0    conda-forge
mkl                       2018.0.2                      1  
mkl_fft                   1.0.6                    py36_0    conda-forge
mkl_random                1.0.1                    py36_0    conda-forge
ncurses                   5.9                          10    conda-forge
netcdf4                   1.3.1                    py36_2    conda-forge
numpy                     1.13.3           py36hdbf6ddf_4  
olefile                   0.46                       py_0    conda-forge
openblas                  0.2.20                        7    conda-forge
openssl                   1.0.2p               h470a237_0    conda-forge
pandas                    0.23.4           py36hf8a1672_0    conda-forge
parso                     0.1.1                      py_0    conda-forge
pcre                      8.41                          1    conda-forge
pexpect                   4.4.0                    py36_0    conda-forge
pickleshare               0.7.4                    py36_0    conda-forge
pillow                    5.2.0            py36hc736899_1    conda-forge
pip                       9.0.1                    py36_1    conda-forge
prompt_toolkit            1.0.15                   py36_0    conda-forge
ptyprocess                0.5.2                    py36_0    conda-forge
pygments                  2.2.0                    py36_0    conda-forge
pyparsing                 2.2.0                    py36_0    conda-forge
pyproj                    1.9.5.1                  py36_0    conda-forge
pyqt                      5.6.0                    py36_4    conda-forge
pyshp                     1.2.12                     py_0    conda-forge
python                    3.6.3                h1284df2_4  
python-dateutil           2.7.2                      py_0    conda-forge
pytz                      2018.3                     py_0    conda-forge
qt                        5.6.2                hf70d934_9    conda-forge
readline                  7.0                           0    conda-forge
scipy                     1.1.0            py36hfc37229_0  
setuptools                39.0.1                   py36_0    conda-forge
simplegeneric             0.8.1                    py36_0    conda-forge
sip                       4.18                     py36_1    conda-forge
six                       1.11.0                   py36_1    conda-forge
sqlite                    3.24.0               h84994c4_0  
tk                        8.6.8                ha92aebf_0    conda-forge
tornado                   5.0.1                    py36_1    conda-forge
traitlets                 4.3.2                    py36_0    conda-forge
wcwidth                   0.1.7                    py36_0    conda-forge
wheel                     0.30.0                   py36_2    conda-forge
xarray                    0.11.0                py36_1000    conda-forge
xorg-libxau               1.0.8                         3    conda-forge
xorg-libxdmcp             1.1.2                         3    conda-forge
xz                        5.2.4                h470a237_1    conda-forge
zlib                      1.2.11                        0    conda-forge
```

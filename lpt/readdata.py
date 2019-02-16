import numpy as np
from netCDF4 import Dataset
import struct
import sys
import os

"""
This module contains functions for reading external data
to use with LPT.
"""



def read_tmpa_hdf(fn):
    """
    DATA = read_tmpa_hdf(fn)

    output:
    list(DATA)
    Out[12]: ['lon', 'lat', 'precip']
    In [21]: DATA['lon'].shape
    Out[21]: (1440,)
    In [22]: DATA['lat'].shape
    Out[22]: (400,)
    In [23]: DATA['precip'].shape
    Out[23]: (400, 1440)
    """

    ## The TMPA HDF files can be read using NetCDF4 Dataset.
    DS = Dataset(fn)
    DATA={}
    DATA['lon'] = np.arange(-179.875, 180.0, 0.25)
    DATA['lat'] = np.arange(-49.875, 50.0, 0.25)
    DATA['precip'] = DS['precipitation'][:].T
    DS.close()

    ## Need to get from (-180, 180) to (0, 360) longitude.
    lon_lt_0, = np.where(DATA['lon'] < -0.0001)
    lon_ge_0, = np.where(DATA['lon'] > -0.0001)
    DATA['lon'][lon_lt_0] += 360.0
    DATA['lon'] = np.concatenate((DATA['lon'][lon_ge_0], DATA['lon'][lon_lt_0]))
    DATA['precip'] = np.concatenate((DATA['precip'][:,lon_ge_0], DATA['precip'][:,lon_lt_0]), axis=1)

    return DATA


def read_tmpa_rt_bin(fn):
    """
    RT = read_tmpa_rt_bin(fn)

    output:
    In [24]: list(RT)
    Out[24]: ['lon', 'lat', 'precip']
    In [25]: RT['lon'].shape
    Out[25]: (1440,)
    In [26]: RT['lat'].shape
    Out[26]: (480,)
    In [27]: RT['precip'].shape
    Out[27]: (480, 1440)

    missing values (stored as -31999) are set to np.NaN.
    """

    ## TMPA RT files are binary.
    dtype=np.dtype([('field1', '<i2')])
    DATA={}
    DATA['lon'] = np.arange(0.125, 360.0, 0.25)
    DATA['lat'] = np.arange(-59.875, 60.0, 0.25)
    fid = open(fn,'rb')

    ## skip header
    fid.seek(2880)
    DATA['precip'] = np.fromfile(fid, dtype=np.int16, count=691200)
    if sys.byteorder == 'little':
        DATA['precip'] = DATA['precip'].byteswap()

    ## Shape and scale the data.
    DATA['precip'] = np.flip(np.reshape(np.double(DATA['precip']) / 100.0, [480, 1440]), axis=0)
    DATA['precip'][DATA['precip'] < -0.001] = 0.0 # Usually, missing high latitude data.
    fid.close()

    return DATA


def read_tmpa_at_datetime(dt, force_rt_tmpa=False, verbose=False):

    YYYY = dt.strftime("%Y")
    MM = dt.strftime("%m")
    DD = dt.strftime("%d")
    HH = dt.strftime("%H")
    YMD = YYYY + MM + DD

    ## First try research product
    fn = ('/home/orca/data/satellite/trmm_global_rainfall/'
       + YYYY+'/'+MM+'/'+YMD+'/3B42.'+YMD+'.'+HH+'.7.HDF')

    DATA=None
    if os.path.exists(fn) and not force_rt_tmpa:
        if verbose:
            print(fn)
        DATA=read_tmpa_hdf(fn)
    else:
        ## If no research grade, use the research product
        fn = ('/home/orca/data/satellite/trmm_global_rainfall/rt/'
           + YYYY+'/'+MM+'/'+YMD+'/3B42RT.'+YMD+HH+'.7.bin')

        ## Sometimes, the file name is "7A" instead of "7".
        if not os.path.exists(fn):
            fn = ('/home/orca/data/satellite/trmm_global_rainfall/rt/'
               + YYYY+'/'+MM+'/'+YMD+'/3B42RT.'+YMD+HH+'.7A.bin')

        if verbose:
            print(fn)
        DATA=read_tmpa_rt_bin(fn)
    return DATA

import numpy as np
from netCDF4 import Dataset
import struct
import sys

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
    DATA['precip'][DATA['precip'] < -0.001] = np.nan
    fid.close()

    return DATA

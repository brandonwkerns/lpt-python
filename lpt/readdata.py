import numpy as np
from netCDF4 import Dataset
import struct
import sys
import os

"""
This module contains functions for reading external data
to use with LPT.
"""

################################################################################
################################################################################
################################################################################

"""
CMORPH reading functions.
"""
def read_cmorph_rt_bin(fn, area=[0,360,-90,90]):

    """
    DATA = read_cmorph_rt_bin(fn)
    DATA is a dict with keys lon, lat, and precip.

    CMORPH RT files are binary.
    The GrADS control file below is used as the basis for this function:

    DSET ^../%y4/%y4%m2/CMORPH_V0.x_RT_8km-30min_%y4%m2%d2%h2
    OPTIONS little_endian template
    UNDEF -999.0
    TITLE CMORPH Rain Rate (Real-Time Version)
    XDEF  4948 LINEAR   0.0363783345 0.072756669
    YDEF  1649 LINEAR -59.963614312  0.072771376
    ZDEF     1 LEVELS   1
    TDEF 99999 LINEAR 00:00z01Jan2017 30mn
    VARS 1
    cmorph  1  99  CMORPH Rain Rate [mm/hr]
    ENDVARS
    """

    dtype=np.dtype([('field1', '<i2')])
    DATA={}
    DATA['lon'] = np.arange(0.0363783345, 360.0, 0.072756669)
    DATA['lat'] = np.arange(-59.963614312, 60.0, 0.072771376)
    fid = open(fn,'rb')

    ## GrADS uses FORTRAN REAL values, which is np.float32 for Python.
    DATA['precip'] = np.fromfile(fid, dtype=np.float32, count=2*4948*1649)
    if sys.byteorder == 'big': # Data is little endian.
        DATA['precip'] = DATA['precip'].byteswap()

    ## Shape and scale the data.
    DATA['precip'] = np.reshape(np.double(DATA['precip']), [2, 1649, 4948])
    DATA['precip'][DATA['precip'] < -0.001] = 0.0 # Usually, missing high latitude data.
    fid.close()

    ## Cut out area.
    keep_lon, = np.where(np.logical_and(DATA['lon'] > area[0], DATA['lon'] < area[1]))
    keep_lat, = np.where(np.logical_and(DATA['lat'] > area[2], DATA['lat'] < area[3]))

    DATA['lon'] = DATA['lon'][keep_lon[0]:keep_lon[-1]]
    DATA['lat'] = DATA['lat'][keep_lat[0]:keep_lat[-1]]
    DATA['precip'] = DATA['precip'][:, keep_lat[0]:keep_lat[-1], keep_lon[0]:keep_lon[-1]]

    return DATA



def read_cmorph_at_datetime(dt, force_rt=False, verbose=False, area=[0,360,-90,90]):

    """
    DATA = read_cmorph_at_datetime(dt, force_rt=False, verbose=False)

    DATA is a dict with keys lon, lat, and precip.

    Based on the provided datetime dt, read in the CMORPH data.
    By default, it will first check for the research product,
    and use the realtime product if the research product was not found.
    However, if force_rt = True, it just uses the realtime product.
    """

    YYYY = dt.strftime("%Y")
    MM = dt.strftime("%m")
    DD = dt.strftime("%d")
    HH = dt.strftime("%H")
    YMD = YYYY + MM + DD

    ## First try research product
    fn = ('/home/orca/data/satellite/cmorph/'
       + YYYY+'/'+MM+'/'+YMD+'/3B42.'+YMD+'.'+HH+'.7.HDF') #TODO: Update this for CMORPH final product.

    DATA=None
    if os.path.exists(fn) and not force_rt:
        if verbose:
            print(fn)
        DATA=read_tmpa_hdf(fn)
    else:
        ## If no research grade, use the research product
        fn = ('/home/orca/data/satellite/cmorph/rt/'
           + YYYY+'/'+MM+'/'+YMD+'/CMORPH_V0.x_RT_8km-30min_'+YMD+HH)

        if verbose:
            print(fn)
        DATA=read_cmorph_rt_bin(fn)
    return DATA


################################################################################
################################################################################
################################################################################


"""
TRMM 3B42/TMPA reading functions.
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


def read_tmpa_at_datetime(dt, force_rt=False, verbose=False):

    """
    DATA = read_tmpa_at_datetime(dt, force_rt=False, verbose=False)

    DATA is a dict with keys lon, lat, and precip.

    Based on the provided datetime dt, read in the TMPA data.
    By default, it will first check for the research product,
    and use the realtime product if the research product was not found.
    However, if force_rt = True, it just uses the realtime product.
    """

    YYYY = dt.strftime("%Y")
    MM = dt.strftime("%m")
    DD = dt.strftime("%d")
    HH = dt.strftime("%H")
    YMD = YYYY + MM + DD

    ## First try research product
    fn = ('/home/orca/data/satellite/trmm_global_rainfall/'
       + YYYY+'/'+MM+'/'+YMD+'/3B42.'+YMD+'.'+HH+'.7.HDF')

    DATA=None

    ## Sometimes, the file name is "7A" instead of "7".
    if not os.path.exists(fn):
        fn = ('/home/orca/data/satellite/trmm_global_rainfall/'
           + YYYY+'/'+MM+'/'+YMD+'/3B42.'+YMD+'.'+HH+'.7A.HDF')

    if os.path.exists(fn) and not force_rt:
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

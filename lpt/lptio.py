import matplotlib; matplotlib.use('agg')
import numpy as np


###################################################
### Output functions
###################################################
def lp_objects_output_ascii(fn, OBJ):
    """
    This function outputs the "bulk" LP object properties (centroid, date, area)
    to an ascii file.
    """
    print('Writing ascii output to: '+fn)
    fmt = '%7.2f%8.2f%7.1f%7.1f%20.1f   %16d\n'
    file = open(fn, 'w')

    for ii in range(len(OBJ['lon'])):

        print(fmt % (OBJ['lat'][ii], OBJ['lon'][ii],
                OBJ['y'][ii], OBJ['x'][ii],
                OBJ['area'][ii], OBJ['id'][ii]))

        file.write(fmt % (OBJ['lat'][ii], OBJ['lon'][ii],
                OBJ['y'][ii], OBJ['x'][ii],
                OBJ['area'][ii], OBJ['id'][ii]))

    file.close()


def lp_objects_output_netcdf(fn, OBJ, label_im):
    pass
    """
    This function outputs the "bulk" LP object properties (centroid, date, area)
    Plus the pixel information to a compressed netcdf file.
    """

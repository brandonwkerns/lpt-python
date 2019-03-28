import matplotlib; matplotlib.use('agg')
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset

###################################################
### Output functions
###################################################
def lp_objects_output_ascii(fn, OBJ):
    """
    This function outputs the "bulk" LP object properties (centroid, date, area)
    to an ascii file.
    """
    print('Writing ascii output to: ' + fn)
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


def lp_objects_output_netcdf(fn, OBJ):
    """
    This function outputs the "bulk" LP object properties (centroid, date, area)
    Plus the pixel information to a compressed netcdf file.
    """
    print('Writing netcdf output to: ' + fn)

    DS = Dataset(fn, 'w', format='NETCDF4_CLASSIC', clobber=True)
    DS.description = 'LP Objects NetCDF file.'

    ##
    ## Dimensions
    ##
    DS.createDimension('nobj', len(OBJ['lon']))
    max_points = np.max(OBJ['n_points'])
    DS.createDimension('npoints', max_points)

    ## Grid stuff.
    DS.createDimension('grid_x', len(OBJ['grid']['lon']))
    DS.createDimension('grid_y', len(OBJ['grid']['lat']))

    ##
    ## Variables
    ##
    var_centroid_lon = DS.createVariable('centroid_lon','f4',('nobj',))
    var_centroid_lat = DS.createVariable('centroid_lat','f4',('nobj',))
    var_area = DS.createVariable('area','f4',('nobj',))
    var_objid = DS.createVariable('objid','d',('nobj',))

    var_pixels_x = DS.createVariable('pixels_x','i4',('nobj','npoints',), zlib=True)
    var_pixels_y = DS.createVariable('pixels_y','i4',('nobj','npoints',), zlib=True)

    var_grid_lon = DS.createVariable('grid_lon','f4',('grid_x',))
    var_grid_lat = DS.createVariable('grid_lat','f4',('grid_y',))
    var_grid_area = DS.createVariable('grid_area','f4',('grid_y','grid_x',), zlib=True)
    var_grid_mask = DS.createVariable('grid_mask','i4',('grid_y','grid_x',), zlib=True)

    ##
    ## Values
    ##
    var_centroid_lon[:] = OBJ['lon']
    var_centroid_lat[:] = OBJ['lat']
    var_area[:] = OBJ['area']
    var_objid[:] = OBJ['id']

    var_grid_lon[:] = OBJ['grid']['lon']
    var_grid_lat[:] = OBJ['grid']['lat']
    var_grid_area[:] = OBJ['grid']['area']
    mask = OBJ['label_im'] - 1
    mask = ma.masked_array(mask, mask = (mask < -0.5))
    var_grid_mask[:] = mask

    for ii in range(len(OBJ['lon'])):
        ypoints, xpoints = np.where(OBJ['label_im'] == ii+1)
        var_pixels_x[ii,:len(xpoints)] = xpoints
        var_pixels_y[ii,:len(ypoints)] = ypoints

    ##
    ## Attributes/Metadata
    ##
    var_centroid_lon.setncatts({'units':'degrees_east','long_name':'centroid longitude (0-360)','standard_name':'longitude','axis':'X'})
    var_centroid_lat.setncatts({'units':'degrees_north','long_name':'centroid latitude (-90-90)','standard_name':'latitude','axis':'Y'})
    var_area.setncatts({'units':'km2','long_name':'LP object enclosed area'})
    var_objid.setncatts({'units':'0','long_name':'LP Object ID'
        ,'description':'A unique ID for each LP object. Convention is YYYYMMDDHHnnnn where nnnn starts at 0000'})

    var_grid_lon.setncatts({'units':'degrees_east','long_name':'grid longitude (0-360)'})
    var_grid_lat.setncatts({'units':'degrees_north','long_name':'grid latitude (-90-90)'})
    var_grid_area.setncatts({'units':'km2','long_name':'area of each grid point'})
    var_grid_area.setncatts({'units':'0','long_name':'mask by nnnn part of LP Object ID'})

    var_pixels_x.setncatts({'units':'0','long_name':'grid point pixel indices in the x direction','note':'zero based'})
    var_pixels_y.setncatts({'units':'0','long_name':'grid point pixel indices in the y direction','note':'zero based'})

    DS.close()

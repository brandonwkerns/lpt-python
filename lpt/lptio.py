import matplotlib; matplotlib.use('agg')
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import os
import os.path

###################################################
### Output functions
###################################################
def lp_objects_output_ascii(fn, OBJ):
    """
    This function outputs the "bulk" LP object properties (centroid, date, area)
    to an ascii file.
    """
    print('Writing LP object ASCII output to: ' + fn)
    fmt = '%7.2f%8.2f%7.1f%7.1f%20.1f   %16d\n'
    file = open(fn, 'w')

    file.write(' lat.__  lon.__    y._    x._         area[km2]._     YYYYMMDDHHnnnn\n') # Header line.
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
    print('Writing LP object NetCDF output to: ' + fn)

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


def lpt_system_tracks_output_ascii(fn, TIMECLUSTERS):
    """
    This function outputs the "bulk" LPT system properties (centroid, date, area)
    to an ascii file.
    """
    print('Writing LPT system track ASCII output to: ' + fn)
    fmt='        %4d%02d%02d%02d %8d %10.2f %10.2f %2d\n'

    os.makedirs(os.path.dirname(fn), exist_ok=True) # Make directory if needed.
    file = open(fn, 'w')

    ## Header
    file.write("LPT nnnn.nn\n")
    file.write("        YYYYMMDDHH _A_[km2] cen_lat.__ cen_lon.__ Nobj\n")

    ## Data
    for ii in range(len(TIMECLUSTERS)):
        file.write("LPT %07.2f\n" % (TIMECLUSTERS[ii]['lpt_id'],))

        for tt in range(len(TIMECLUSTERS[ii]['datetime'])):
            year,month,day,hour = TIMECLUSTERS[ii]['datetime'][tt].timetuple()[0:4]
            file.write(fmt % (year,month,day,hour
                                , TIMECLUSTERS[ii]['area'][tt]
                                , TIMECLUSTERS[ii]['centroid_lat'][tt]
                                , TIMECLUSTERS[ii]['centroid_lon'][tt]
                                , TIMECLUSTERS[ii]['nobj'][tt]))

    file.close()


def lpt_system_tracks_output_netcdf(fn, TIMECLUSTERS):
    """
    This function outputs the "bulk" LPT system properties (centroid, date, area)
    plus the LP Objects belonging to each "TIMECLUSTER" to an ascii file.
    """
    print('Writing LPT system track NetCDF output to: ' + fn)

    os.makedirs(os.path.dirname(fn), exist_ok=True) # Make directory if needed.

    DS = Dataset(fn, 'w', format='NETCDF4_CLASSIC', clobber=True)
    DS.description = 'LPT Systems "timeclusters" NetCDF file.'

    MISSING = np.nan
    FILL_VALUE = -990.0
    ##
    ## Dimensions
    ##
    DS.createDimension('nlpt', len(TIMECLUSTERS))
    max_points = 1
    max_times = 1
    timestamp_collect = np.array([MISSING])
    centroid_lon_collect = np.array([MISSING])
    max_lon_collect = np.array([MISSING])
    centroid_lat_collect = np.array([MISSING])
    area_collect = np.array([MISSING])
    for ii in range(len(TIMECLUSTERS)):
        max_points = max(max_points, len(TIMECLUSTERS[ii]['objid']))
        max_times = max(max_times, len(TIMECLUSTERS[ii]['datetime']))
        timestamp_collect = np.append(np.append(timestamp_collect, TIMECLUSTERS[ii]['timestamp']/3600.0),MISSING)
        centroid_lon_collect = np.append(np.append(centroid_lon_collect, TIMECLUSTERS[ii]['centroid_lon']),MISSING)
        max_lon_collect = np.append(np.append(max_lon_collect, TIMECLUSTERS[ii]['max_lon']),MISSING)
        centroid_lat_collect = np.append(np.append(centroid_lat_collect, TIMECLUSTERS[ii]['centroid_lat']),MISSING)
        area_collect = np.append(np.append(area_collect, TIMECLUSTERS[ii]['area']),MISSING)

    DS.createDimension('nobj', max_points)
    DS.createDimension('ntimes', max_times)
    DS.createDimension('nall', len(timestamp_collect))

    ##
    ## Variables
    ##

    ## Stitchec "bulk" variables.
    var_centroid_lon_all = DS.createVariable('centroid_lon_stitched','f4',('nall',))
    #var_max_lon_all = DS.createVariable('max_lon_stitched','f4',('nall',))
    var_centroid_lat_all = DS.createVariable('centroid_lat_stitched','f4',('nall',))
    var_area_all = DS.createVariable('area_stitched','d',('nall',))
    var_timestamp_all = DS.createVariable('timestamp_stitched','d',('nall',))

    ##
    ## Values
    ##
    var_centroid_lon_all[:] = centroid_lon_collect
    #var_max_lon_all[:] = max_lon_collect
    var_centroid_lat_all[:] = centroid_lat_collect
    var_area_all[:] = area_collect
    var_timestamp_all[:] = timestamp_collect

    ##
    ## Attributes/Metadata
    ##
    var_centroid_lon_all.setncatts({'units':'degrees_east','long_name':'centroid longitude (0-360) -- stitched','standard_name':'longitude'})
    #var_max_lon_all.setncatts({'units':'degrees_east','long_name':'eastmost longitude (0-360) -- stitched','standard_name':'longitude'})
    var_centroid_lat_all.setncatts({'units':'degrees_north','long_name':'centroid latitude (-90-90) -- stitched','standard_name':'latitude'})
    var_area_all.setncatts({'units':'km2','long_name':'LPT System enclosed area -- stitched'})
    var_timestamp_all.setncatts({'units':'hours since 1970-1-1 0:0','long_name':'LPT System time stamp -- stitched'})

    DS.close()

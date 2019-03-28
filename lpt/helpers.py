import matplotlib; matplotlib.use('agg')
import numpy as np
from scipy.signal import convolve2d
from scipy import ndimage
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

## These functions are used for LPT.

############################################################################
######################  Gaussian Smoothing Functions  ######################
############################################################################


def calc_scaled_average(data_in_accumulation_period, factor):
    """
    accumulated_data = calc_accumulation(data_in accumulation_period, factor)

    Calculate the sum and multiply by the data time interval to get the accumulation.
    -- data_in_accumulation_period[t,y,x] is a 3D array.
    -- factor gets multiplied by the mean. E.g., if the data is rain rate in mm/h,
       using factor of 24 would be in mm/day.
    """

    return factor * np.nanmean(data_in_accumulation_period, axis=0)


def identify_lp_objects(field, threshold
                        , object_is_gt_threshold=True, verbose=False):

    """
    label_im = identify_lp_objects(lon, lat, field, threshold
                            , object_minimum_gridpoints=0
                            , object_is_gt_threshold=True)

    Given an input data field (e.g., already accumulated and filtered),
    identify the LP Objects in that field. Return an array the same size
    as field, but with values indexed by object IDs.
    """

    field_bw = 0 * field
    if object_is_gt_threshold:
        field_bw[(field > threshold)] = 1
    else:
        field_bw[(field < threshold)] = 1


    label_im, nb_labels = ndimage.label(field_bw)
    if verbose:
        print('Found '+str(nb_labels)+' objects.', flush=True) # how many regions?

    return label_im


def calc_grid_cell_area(lon, lat):

    """
    area = calc_grid_cell_area(lon, lat)

    Given lon and lat arrays, calculate the area of each grid cell.
    - lon and lat don't need to be a uniform grid, but they need to be increasing
      in both the x and y direction for this function to work.
    - If 1-D arrays are given, they will be converted to 2D using np.meshgrid.
    """

    area = None
    if lon.ndim == 1:
        print('ERROR: lon and lat must be 2D arrays for function calc_grid_cell_area.', flush=True)
    else:
        ny,nx = lon.shape
        dlon = 0.0*lon
        dlat = 0.0*lat

        dlon[:,1:nx-1] = 0.5*(lon[:,1:nx-1] + lon[:,2:nx]) - 0.5*(lon[:,0:nx-2] + lon[:,1:nx-1])
        dlon[:,0] = dlon[:,1]
        dlon[:,nx-1] = dlon[:,nx-2]
        dlat[1:ny-1,:] = 0.5*(lat[1:ny-1,:] + lat[2:ny,:]) - 0.5*(lat[0:ny-2,:] + lat[1:ny-1,:])
        dlat[0,:] = dlat[1,:]
        dlat[ny-1,:] = dlat[ny-2,:]

        area = (dlat*111.195) * (dlon*111.195*np.cos(np.pi*lat/180.0))

    return area



def calculate_lp_object_properties(lon, lat, field, field_accum, label_im
                        , object_minimum_gridpoints, end_of_accumulation_time
                        , verbose=False):

    nb_labels = np.max(label_im)
    mask = 1*label_im
    mask[label_im > 0] = 1

    ## If lon and lat not in 2d arrays, put them through np.meshgrid.
    if lon.ndim == 1:
        if verbose:
            print('Detected 1-D lat/lon. Using np.meshgrid to get 2d lat/lon.', flush=True)
        lon2, lat2 = np.meshgrid(lon, lat)
    else:
        lon2 = lon
        lat2 = lat

    X2, Y2 = np.meshgrid(np.arange(lon2.shape[1]), np.arange(lon2.shape[0]))

    area2d = calc_grid_cell_area(lon2, lat2)

    sizes = ndimage.sum(mask, label_im, range(1, nb_labels + 1))
    mean_instantaneous_rain_rate = ndimage.mean(field, label_im, range(1, nb_labels + 1))
    mean_accum_rain_rate = ndimage.mean(field_accum, label_im, range(1, nb_labels + 1))
    centroid_lon = ndimage.mean(lon2, label_im, range(1, nb_labels + 1))
    centroid_lat = ndimage.mean(lat2, label_im, range(1, nb_labels + 1))
    centroid_x = ndimage.mean(X2, label_im, range(1, nb_labels + 1))
    centroid_y = ndimage.mean(Y2, label_im, range(1, nb_labels + 1))
    area = ndimage.sum(area2d, label_im, range(1, nb_labels + 1))


    ## Assign LPT IDs. Order is by longitude. Use zero-base indexing.
    id0 = 1e10 * end_of_accumulation_time.year + 1e8 * end_of_accumulation_time.month + 1e6 * end_of_accumulation_time.day + 1e4 * end_of_accumulation_time.hour

    id = id0 + np.arange(len(centroid_lon))

    ## Prepare output dict.
    OBJ={}
    OBJ['id'] = id
    OBJ['label_im'] = label_im
    OBJ['lon'] = centroid_lon
    OBJ['lat'] = centroid_lat
    OBJ['x'] = centroid_x
    OBJ['y'] = centroid_y
    OBJ['n_points'] = sizes
    OBJ['area'] = area
    OBJ['mean_inst'] = mean_instantaneous_rain_rate
    OBJ['mean_accum'] = mean_accum_rain_rate

    # Grid stuff.
    OBJ['grid'] = {}
    OBJ['grid']['lon'] = lon
    OBJ['grid']['lat'] = lat
    OBJ['grid']['area'] = area2d

    return OBJ



###################################################
### Plotting functions. ###########################
###################################################
def plot_map_background(plotArea=[0,360,-60,60], lon_labels = [1,0,0,0], lat_labels = [0,0,0,1], res='c', anchor='C'
                        , coast_color = 'k', fontsize=10):
    map = Basemap(projection='cyl',resolution=res, anchor=anchor,
                  llcrnrlat = plotArea[2],
                  urcrnrlat = plotArea[3],
                  llcrnrlon = plotArea[0],
                  urcrnrlon = plotArea[1]);

    #map.fillcontinents(color='lightgray', lake_color='lightgray')
    map.drawcoastlines(linewidth=0.5, color=coast_color) ;
    map.drawparallels(np.arange(-80,81,20), linewidth=0.5, labels = lon_labels, fontsize = fontsize);
    map.drawmeridians(np.arange(0,361,20), linewidth=0.5, labels = lat_labels, fontsize = fontsize);

    return map ;

def print_and_save(file_out_base):
    print(file_out_base + '.png')
    plt.savefig(file_out_base + '.png' ,bbox_inches='tight', dpi=150)

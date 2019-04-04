import matplotlib; matplotlib.use('agg')
import numpy as np
import datetime as dt
from scipy.signal import convolve2d
from scipy import ndimage
import matplotlib.pylab as plt
from netCDF4 import Dataset
import glob

## These functions are used for LPT.

###################################################################
######################  LP Object Functions  ######################
###################################################################


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


def get_objid_datetime(objid):
    """
    usge: this_datetime = get_objid_datetime(this_objid)

    Get the datetime from an objid of form YYYYMMDDHHnnnn.
    """
    ymdh_int = int(np.floor(objid/1e4))
    ymdh_str = str(ymdh_int)
    return dt.datetime.strptime(ymdh_str, "%Y%m%d%H")


def read_lp_object_properties(objid, objdir, property_list, verbose=False):

    dt1 = get_objid_datetime(objid)
    fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")
    fn1 = (objdir + dt1.strftime(fmt))

    if verbose:
        print(fn1)

    DS1 = Dataset(fn1)
    id1 = DS1['objid'][:]
    idx1, = np.where(np.abs(id1 - objid) < 0.1)

    out_dict = {}
    for property in property_list:
        out_dict[property] = to1d(DS1[property][:][idx1])
    DS1.close()

    return out_dict


def get_latest_lp_object_time(objdir):
    obj_file_list = sorted(glob.glob((objdir + "/????/??/????????/*.nc")))
    last_obj_file = obj_file_list[-1]
    return dt.datetime.strptime(last_obj_file[-13:-3], "%Y%m%d%H")


##################################################################
######################  Tracking Functions  ######################
##################################################################

def to1d(ndarray_or_ma):
    try:
        fout = ndarray_or_ma.compressed()
    except:
        fout = ndarray_or_ma.flatten()
    return fout


def calc_overlapping_points(objid1, objid2, objdir):

    dt1 = get_objid_datetime(objid1)
    dt2 = get_objid_datetime(objid2)

    fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")
    fn1 = (objdir + dt1.strftime(fmt))
    fn2 = (objdir + dt2.strftime(fmt))

    DS1 = Dataset(fn1)
    id1 = DS1['objid'][:]
    idx1, = np.where(np.abs(id1 - objid1) < 0.1)
    x1 = to1d(DS1['pixels_x'][:][idx1])
    y1 = to1d(DS1['pixels_y'][:][idx1])
    DS1.close()

    DS2 = Dataset(fn2)
    id2 = DS2['objid'][:]
    idx2, = np.where(np.abs(id2 - objid2) < 0.1)
    x2 = to1d(DS2['pixels_x'][:][idx2])
    y2 = to1d(DS2['pixels_y'][:][idx2])
    DS2.close()

    xy1 = set(zip(x1,y1))
    xy2 = set(zip(x2,y2))
    overlap = [x in xy2 for x in xy1]

    return (len(x1), len(x2), np.sum(overlap))


def init_lpt_group_array(dt_list, objdir):
    """
    "LPT" is a 2-D array with columns: [timestamp, objid, lpt_group_id, begin_point, end_point, split_point]
    -- timestamp = Linux time stamp (e.g., seconds since 00 UTC 1970-1-1)
    -- objid = LP object id (YYYYMMDDHHnnnn)
    -- lpt_group_id = LPT group id, connected LP objects have a common LPT group id.
    -- begin point = 1 if it is the beginning of a track. 0 otherwise.
    -- end point = 1 if no tracks were connected to it, 0 otherwise.
    -- split point = 1 if split detected, 0 otherwise.
    """

    fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")

    LPT = []

    for this_dt in dt_list:

        fn = (objdir + this_dt.strftime(fmt))
        DS = Dataset(fn)
        id_list = DS['objid'][:]
        DS.close()

        for ii in range(len(id_list)):
            LPT.append([this_dt.timestamp(), id_list[ii], -1, 0, 1, 0])

    return np.array(LPT)


def lpt_group_array_remove_small_objects(LPT, options, verbose=False):
    """
    LPT comes from the init_lpt_group_array function
    options needs:
    options['min_lp_objects_points']
    """
    objdir = options['objdir']
    keep_list = np.full(len(LPT[:,1]), True, dtype=bool)

    for ii in range(len(LPT[:,1])):
        this_objid = LPT[ii,1]

        dt1 = get_objid_datetime(this_objid)
        fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")
        fn1 = (objdir + dt1.strftime(fmt))

        if verbose:
            print(fn1)

        DS1 = Dataset(fn1)
        id1 = DS1['objid'][:]
        idx1, = np.where(np.abs(id1 - this_objid) < 0.1)
        x1 = to1d(DS1['pixels_x'][:][idx1])
        DS1.close()

        if (len(x1) < options['min_lp_objects_points']):
            keep_list[ii] = False

    return LPT[keep_list,:]


def calc_lpt_group_array(LPT, options, verbose=False, reversed=False):

    """
    usage: LPT = calc_lpt_group_array(LPT, objdir, options)
    Calculate the simple LPT groups.

    options dictionary entries needed:
    options['objdir']
    options['min_overlap_points']
    options['min_overlap_frac']

    "LPT" is a 2-D array with columns: [timestamp, objid, lpt_group_id, begin_point, end_point, split_point]
    -- timestamp = Linux time stamp (e.g., seconds since 00 UTC 1970-1-1)
    -- objid = LP object id (YYYYMMDDHHnnnn)
    -- lpt_group_id = LPT group id, connected LP objects have a common LPT group id.
    -- begin point = 1 if it is the beginning of a track. 0 otherwise.
    -- end point = 1 if no tracks were connected to it, 0 otherwise.
    -- split point = 1 if split detected, 0 otherwise.
    """

    objdir = options['objdir']
    next_lpt_group_id = 0.0

    time_list = np.unique(LPT[:,0])
    if reversed:
        time_list = time_list[::-1]
        LPT = LPT[::-1,:]

    first_time = time_list[0]
    first_time_idx, = np.where(np.abs(LPT[:,0] - first_time) < 0.1)

    for ii in range(len(first_time_idx)):
        LPT[first_time_idx[ii],2] = next_lpt_group_id
        LPT[first_time_idx[ii],3] = 1
        next_lpt_group_id += 1



    for tt in range(1,len(time_list)):
        this_time = time_list[tt]
        prev_time = time_list[tt-1]

        if verbose:
            print(dt.datetime.fromtimestamp(this_time), flush=True)

        this_time_idx, = np.where(np.abs(LPT[:,0] - this_time) < 0.1)
        prev_time_idx, = np.where(np.abs(LPT[:,0] - prev_time) < 0.1)

        already_connected_objid_list = [] # Keep track of previously connected objids, for split detection.

        for ii in this_time_idx:
            this_objid = LPT[ii,1]
            match = -1
            for jj in prev_time_idx:
                prev_objid = LPT[jj,1]

                n_this, n_prev, n_overlap = calc_overlapping_points(this_objid,prev_objid,objdir)

                if n_overlap >= options['min_overlap_points']:
                    match=jj
                if 1.0*n_overlap/n_this > options['min_overlap_frac']:
                    match=jj
                if 1.0*n_overlap/n_prev > options['min_overlap_frac']:
                    match=jj

            if match>-1:
                LPT[ii,2] = LPT[match,2] # assign the same group.
                LPT[match,4] = 0 # this one got matched, so not a track end point.
                if LPT[match,1] in already_connected_objid_list:
                    LPT[match,5] = 1
                else:
                    already_connected_objid_list.append(LPT[match,1])
            else:
                LPT[ii,2] = next_lpt_group_id
                LPT[ii,3] = 1 # This is the beginning of a track.
                next_lpt_group_id += 1

    if reversed:
        LPT = LPT[::-1,:]

    return LPT

def overlap_forward_backward(LPT1, LPT2, options, verbose=True):
    """
    LPT1 is meant to be forward.
    LPT2 is meant to be backwards.
    Begin and end points are from LPT1 (forward).
    """

    LPT3 = LPT1.copy()
    unique_lpt_groups2 = np.unique(LPT2[:,2]) # Does not change.


    more_to_do = True
    while more_to_do:
        more_to_do = False

        unique_lpt_groups3 = np.unique(LPT3[:,2]) # Changes each iteration.


        for this_lpt_group3 in unique_lpt_groups3:
            print((str(this_lpt_group3) + ' of ' + str(np.max(unique_lpt_groups3)-1)), flush=True)
            idx_this_lpt_group3 = np.where(LPT3[:,2] == this_lpt_group3)[0]
            objid_in_this_lpt_group3 = LPT3[idx_this_lpt_group3 , 1]

            for this_lpt_group2 in unique_lpt_groups2:
                idx_this_lpt_group2 = np.where(LPT2[:,2] == this_lpt_group2)[0]
                objid_in_this_lpt_group2 = LPT2[idx_this_lpt_group2 , 1]

                if (len(np.intersect1d(objid_in_this_lpt_group3, objid_in_this_lpt_group2)) > 0
                        and len(np.setdiff1d(objid_in_this_lpt_group2, objid_in_this_lpt_group3)) > 0):

                    ## Add the ones from overlapping reversed track
                    LPT3[idx_this_lpt_group2,2] = this_lpt_group3


                    unique_lpt_groups1 = np.unique(LPT1[idx_this_lpt_group2,2])
                    for this_lpt_group1 in unique_lpt_groups1:
                        idx_this_lpt_group1 = np.where(LPT1[:,2] == this_lpt_group1)[0]
                        LPT3[idx_this_lpt_group1,2] = this_lpt_group3

                    #unique_lpt_groups2 = np.delete(unique_lpt_groups2, np.where(unique_lpt_groups2==this_lpt_group2)[0])
                    print((' ---> '+ str(this_lpt_group2) + ' of ' + str(np.max(unique_lpt_groups2)-1)), flush=True)
                    print('New Iteration!', flush=True)
                    more_to_do = True
                    break
            if more_to_do:
                break

    LPT3[:,3] = LPT1[:,3] * LPT2[:,4]
    LPT3[:,4] = LPT1[:,4] * LPT2[:,3]
    LPT3[:,5] = LPT1[:,5] + LPT2[:,5]
    return reorder_LPT_group_id(LPT3)


def lpt_group_array_allow_center_jumps(LPT, center_jump_max_hour):
    """
    Check duration of "end" to "start" points, and connect if less than
    center_jump_max_hour.
    """
    pass


def reorder_LPT_group_id(LPT):

    """
    re-order and relabel the LPT group IDs to 0 to N
    where N is the number of unique LPT group IDs.
    """

    LPT2 = LPT.copy()
    ## Re-order LPT system groups
    unique_lpt_groups = np.unique(LPT[:,2])

    for jjj in range(len(unique_lpt_groups)):
        LPT2[(LPT2[:,2] == unique_lpt_groups[jjj]), 2] = jjj

    return LPT2



def calc_lpt_system_group_properties(LPT, options):

    unique_lpt_groups = np.unique(LPT[:,2])

    TC_all = []

    for this_group in unique_lpt_groups:
        TC_this = {}
        TC_this['lpt_group_id'] = this_group
        TC_this['lpt_id'] = this_group

        this_lpt_group_idx = np.where(LPT[:,2] == this_group)[0]
        TC_this['objid'] = LPT[this_lpt_group_idx,1]
        TC_this['timestamp'] = np.unique(LPT[this_lpt_group_idx,0])
        TC_this['datetime'] = [dt.datetime.fromtimestamp(x) for x in TC_this['timestamp']]

        ##
        ## Sum/average the LPTs to get bulk/mean properties at each time.
        ##

        ## Initialize
        TC_this['nobj'] = np.zeros(len(TC_this['timestamp']))
        TC_this['area'] = np.zeros(len(TC_this['timestamp']))
        TC_this['centroid_lon'] = np.zeros(len(TC_this['timestamp']))
        TC_this['centroid_lat'] = np.zeros(len(TC_this['timestamp']))
        TC_this['min_lon'] =  999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['max_lon'] = -999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['min_lat'] =  999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['max_lat'] = -999.0 * np.ones(len(TC_this['timestamp']))

        ## Loop over time.
        for tt in range(len(TC_this['timestamp'])):
            idx_for_this_time = np.where(np.logical_and(
                LPT[:,0] == TC_this['timestamp'][tt],
                LPT[:,2] == this_group))[0]

            for this_objid in LPT[idx_for_this_time,1]:

                OBJ = read_lp_object_properties(this_objid, options['objdir'], ['centroid_lon','centroid_lat','area','pixels_x','pixels_y'])

                TC_this['nobj'][tt] += 1
                TC_this['area'][tt] += OBJ['area']
                TC_this['centroid_lon'][tt] += OBJ['centroid_lon'] * OBJ['area']
                TC_this['centroid_lat'][tt] += OBJ['centroid_lat'] * OBJ['area']

            TC_this['centroid_lon'][tt] /= TC_this['area'][tt]
            TC_this['centroid_lat'][tt] /= TC_this['area'][tt]

        TC_all.append(TC_this)

    return TC_all


def separate_lpt_system_branches(LPTfb, LPTf, LPTb, options):
    LPT_with_branches = LPTfb.copy() # Start with forward/backwards merged system.


def calc_lpt_system_group_properties_with_branches(LPT_with_branches, options):
    pass



###################################################
### Other processing functions. ###########################
###################################################


def get_lpo_mask(objid, objdir):

    dt1 = get_objid_datetime(objid)

    fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")
    fn1 = (objdir + dt1.strftime(fmt))

    DS1 = Dataset(fn1)
    id1 = DS1['objid'][:]
    idx1, = np.where(np.abs(id1 - objid) < 0.1)

    x1 = DS1['pixels_x'][:][idx1].compressed()
    y1 = DS1['pixels_y'][:][idx1].compressed()
    lon = DS1['grid_lon'][:]
    lat = DS1['grid_lat'][:]

    DS1.close()

    mask = np.zeros([len(lat), len(lon)])
    mask[y1,x1] = 1

    return (lon, lat, mask)

def plot_lpt_groups_time_lon_text(ax, LPT, options, text_color='k'):

    objdir = options['objdir']
    dt_min = dt.datetime.fromtimestamp(np.min(LPT[:,0]))
    dt_max = dt.datetime.fromtimestamp(np.max(LPT[:,0]))

    for ii in range(len(LPT[:,0])):
        objid = LPT[ii,1]
        dt1 = get_objid_datetime(objid)
        fmt = ("/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc")
        fn1 = (objdir + dt1.strftime(fmt))

        DS1 = Dataset(fn1)
        id1 = DS1['objid'][:]
        idx1, = np.where(np.abs(id1 - objid) < 0.1)
        lon = DS1['centroid_lon'][:][idx1]
        DS1.close()

        this_text_color = text_color
        this_zorder = 10
        if (LPT[ii,3] == 1):
            this_text_color = 'b'
            this_zorder = 20
        if (LPT[ii,4] == 1):
            this_text_color = 'm'
            this_zorder = 20
        if (LPT[ii,5] == 1):
            this_text_color = 'g'
            this_zorder = 20

        plt.text(lon, dt1, str(int(LPT[ii,2])), color=this_text_color, zorder=this_zorder)

    ax.set_xlim([0.0, 360.0])
    ax.set_ylim([dt_min, dt_max])


def plot_timeclusters_time_lon(ax, TIMECLUSTERS, linewidth=2.0):

    for ii in range(len(TIMECLUSTERS)):
        x = TIMECLUSTERS[ii]['centroid_lon']
        y = TIMECLUSTERS[ii]['datetime']
        ax.plot(x, y, 'k', linewidth=linewidth)

        plt.text(x[0], y[0], str(int(ii)), fontweight='bold', color='red')
        plt.text(x[-1], y[-1], str(int(ii)), fontweight='bold', color='red')

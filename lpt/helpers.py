import matplotlib; matplotlib.use('agg')
import numpy as np
import datetime as dt
from scipy.signal import convolve2d
from scipy import ndimage
import matplotlib.pylab as plt
from netCDF4 import Dataset
import glob


## These functions are used for LPT.

"""
def get_col(list2d, col_num):
    return [list2d[x][col_num] for x in range(len(list2d))]

def get_rows(list2d, rows_list):
    return [list2d[x] for x in rows_list]

def get_row_indices_of_column_value()
[x for x in range(len(LPT3)) if LPT3[x][2] == this_lpt_group3]
"""

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

        dlon[:,1:nx-1] = abs(0.5*(lon[:,1:nx-1] + lon[:,2:nx]) - 0.5*(lon[:,0:nx-2] + lon[:,1:nx-1]))
        dlon[:,0] = dlon[:,1]
        dlon[:,nx-1] = dlon[:,nx-2]
        dlat[1:ny-1,:] = abs(0.5*(lat[1:ny-1,:] + lat[2:ny,:]) - 0.5*(lat[0:ny-2,:] + lat[1:ny-1,:]))
        dlat[0,:] = dlat[1,:]
        dlat[ny-1,:] = dlat[ny-2,:]

        area = (dlat*111.195) * (dlon*111.195*np.cos(np.pi*lat/180.0))

    return area



def calculate_lp_object_properties(lon, lat, field, field_running, field_filtered, label_im
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

    amean_instantaneous_field = ndimage.mean(field, label_im, range(1, nb_labels + 1))
    amean_running_field = ndimage.mean(field_running, label_im, range(1, nb_labels + 1))
    amean_filtered_running_field = ndimage.mean(field_filtered, label_im, range(1, nb_labels + 1))
    min_instantaneous_field = ndimage.minimum(field, label_im, range(1, nb_labels + 1))
    min_running_field = ndimage.minimum(field_running, label_im, range(1, nb_labels + 1))
    min_filtered_running_field = ndimage.minimum(field_filtered, label_im, range(1, nb_labels + 1))
    max_instantaneous_field = ndimage.maximum(field, label_im, range(1, nb_labels + 1))
    max_running_field = ndimage.maximum(field_running, label_im, range(1, nb_labels + 1))
    max_filtered_running_field = ndimage.maximum(field_filtered, label_im, range(1, nb_labels + 1))

    centroid_lon = ndimage.mean(lon2, label_im, range(1, nb_labels + 1))
    centroid_lat = ndimage.mean(lat2, label_im, range(1, nb_labels + 1))
    centroid_x = ndimage.mean(X2, label_im, range(1, nb_labels + 1))
    centroid_y = ndimage.mean(Y2, label_im, range(1, nb_labels + 1))
    area = ndimage.sum(area2d, label_im, range(1, nb_labels + 1))
    max_lon = ndimage.maximum(lon2, label_im, range(1, nb_labels + 1))
    min_lon = ndimage.minimum(lon2, label_im, range(1, nb_labels + 1))
    max_lat = ndimage.maximum(lat2, label_im, range(1, nb_labels + 1))
    min_lat = ndimage.minimum(lat2, label_im, range(1, nb_labels + 1))


    ## Assign LPT IDs. Order is by longitude. Use zero-base indexing.
    id0 = 1e10 * end_of_accumulation_time.year + 1e8 * end_of_accumulation_time.month + 1e6 * end_of_accumulation_time.day + 1e4 * end_of_accumulation_time.hour

    id = id0 + np.arange(len(centroid_lon))

    ## Prepare output dict.
    OBJ={}
    OBJ['id'] = id
    OBJ['label_im'] = label_im
    OBJ['lon'] = centroid_lon
    OBJ['lat'] = centroid_lat
    OBJ['min_lon'] = min_lon
    OBJ['min_lat'] = min_lat
    OBJ['max_lon'] = max_lon
    OBJ['max_lat'] = max_lat
    OBJ['x'] = centroid_x
    OBJ['y'] = centroid_y
    OBJ['n_points'] = sizes
    OBJ['area'] = area

    OBJ['amean_inst_field'] = amean_instantaneous_field
    OBJ['amean_running_field'] = amean_running_field
    OBJ['amean_filtered_running_field'] = amean_filtered_running_field
    OBJ['min_inst_field'] = min_instantaneous_field
    OBJ['min_running_field'] = min_running_field
    OBJ['min_filtered_running_field'] = min_filtered_running_field
    OBJ['max_inst_field'] = max_instantaneous_field
    OBJ['max_running_field'] = max_running_field
    OBJ['max_filtered_running_field'] = max_filtered_running_field

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


def read_lp_object_properties(objid, objdir, property_list, verbose=False, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):

    dt1 = get_objid_datetime(objid)
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


def get_latest_lp_object_time(objdir, level=3):
    #obj_file_list = sorted(glob.glob((objdir + "/????/??/????????/*.nc")))
    obj_file_list = sorted(glob.glob((objdir + "/*"*level + "/*.nc")))

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


def calc_overlapping_points(objid1, objid2, objdir, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):

    dt1 = get_objid_datetime(objid1)
    dt2 = get_objid_datetime(objid2)

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


def init_lpt_group_array(dt_list, objdir, fmt = "/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):
    """
    "LPT" is a 2-D "group" array (np.int64) with columns: [timestamp, objid, lpt_group_id, begin_point, end_point, split_point]
    -- timestamp = Linux time stamp (e.g., seconds since 0000 UTC 1970-1-1)
    -- objid = LP object id (YYYYMMDDHHnnnn)
    -- lpt_group_id = LPT group id, connected LP objects have a common LPT group id.
    -- begin point = 1 if it is the beginning of a track. 0 otherwise.
    -- end point = 1 if no tracks were connected to it, 0 otherwise.
    -- split point = 1 if split detected, 0 otherwise.

    BRANCHES is a 1-D native Python list with native Python int values.
    This is needed because BRANCHES is bitwise, and there can be more than 64 branches in a group.
    -- branches = bitwise binary starts from 1 at each branch. Mergers will have separate branch numbers.
                   overlapping portions will have multiple branch numbers associated with them.
    """

    LPT = []  # Create this as a list, then conert to np.int64 array.
    BRANCHES = []  # This will stay as a native Python list.

    for this_dt in dt_list:

        print(this_dt)
        fn = (objdir + this_dt.strftime(fmt))
        print(fn)
        try:
            DS = Dataset(fn)
            try:
                id_list = DS['objid'][:]
            except IndexError:
                id_list = [] # In case of no LPOs at this time.
            DS.close()

            for ii in range(len(id_list)):
                LPT.append([int((this_dt - dt.datetime(1970,1,1,0,0,0)).total_seconds()), int(id_list[ii]), -1, 0, 1, 0])
                BRANCHES.append(int(0))

        except FileNotFoundError:
            print('WARNING: Missing this file!')

    return (np.int_(LPT), BRANCHES)


def lpt_group_array_remove_small_objects(LPT, BRANCHES, options, verbose=False, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):
    """
    LPT comes from the init_lpt_group_array function
    options needs:
    options['min_lp_objects_points']
    """

    # Make copies to avoid immutability weirdness.
    LPT2 = LPT.copy()
    BRANCHES2 = BRANCHES.copy()

    objdir = options['objdir']
    keep_list = np.full(len(LPT2), True, dtype=bool)

    for ii in range(len(LPT2)):
        this_objid = LPT[ii,1]

        dt1 = get_objid_datetime(this_objid)
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

    return (LPT2[keep_list,:], [BRANCHES2[ii]for ii in range(len(BRANCHES2)) if keep_list[ii]])

def calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=False, reversed=False, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):

    """
    usage: LPT = calc_lpt_group_array(LPT, objdir, options)
    Calculate the simple LPT groups.

    options dictionary entries needed:
    options['objdir']
    options['min_overlap_points']
    options['min_overlap_frac']

    "LPT" is a 2-D "group" array (np.int64) with columns: [timestamp, objid, lpt_group_id, begin_point, end_point, split_point]
    -- timestamp = Linux time stamp (e.g., seconds since 00 UTC 1970-1-1)
    -- objid = LP object id (YYYYMMDDHHnnnn)
    -- lpt_group_id = LPT group id, connected LP objects have a common LPT group id.
    -- begin point = 1 if it is the beginning of a track. 0 otherwise.
    -- end point = 1 if no tracks were connected to it, 0 otherwise.
    -- split point = 1 if split detected, 0 otherwise.

    BRANCHES is a 1-D native Python list with native Python int values.
    This is needed because BRANCHES is bitwise, and there can be more than 64 branches in a group.
    -- branches = bitwise binary starts from 1 at each branch. Mergers will have separate branch numbers.
                   overlapping portions will have multiple branch numbers associated with them.
    """

    # Make copies to avoid immutability weirdness.
    LPT = LPT0.copy()
    BRANCHES = BRANCHES0.copy()

    objdir = options['objdir']
    next_lpt_group_id = 0

    time_list = np.unique(LPT[:,0])
    if reversed:
        time_list = time_list[::-1]
        LPT = LPT[::-1,:]
        BRANCHES = BRANCHES[::-1]

    first_time = time_list[0]
    first_time_idx, = np.where(np.abs(LPT[:,0] - first_time) < 0.1)

    for ii in range(len(first_time_idx)):
        LPT[first_time_idx[ii]][2] = next_lpt_group_id
        LPT[first_time_idx[ii]][3] = 1
        BRANCHES[first_time_idx[ii]] = int(1) # This is the first branch of the LPT. (e.g., 2**0)
        next_lpt_group_id += 1

    for tt in range(1,len(time_list)):
        this_time = time_list[tt]
        prev_time = time_list[tt-1]

        if verbose:
            print(dt.datetime(1970,1,1,0,0,0) + dt.timedelta(seconds=int(this_time)), flush=True)

        this_time_idx, = np.where(np.abs(LPT[:,0] - this_time) <= 0)
        prev_time_idx, = np.where(np.abs(LPT[:,0] - prev_time) <= 0)

        already_connected_objid_list = [] # Keep track of previously connected objids, for split detection.


        append_branch_list = [] # For splits, I will have to add branches to previous points.
                                # This will be a list of tuples, each like (jj, branch_to_append).
        for ii in this_time_idx:
            this_objid = LPT[ii][1]
            match = -1
            for jj in prev_time_idx:
                prev_objid = LPT[jj][1]

                n_this, n_prev, n_overlap = calc_overlapping_points(this_objid,prev_objid,objdir, fmt=fmt)

                if n_overlap >= options['min_overlap_points']:
                    match=jj
                if 1.0*n_overlap/n_this > options['min_overlap_frac']:
                    match=jj
                if 1.0*n_overlap/n_prev > options['min_overlap_frac']:
                    match=jj

            if match > -1:  # I found a match! Note: for merging cases,
                            # this will use the *last* one to match.

                LPT[ii][2] = LPT[match][2] # assign the same group.
                LPT[match][4] = 0 # this one got matched, so not a track end point.

                if LPT[match][1] in already_connected_objid_list:

                    ## OK, This is a split situation!

                    LPT[match][5] = 1 # Identify the split point
                    new_branch = get_group_max_branch(LPT, BRANCHES, LPT[match][2]) + 1
                    #LPT[ii][6] = append_branch(LPT[ii][6], new_branch)  # Start with a fresh, new branch.
                    BRANCHES[ii] = append_branch(BRANCHES[ii], new_branch)  # Start with a fresh, new branch.
                                                                      # Note: I start with branches = 0 here.
                                                                      # LPT[ii,6] should be 0 before appending.

                    append_branch_list.append((match,new_branch)) # Append new branch to splitting point.
                    ##############################################################
                    ## Now I need to append the branch to the previous section of track.
                    ## Basically, everything that is:
                    ##  1) In the same group.
                    ##  2) Has a time stamp prior to prev_time (after for reversed tracking case)
                    ##  3) Has overlapping branches with the splitting point.
                    ##############################################################
                    for dddd in np.where(LPT[:,2] == LPT[ii][2])[0]:
                        if reversed:
                            if LPT[dddd][0] > LPT[match][0]:
                                if BRANCHES[dddd] & BRANCHES[match]: # & with int is *bitwise and*.
                                    append_branch_list.append((dddd,new_branch)) # Append new branch.

                        else:
                            if LPT[dddd][0] < LPT[match][0]:
                                if BRANCHES[dddd] & BRANCHES[match]: # & with int is *bitwise and*.
                                    append_branch_list.append((dddd,new_branch)) # Append new branch.

                else:

                    ## Split situation not detected. Continue the previous track.
                    ## (Subsequent passes through the loop may reveal splits)

                    already_connected_objid_list.append(LPT[match][1])
                    BRANCHES[ii] = BRANCHES[match] # Inherit all branches from the previous one.

            else:   # No overlaps with prior tracks. Therefore, this begins a new group.
                LPT[ii][2] = next_lpt_group_id
                LPT[ii][3] = 1 # This is the beginning of a track.
                BRANCHES[ii] = 1 # This is the first branch of the LPT. (e.g., 2**0)
                next_lpt_group_id += 1

        for this_append_tuple in append_branch_list:
            jj = this_append_tuple[0]
            this_branch_to_append = this_append_tuple[1]
            BRANCHES[jj] = append_branch(BRANCHES[jj], this_branch_to_append) #split point inherits matches.

    if reversed:
        LPT = LPT[::-1,:]
        BRANCHES = BRANCHES[::-1]

    return (LPT, BRANCHES)

def append_branch(branches_binary_int, new_branch_to_append):
    """
    Append branch number IF is it not already there.
    Note: Python 3 or operator in integers is a binary or.
    """
    return int(branches_binary_int) | int(2**(new_branch_to_append-1))


def append_branches_binary(branches_binary_int, new_branches_to_append):
    """
    Append branch number IF is it not already there.
    Note: Python 3 or operator in integers is a binary or.
    """
    return int(branches_binary_int) | int(new_branches_to_append)



def max_branch(branches_binary_int):
    return len(bin(int(branches_binary_int))) - 2 # Subtract out the initial "0b" header.


def has_branch(branches_binary_int, branch_num):
    """
    Return true if the branch_num is a part of the branches_binary_int
    e.g., has_branch(3,2) is True. has_branch(3,3) is False (2**2 = 4 is not included)
    NOTE: branch numbers are ONE based.
    """
    branches_binary_str = str(bin(int(branches_binary_int)))
    if (branches_binary_str[-1*int(branch_num)] == '1'):
        return True
    else:
        return False


def remove_branch_from_group(LPT0, BRANCHES0, lpt_group_id, branch_id_int):

    LPT = LPT0.copy()
    BRANCHES = BRANCHES0.copy()

    indices = sorted([x for x in range(len(LPT[:,0])) if LPT[x,2] == lpt_group_id
        and (BRANCHES[x] & 2**int(branch_id_int-1))])

    for ii in indices:
        BRANCHES[ii] -= 2**int(branch_id_int-1)

    return (LPT, BRANCHES)


def get_group_max_branch(LPT, BRANCHES, lpt_group_id):
    """
    branches is of form "0000 ... 001111"
    What I want is the position of the highest "1"
    in the group. e.g., 4 in this case.
    """
    current_max_branch = 0;

    idx_this_lpt_group = np.where(LPT[:,2] == lpt_group_id)[0]
    for this_idx in idx_this_lpt_group:
        current_max_branch = max(current_max_branch, max_branch(BRANCHES[this_idx]))

    return current_max_branch


def get_group_branches_as_list(branches_binary_int):

    branches_binary_str = str(bin(int(branches_binary_int)))[2:]
    branches_binary_str_r = branches_binary_str[::-1]
    return [x + 1 for x in range(len(branches_binary_str_r)) if branches_binary_str_r[x] == '1']


def get_group_branches_as_int_list(branches_binary_int):

    branches_binary_str = str(bin(int(branches_binary_int)))[2:]
    branches_binary_str_r = branches_binary_str[::-1]
    return [2**x for x in range(len(branches_binary_str_r)) if branches_binary_str_r[x] == '1']


def branches_binary_str4(branches_binary_int):
    """
    String of binary 0 and 1s with break points every 4 digits.
    """
    str0 = str(bin(int(branches_binary_int)))[2:]
    str_length0 = len(str0)
    str_length1 = int(4 * np.ceil(str_length0/4)) # Get new length -- a multiple of 4.
    str1 = str0.zfill(str_length1) # Append zeros to the left if needed.
    str_pieces = [str1[x:x+4] for x in range(0,len(str1),4)]
    return " ".join(str_pieces)


def overlap_forward_backward(LPT1, LPT2, BRANCHES1, BRANCHES2, options, verbose=True):
    """
    LPT1 is meant to be forward.
    LPT2 is meant to be backwards.
    Begin and end points are from LPT1 (forward).
    """

    LPT3 = LPT1.copy()
    BRANCHES3 = BRANCHES1.copy()
    unique_lpt_groups2 = np.unique(LPT2[:,2]) # Does not change.

    more_to_do = True
    while more_to_do:
        more_to_do = False

        unique_lpt_groups3 = np.unique(LPT3[:,2]) # Changes each iteration.

        for this_lpt_group3 in unique_lpt_groups3:
            print((str(this_lpt_group3) + ' of ' + str(np.max(unique_lpt_groups3))), flush=True)
            #idx_this_lpt_group3 = np.where(get_col(LPT3,2) == this_lpt_group3)[0]
            idx_this_lpt_group3 = [x for x in range(len(LPT3)) if LPT3[x][2] == this_lpt_group3]

            #objid_in_this_lpt_group3 = LPT3[idx_this_lpt_group3 , 1]
            objid_in_this_lpt_group3 = LPT3[idx_this_lpt_group3, 1]
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

    LPT3[:,3] = LPT1[:,3] * LPT2[:,4]  ## Begin - begin (end) points in forward (backwards)
    LPT3[:,4] = LPT1[:,4] * LPT2[:,3]  ## End - end (begin) poings in forward (backwards)
    LPT3[:,5] = LPT1[:,5] + LPT2[:,5]  ## Split -- pick up splits identified by forward AND backwards.
    BRANCHES3 = [0 for ii in range(len(BRANCHES3))]

    print('Re-order and return.')
    return (reorder_LPT_group_id(LPT3), BRANCHES3)


def lpt_group_id_separate_branches(LPT, BRANCHES, options, verbose=True, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):

    ##########################################################
    ## Overlapping branches stuff.
    ## Using the "brute force" method of going through time stamp by time stamp.
    ##########################################################

    LPT3 = LPT.copy()
    BRANCHES3 = BRANCHES.copy()

    print('Managing overlapping branches.')

    for this_lpt_group in np.unique(LPT3[:,2]):

        print('Group: ', str(this_lpt_group))
        this_lpt_group_idx_all = np.where(LPT3[:,2] == this_lpt_group)[0]
        this_lpt_time_stamps_all = np.sort(np.unique(LPT3[this_lpt_group_idx_all,0]))

        next_new_branch = 0

        this_lpt_group_idx_begin_points = np.where(np.logical_and(LPT3[:,2] == this_lpt_group, LPT3[:,3] == 1))[0]

        nn = 0
        for begin_idx in this_lpt_group_idx_begin_points:
            nn += 1
            print(("- Begin node #" + str(nn) + " of "
                    + str(len(this_lpt_group_idx_begin_points))))

            ## Iterate until I have connected this initial node point
            ##   with each of the ending node points.
            more_to_do = True # If I find a split, then there is more to do for this starting node.
                              # (mergers are handled as separate starting nodes.)
            already_matched_indices = [] # Keep track of which indices have already been assigned a branch.
            niter = 0
            while more_to_do:
                more_to_do = False
                niter += 1

                ## Seed a new branch.
                BRANCHES3[begin_idx] = BRANCHES3[begin_idx] | 2**next_new_branch # Bitwise or

                current_branch = 1 * next_new_branch # Make sure assignment is by value.
                next_new_branch += 1
                print("-- Branch #" + str(next_new_branch))

                prev_objid = LPT3[begin_idx,1]

                for this_time_stamp in this_lpt_time_stamps_all:
                    if this_time_stamp > LPT3[begin_idx,0] + 0.01:

                        ## Looping forward in time. search for matches with previous branches.
                        lpt_indices_at_this_time = np.where(np.logical_and(LPT3[:,2] == this_lpt_group, LPT3[:,0] == this_time_stamp))[0]
                        matches_this_branch = 0 * lpt_indices_at_this_time


                        lpt_indices_at_this_time_matches = []


                        for ii in range(len(lpt_indices_at_this_time)):
                            this_objid = LPT3[lpt_indices_at_this_time[ii],1]
                            n_this, n_prev, n_overlap = calc_overlapping_points(this_objid,prev_objid,options['objdir'], fmt=fmt)
                            match = False
                            if n_overlap >= options['min_overlap_points']:
                                match = True
                            if 1.0*n_overlap/n_this > options['min_overlap_frac']:
                                match = True
                            if 1.0*n_overlap/n_prev > options['min_overlap_frac']:
                                match = True

                            if match:
                                matches_this_branch[ii] = 1 # I should always find at least 1 match here.
                                if not lpt_indices_at_this_time[ii] in already_matched_indices:
                                    lpt_indices_at_this_time_matches.append(lpt_indices_at_this_time[ii])

                        if len(lpt_indices_at_this_time_matches) < 1:
                            continue # Should only happen if it is a center jump situation.

                        move_forward_idx = lpt_indices_at_this_time_matches[0] # Always take the first (e.g., leftmost) available direction.

                        BRANCHES3[move_forward_idx] = BRANCHES3[move_forward_idx] | int(2**current_branch) # Bitwise or

                        ## If I hit splitting node, take this branch out of "already matched indices".
                        if len(lpt_indices_at_this_time_matches) > 1 and LPT3[prev_idx,5] == 1:
                            idx_with_this_branch = [ x for x in this_lpt_group_idx_all if BRANCHES3[x] & int(2**current_branch) ]
                            for kk in idx_with_this_branch:
                                if kk in already_matched_indices:
                                    already_matched_indices.remove(kk)
                            more_to_do = True

                        if len(lpt_indices_at_this_time_matches) > 1:
                            already_matched_indices.append(move_forward_idx)

                        prev_objid = 1 * LPT3[move_forward_idx,1]
                        prev_idx = 1 * move_forward_idx

                        ## Break out of TIME loop if I have found an end node point.
                        if LPT3[move_forward_idx,4] == 1:
                            break

                #if niter > 0:
                #    more_to_do = False
                #print(more_to_do)

    return (LPT3, BRANCHES3)


def lpt_group_array_allow_center_jumps(LPT, options, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):
    """
    Check duration of "end" (e.g., "this") to "start" points (e.g., "other"), and connect if less than
    center_jump_max_hours.

    NOTE: This does not deal with branches.
    """
    LPT2 = LPT.copy() # Make a copy so I don't inadvertantly over-write the input LPT!

    more_to_do = True
    while more_to_do:
        more_to_do = False

        unique_lpt_groups = np.unique(LPT2[:,2])
        lpt_indices_to_keep = np.array([])

        for this_lpt_group in range(len(unique_lpt_groups)):                # 1
            this_group_all_idx = np.where(LPT2[:,2] == this_lpt_group)[0]
            this_group_end_idx = np.where(np.logical_and(LPT2[:,2] == this_lpt_group, LPT2[:,4] > 0.5))[0]

            other_lpt_group_list = unique_lpt_groups[np.logical_not(unique_lpt_groups == this_lpt_group)]

            for other_lpt_group in other_lpt_group_list:                    # 2
                other_group_all_idx = np.where(LPT2[:,2] == other_lpt_group)[0]
                other_group_begin_idx = np.where(np.logical_and(LPT2[:,2] == other_lpt_group, LPT2[:,3] > 0.5))[0]

                for other_idx in other_group_begin_idx:                     # 3
                    other_objid = LPT2[other_idx,1]
                    other_begin_timestamp = LPT2[other_idx, 0]

                    for this_idx in this_group_end_idx:                     # 4
                        this_objid = LPT2[this_idx,1]
                        this_end_timestamp = LPT2[this_idx, 0]

                        tsdiff = (other_begin_timestamp - this_end_timestamp)

                        if (tsdiff > 1.0 and tsdiff <  options['center_jump_max_hours'] * 3600.0 + 1.0):

                            ## If I got here, the timing is OK for a center jump.
                            ## Now, check the overlapping criteria.
                            n_this, n_prev, n_overlap = calc_overlapping_points(this_objid,other_objid,options['objdir'], fmt=fmt)
                            match = False
                            if n_overlap >= options['min_overlap_points']:
                                match=True
                            if 1.0*n_overlap/n_this > options['min_overlap_frac']:
                                match=True
                            if 1.0*n_overlap/n_prev > options['min_overlap_frac']:
                                match=True

                            if match:
                                ## I found a center jump! Add it to "this group"
                                LPT2[other_group_all_idx, 2] = this_lpt_group
                                LPT2[other_idx, 3] = 0 # No longer the beginning of the track.
                                LPT2[this_idx, 4] = 0  # No longer the end of the track.

                                more_to_do = True
                                break                                       # 4

                    if more_to_do:
                        break                                               # 3

                if more_to_do:
                    break                                                   # 2

            if more_to_do:
                break                                                       # 1

    return reorder_LPT_group_id(LPT2)



def lpt_branches_indices_list(LPT, BRANCHES, group_id, branch_id_int):
    """
    Return the ROW indices that are in branch_id_int for group group_id.
    branch IDs are INTEGER values, not binary.
    Needs LPT (2-d LPT array) and group_id as input.
    """
    return sorted([x for x in range(len(LPT[:,0])) if LPT[x,2] == group_id
        and (BRANCHES[x] & 2**int(branch_id_int-1))])


def lpt_branches_difference(LPT, BRANCHES, group_id, branch_id_int1, branch_id_int2):
    """
    Return the ROW indices that are in branch_id2 but NOT in branch_id1.
    branch IDs are INTEGER values, not binary.
    Needs LPT (2-d LPT array) and group_id as input.
    """
    idx1 = [x for x in range(len(LPT[:,0])) if LPT[x,2] == group_id and (BRANCHES[x] & 2**int(branch_id_int1-1))]
    idx2 = [x for x in range(len(LPT[:,0])) if LPT[x,2] == group_id and (BRANCHES[x] & 2**int(branch_id_int2-1))]
    idx_diff = sorted(list(set(idx2) - set(idx1)))
    return idx_diff



def lpt_branches_intersection(LPT, BRANCHES, group_id, branch_id_int1, branch_id_int2):
    """
    Return the ROW indices that are both branch_id_int1 and branch_id_int2.
    branch IDs are INTEGER values, not binary.
    Needs LPT (2-d LPT array) and group_id as input.
    """
    idx1 = [x for x in range(len(LPT[:,0])) if LPT[x,2] == group_id and ((BRANCHES[x] & 2**int(branch_id_int1-1)) and (BRANCHES[x] & 2**int(branch_id_int2-1)))]
    idx_intersect = sorted(list(set(idx1)))
    return idx_intersect




def lpt_split_and_merge(LPT0, BRANCHES0, merge_split_options):

    LPT = LPT0.copy()
    BRANCHES = BRANCHES0.copy()

    for this_group in np.unique(LPT[:,2]):
        print("!!!!!!!!!!!!!!!!!!!!!!!!!  Group #" + str(this_group) + "  !!!!!!!!!!!!!!!!!!!!!!!!!", flush=True)

        more_to_do = True
        niter=0

        ########################################################################
        ## Split and recombine cases. ##########################################
        ########################################################################
        print('------------------------------------------------')
        print('------------------------------------------------')
        print("Split and recombine (Rejoin 1).")
        print('------------------------------------------------')
        print('------------------------------------------------', flush=True)

        while more_to_do:
            more_to_do = False
            niter+=1

            #if niter > 5:
            #    break

            print('------------------------')
            print('Iteration #' + str(niter))
            print('------------------------', flush=True)

            branch_list = get_branches_in_lpt_group(LPT, BRANCHES, this_group)
            print("Unsorted branch list: " + str(branch_list))
            if len(branch_list) > 1:

                ## Put in order by duration.
                branch_durations = []

                for this_branch in branch_list:

                    lpt_all1 = lpt_branches_indices_list(LPT, BRANCHES, this_group, this_branch)
                    if len(lpt_all1) < 2:
                        continue
                    dt1all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all1]))
                    dt1all_begin = dt1all_list[0]
                    dt1all_end = dt1all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])

                    branch_durations += [(dt1all_end - dt1all_begin).total_seconds()/3600.0]

                branch_list_sorted = [branch_list[x] for x in np.argsort(branch_durations)]
                print("Sorted branch list: " + str(branch_list_sorted))

                for this_branch in branch_list_sorted:

                    #Which branch do you have the most intersection with?
                    other_branch = -1
                    max_intersect_duration = -1
                    for other_branch_test in branch_list_sorted:

                        if this_branch == other_branch_test:
                            continue

                        intersection_with_other_branch_test = lpt_branches_intersection(LPT,  BRANCHES, this_group, this_branch, other_branch_test)
                        if len(intersection_with_other_branch_test) > max_intersect_duration:
                            max_intersect_duration = len(intersection_with_other_branch_test)
                            other_branch = other_branch_test

                    if max_intersect_duration > 0:

                        print(str(this_branch) + ' with ' + str(other_branch) + '.')

                        lpt_all1 = lpt_branches_indices_list(LPT,  BRANCHES, this_group, this_branch)
                        if len(lpt_all1) < 2:
                            continue
                        dt1all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all1]))
                        dt1all_begin = dt1all_list[0]
                        dt1all_end = dt1all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])

                        lpt_all2 = lpt_branches_indices_list(LPT,  BRANCHES, this_group, other_branch)
                        if len(lpt_all2) < 2:
                            continue
                        dt2all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all2]))
                        dt2all_begin = dt2all_list[0]
                        dt2all_end = dt2all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])


                        lpt_diff1 = lpt_branches_difference(LPT, BRANCHES, this_group, this_branch, other_branch)
                        if len(lpt_diff1) > 0:

                            dt1_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_diff1]))
                            dt1_begin = dt1_list[0]
                            dt1_end = dt1_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])
                            dur1 = (dt1_end - dt1_begin).total_seconds()/3600.0

                        else:

                            dt1_begin = None
                            dt1_end = None
                            dur1 = 0.0

                        lpt_diff2 = lpt_branches_difference(LPT, BRANCHES, this_group, other_branch, this_branch)
                        if len(lpt_diff2) > 0:

                            dt2_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_diff2]))
                            dt2_begin = dt2_list[0]
                            dt2_end = dt2_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])
                            dur2 = (dt2_end - dt2_begin).total_seconds()/3600.0

                        else:

                            dt2_begin = None
                            dt2_end = None
                            dur2 = 0.0

                        ## Check to see if I have any difference in the branches.
                        ## NOTE: If branches have already been merged, this should not get triggered again.
                        if (len(lpt_diff1) > 0 or len(lpt_diff2) > 0):
                            ## If the difference is embeded in the intersection, it is a split-them-recombine case.
                            ## The two LPT branches are to be merged in to one.

                            if dt1_begin is None:
                                dt1_begin = dt2_begin #Note: datetimes are immutable.
                                dt1_end = dt2_end #Note: datetimes are immutable.
                            if dt2_begin is None:
                                dt2_begin = dt1_begin #Note: datetimes are immutable.
                                dt2_end = dt1_end #Note: datetimes are immutable.


                            ## Make sure the separated portion is *not* at the beginning or end of either track.
                            if (dt1_begin > dt1all_begin and dt1_end < dt1all_end) and (dt2_begin > dt2all_begin and dt2_end < dt2all_end):
                                print("Split and Re-combine.")
                                print("--> Combine these LPT branches.")


                                # Remove the smaller branch.
                                branches_to_remove = get_group_branches_as_list(BRANCHES[lpt_diff2[0]])
                                for branch_to_remove in branches_to_remove:
                                    LPT, BRANCHES = remove_branch_from_group(LPT, BRANCHES, this_group, branch_to_remove)

                                # Assign those LPT group array indices to the larger branch.
                                for jj in lpt_diff2:
                                    for kk in lpt_diff1:
                                        ## Only "inherit" larger branches for the relevant time span.
                                        kkdt = get_objid_datetime(LPT[kk,1])
                                        if kkdt >= dt2_begin and kkdt <= dt2_end:
                                            #LPT[jj,6] = int(LPT[jj,6]) | int(LPT[kk,6])
                                            BRANCHES[jj] = BRANCHES[jj] | BRANCHES[kk]

                                more_to_do = True


                    if more_to_do:
                        #more_to_do = False
                        break


        ########################################################################
        ## Splits and mergers. #################################################
        ########################################################################
        print('------------------------------------------------')
        print('------------------------------------------------')
        print('Splits and mergers (Rejoin 2).')
        print('------------------------------------------------')
        print('------------------------------------------------', flush=True)

        more_to_do = True

        while more_to_do:
            more_to_do = False
            niter+=1

            #if niter > 5:
            #    break

            print('------------------------')
            print('Iteration #' + str(niter))
            print('------------------------', flush=True)

            branch_list = get_branches_in_lpt_group(LPT, BRANCHES, this_group)
            print("Unsorted branch list: " + str(branch_list))
            if len(branch_list) > 1:

                ## Put in order by duration.
                branch_durations = []

                for this_branch in branch_list:

                    lpt_all1 = lpt_branches_indices_list(LPT,  BRANCHES, this_group, this_branch)
                    if len(lpt_all1) < 2:
                        continue
                    dt1all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all1]))
                    dt1all_begin = dt1all_list[0]
                    dt1all_end = dt1all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])

                    branch_durations += [(dt1all_end - dt1all_begin).total_seconds()/3600.0]

                branch_list_sorted = [branch_list[x] for x in np.argsort(branch_durations)]
                print("Sorted branch list: " + str(branch_list_sorted))

                for this_branch in branch_list_sorted:

                    #Which branch do you have the most intersection with?
                    other_branch = -1
                    max_intersect_duration = -1
                    for other_branch_test in branch_list_sorted:

                        if this_branch == other_branch_test:
                            continue

                        intersection_with_other_branch_test = lpt_branches_intersection(LPT,  BRANCHES, this_group, this_branch, other_branch_test)
                        if len(intersection_with_other_branch_test) > max_intersect_duration:
                            max_intersect_duration = len(intersection_with_other_branch_test)
                            other_branch = other_branch_test

                    if max_intersect_duration > 0:

                        print(str(this_branch) + ' with ' + str(other_branch) + '.')

                        lpt_all1 = lpt_branches_indices_list(LPT,  BRANCHES, this_group, this_branch)
                        if len(lpt_all1) < 2:
                            continue
                        dt1all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all1]))
                        dt1all_begin = dt1all_list[0]
                        dt1all_end = dt1all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])

                        lpt_all2 = lpt_branches_indices_list(LPT,  BRANCHES, this_group, other_branch)
                        if len(lpt_all2) < 2:
                            continue
                        dt2all_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_all2]))
                        dt2all_begin = dt2all_list[0]
                        dt2all_end = dt2all_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])


                        lpt_diff1 = lpt_branches_difference(LPT, BRANCHES, this_group, this_branch, other_branch)
                        if len(lpt_diff1) > 0:

                            #print('Diff1: ' + str(lpt_diff1))
                            dt1_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_diff1]))
                            dt1_begin = dt1_list[0]
                            dt1_end = dt1_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])
                            dur1 = (dt1_end - dt1_begin).total_seconds()/3600.0

                        else:

                            dt1_begin = None
                            dt1_end = None
                            dur1 = 0.0

                        lpt_diff2 = lpt_branches_difference(LPT, BRANCHES, this_group, other_branch, this_branch)
                        if len(lpt_diff2) > 0:
                            #print('Diff2: ' + str(lpt_diff2))
                            dt2_list = sorted(np.unique([get_objid_datetime(LPT[xx,1]) for xx in lpt_diff2]))
                            dt2_begin = dt2_list[0]
                            dt2_end = dt2_list[-1] #get_objid_datetime(LPT[np.max(lpt_diff1),1])
                            dur2 = (dt2_end - dt2_begin).total_seconds()/3600.0

                        else:

                            dt2_begin = None
                            dt2_end = None
                            dur2 = 0.0

                        ## Check to see if I have any difference in the branches.
                        ## NOTE: If branches have already been merged, this should not get triggered again.
                        if (len(lpt_diff1) > 0 or len(lpt_diff2) > 0):
                        #if (len(lpt_diff1) > 0 and len(lpt_diff2) > 0):

                            ## If the difference is embeded in the intersection, it is a split-them-recombine case.
                            ## The two LPT branches are to be merged in to one.

                            if dt1_begin is None:
                                dt1_begin = dt2_begin #Note: datetimes are immutable.
                                dt1_end = dt2_end #Note: datetimes are immutable.
                            if dt2_begin is None:
                                dt2_begin = dt1_begin #Note: datetimes are immutable.
                                dt2_end = dt1_end #Note: datetimes are immutable.


                            print("Merger or Split.")
                            if min(dur1,dur2) > merge_split_options['split_merger_min_hours'] + 0.1:
                                print("--> Retain both branches.")
                            else:
                                print("--> Combine these LPT branches.")
                                if dur1 > dur2:

                                    print('1 > 2')

                                    # Remove the smaller branch.
                                    branches_to_remove = get_group_branches_as_list(BRANCHES[lpt_diff2[0]])
                                    for branch_to_remove in branches_to_remove:
                                        LPT, BRANCHES = remove_branch_from_group(LPT, BRANCHES, this_group, branch_to_remove)

                                    # Assign those LPT group array indices to the larger branch.
                                    for jj in lpt_diff2:
                                        for kk in lpt_diff1:
                                            ## Only "inherit" larger branches for the relevant time span.
                                            kkdt = get_objid_datetime(LPT[kk,1])
                                            if kkdt >= dt2_begin and kkdt <= dt2_end:
                                                #LPT[jj,6] = int(LPT[jj,6]) | int(LPT[kk,6])
                                                BRANCHES[jj] = BRANCHES[jj] | BRANCHES[kk]

                                else:

                                    print('2 > 1')

                                    # Remove the smaller branch.
                                    branches_to_remove = get_group_branches_as_list(BRANCHES[lpt_diff1[0]])
                                    for branch_to_remove in branches_to_remove:
                                        LPT, BRANCHES = remove_branch_from_group(LPT, BRANCHES, this_group, branch_to_remove)

                                    # Assign those LPT group array indices to the larger branch.
                                    for jj in lpt_diff1:
                                        for kk in lpt_diff2:
                                            ## Only "inherit" larger branches for the relevant time span.
                                            kkdt = get_objid_datetime(LPT[kk,1])
                                            if kkdt >= dt1_begin and kkdt <= dt1_end:
                                                #LPT[jj,6] = int(LPT[jj,6]) | int(LPT[kk,6])
                                                BRANCHES[jj] = BRANCHES[jj] | BRANCHES[kk]

                                more_to_do = True
                                break


                    if more_to_do:
                        #more_to_do = False
                        break


    return (LPT, BRANCHES)


def remove_short_lived_systems(LPT, BRANCHES, minimum_duration_hours, latest_datetime = dt.datetime(2063,4,5,0,0,0)):

    """
    Remove short duration LPT groups.
    But only up to the latest_datetime.
    """

    latest_timestamp = latest_datetime.timestamp()
    unique_lpt_groups = np.unique(LPT[:,2])
    lpt_indices_to_keep = np.array([])

    for this_lpt_group in range(len(unique_lpt_groups)):
        this_lpt_group_idx = np.where(LPT[:,2] == this_lpt_group)[0]

        timestamp_this_group = LPT[this_lpt_group_idx,0]
        min_timestamp_this_group = np.min(timestamp_this_group)
        max_timestamp_this_group = np.max(timestamp_this_group)

        if max_timestamp_this_group > latest_timestamp:
            lpt_indices_to_keep = np.append(lpt_indices_to_keep, this_lpt_group_idx)
        else:
            duration_this_group_hours = (max_timestamp_this_group - min_timestamp_this_group) / 3600.0
            print(duration_this_group_hours)
            if duration_this_group_hours > (minimum_duration_hours - 0.01):
                lpt_indices_to_keep = np.append(lpt_indices_to_keep, this_lpt_group_idx)

    LPT2 = LPT[lpt_indices_to_keep.astype('int').tolist(),:]
    BRANCHES2 = [BRANCHES[x] for x in lpt_indices_to_keep.astype('int').tolist()]

    return (reorder_LPT_group_id(LPT2), BRANCHES2)



def reorder_LPT_group_id(LPT0):

    """
    re-order and relabel the LPT group IDs to 0 to N
    where N is the number of unique LPT group IDs.
    """

    LPT = LPT0.copy()
    ## Re-order LPT system groups
    unique_lpt_groups = np.unique(LPT[:,2])

    for jjj in range(len(unique_lpt_groups)):
        LPT[(LPT[:,2] == unique_lpt_groups[jjj]), 2] = jjj

    return LPT


def reorder_LPT_branches(LPT0, BRANCHES0):

    """
    re-order and relabel the LPT branches 1 to N
    Useful when mergers/splits eliminated some of the branches.
    NOTE: This returns a new BRANCHES list. It only uses LPT to
    get the group information, and doesn't touch or return LPT.
    """
    LPT = LPT0.copy()
    BRANCHES = BRANCHES0.copy()
    ## Re-order LPT system groups
    unique_lpt_groups = np.unique(LPT[:,2])

    for this_group in unique_lpt_groups:
        idx_this_lpt_group = [x for x in len(BRANCHES) if LPT[x,2] == this_group]

    return LPT



def calc_lpt_system_group_properties(LPT, options, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):

    unique_lpt_groups = np.unique(LPT[:,2])

    TC_all = []

    for this_group in unique_lpt_groups:
        TC_this = {}
        TC_this['lpt_group_id'] = this_group
        TC_this['lpt_id'] = this_group

        this_lpt_group_idx = np.where(LPT[:,2] == this_group)[0]
        TC_this['objid'] = LPT[this_lpt_group_idx,1]
        TC_this['timestamp'] = np.unique(LPT[this_lpt_group_idx,0])
        TC_this['datetime'] = [dt.datetime(1970,1,1,0,0,0) + dt.timedelta(seconds=int(x)) for x in TC_this['timestamp']]

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

        TC_this['amean_inst_field'] = np.zeros(len(TC_this['timestamp']))
        TC_this['amean_running_field'] = np.zeros(len(TC_this['timestamp']))
        TC_this['amean_filtered_running_field'] = np.zeros(len(TC_this['timestamp']))
        TC_this['min_inst_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['min_running_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['min_filtered_running_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['max_inst_field'] = -999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['max_running_field'] = -999.0 * np.ones(len(TC_this['timestamp']))
        TC_this['max_filtered_running_field'] = -999.0 * np.ones(len(TC_this['timestamp']))


        ## Loop over unique time stamps.
        for tt in range(len(TC_this['timestamp'])):
            idx_for_this_time = np.where(np.logical_and(
                LPT[:,0] == TC_this['timestamp'][tt],
                LPT[:,2] == this_group))[0]

            for this_objid in LPT[idx_for_this_time,1]:

                OBJ = read_lp_object_properties(this_objid, options['objdir']
                        , ['centroid_lon','centroid_lat','area','pixels_x','pixels_y'
                        ,'min_lon','max_lon','min_lat','max_lat'
                        ,'amean_inst_field','amean_running_field','max_inst_field','max_running_field'
                        ,'min_inst_field','min_running_field','min_filtered_running_field'
                        ,'amean_filtered_running_field','max_filtered_running_field'], fmt=fmt)

                TC_this['nobj'][tt] += 1
                TC_this['area'][tt] += OBJ['area']
                TC_this['centroid_lon'][tt] += OBJ['centroid_lon'] * OBJ['area']
                TC_this['centroid_lat'][tt] += OBJ['centroid_lat'] * OBJ['area']

                TC_this['min_lon'][tt] = min((TC_this['min_lon'][tt], OBJ['min_lon']))
                TC_this['min_lat'][tt] = min((TC_this['min_lat'][tt], OBJ['min_lat']))
                TC_this['max_lon'][tt] = max((TC_this['max_lon'][tt], OBJ['max_lon']))
                TC_this['max_lat'][tt] = max((TC_this['max_lat'][tt], OBJ['max_lat']))

                TC_this['amean_inst_field'][tt] += OBJ['amean_inst_field'] * OBJ['area']
                TC_this['amean_running_field'][tt] += OBJ['amean_running_field'] * OBJ['area']
                TC_this['amean_filtered_running_field'][tt] += OBJ['amean_filtered_running_field'] * OBJ['area']
                TC_this['min_inst_field'][tt] = min((TC_this['min_inst_field'][tt], OBJ['min_inst_field']))
                TC_this['min_running_field'][tt] = min((TC_this['min_running_field'][tt], OBJ['min_running_field']))
                TC_this['min_filtered_running_field'][tt] = min((TC_this['min_filtered_running_field'][tt], OBJ['min_filtered_running_field']))
                TC_this['max_inst_field'][tt] = max((TC_this['max_inst_field'][tt], OBJ['max_inst_field']))
                TC_this['max_running_field'][tt] = max((TC_this['max_running_field'][tt], OBJ['max_running_field']))
                TC_this['max_filtered_running_field'][tt] = max((TC_this['max_filtered_running_field'][tt], OBJ['max_filtered_running_field']))

            TC_this['centroid_lon'][tt] /= TC_this['area'][tt]
            TC_this['centroid_lat'][tt] /= TC_this['area'][tt]

            TC_this['amean_inst_field'][tt] /= TC_this['area'][tt]
            TC_this['amean_running_field'][tt] /= TC_this['area'][tt]
            TC_this['amean_filtered_running_field'][tt] /= TC_this['area'][tt]

        ## Least squares linear fit for propagation speed.
        Pzonal = np.polyfit(TC_this['timestamp'],TC_this['centroid_lon'],1)
        TC_this['zonal_propagation_speed'] = Pzonal[0] * 111000.0  # deg / s --> m / s

        Pmeridional = np.polyfit(TC_this['timestamp'],TC_this['centroid_lat'],1)
        TC_this['meridional_propagation_speed'] = Pmeridional[0] * 111000.0  # deg / s --> m / s

        TC_all.append(TC_this)

    return TC_all


def separate_lpt_system_branches(LPTfb, LPTf, LPTb, options):
    """
    TODO: NOT YET IMPLEMENTED.
    """
    LPT_with_branches = LPTfb.copy() # Start with forward/backwards merged system.


def get_branches_in_lpt_group(LPT, BRANCHES, lpt_group_id):
    branches = 0
    for idx in [x for x in range(len(LPT[:,0])) if LPT[x,2] == lpt_group_id]:
        branches = branches | BRANCHES[idx]
    return  get_group_branches_as_list(branches)


def calc_lpt_system_group_properties_with_branches(LPT, BRANCHES, options, fmt="/%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc"):
    unique_lpt_groups = np.unique(LPT[:,2])

    TC_all = []

    for this_group in unique_lpt_groups:

        this_branch_list = get_branches_in_lpt_group(LPT, BRANCHES, this_group)

        for this_branch in this_branch_list:
            TC_this = {}
            TC_this['lpt_group_id'] = this_group
            TC_this['lpt_id'] = this_group + this_branch / 100.0

            this_lpt_group_idx = [x for x in range(len(LPT[:,0])) if LPT[x,2] == this_group and (BRANCHES[x] & 2**(this_branch-1))]
            TC_this['objid'] = LPT[this_lpt_group_idx,1]
            TC_this['timestamp'] = np.unique(LPT[this_lpt_group_idx,0])
            TC_this['datetime'] = [dt.datetime(1970,1,1,0,0,0) + dt.timedelta(seconds=int(x)) for x in TC_this['timestamp']]

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

            TC_this['amean_inst_field'] = np.zeros(len(TC_this['timestamp']))
            TC_this['amean_running_field'] = np.zeros(len(TC_this['timestamp']))
            TC_this['amean_filtered_running_field'] = np.zeros(len(TC_this['timestamp']))
            TC_this['min_inst_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
            TC_this['min_running_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
            TC_this['min_filtered_running_field'] = 999.0 * np.ones(len(TC_this['timestamp']))
            TC_this['max_inst_field'] = -999.0 * np.ones(len(TC_this['timestamp']))
            TC_this['max_running_field'] = -999.0 * np.ones(len(TC_this['timestamp']))
            TC_this['max_filtered_running_field'] = -999.0 * np.ones(len(TC_this['timestamp']))



            ## Loop over time.
            for tt in range(len(TC_this['timestamp'])):
                #idx_for_this_time = np.where(np.logical_and(
                #    LPT[:,0] == TC_this['timestamp'][tt],
                #    LPT[:,2] == this_group))[0]

                idx_for_this_time = [x for x in range(len(LPT[:,0])) if (LPT[x,0] == TC_this['timestamp'][tt]) and (LPT[x,2] == this_group) and (BRANCHES[x] & 2**(this_branch-1))]

                for this_objid in LPT[idx_for_this_time,1]:

                    OBJ = read_lp_object_properties(this_objid, options['objdir']
                            , ['centroid_lon','centroid_lat','area','pixels_x','pixels_y'
                            ,'min_lon','max_lon','min_lat','max_lat'
                            ,'amean_inst_field','amean_running_field','max_inst_field','max_running_field'
                            ,'min_inst_field','min_running_field','min_filtered_running_field'
                            ,'amean_filtered_running_field','max_filtered_running_field'], fmt=fmt)

                    TC_this['nobj'][tt] += 1
                    TC_this['area'][tt] += OBJ['area']
                    TC_this['centroid_lon'][tt] += OBJ['centroid_lon'] * OBJ['area']
                    TC_this['centroid_lat'][tt] += OBJ['centroid_lat'] * OBJ['area']

                    TC_this['min_lon'][tt] = min((TC_this['min_lon'][tt], OBJ['min_lon']))
                    TC_this['min_lat'][tt] = min((TC_this['min_lat'][tt], OBJ['min_lat']))
                    TC_this['max_lon'][tt] = max((TC_this['max_lon'][tt], OBJ['max_lon']))
                    TC_this['max_lat'][tt] = max((TC_this['max_lat'][tt], OBJ['max_lat']))

                    TC_this['amean_inst_field'][tt] += OBJ['amean_inst_field'] * OBJ['area']
                    TC_this['amean_running_field'][tt] += OBJ['amean_running_field'] * OBJ['area']
                    TC_this['amean_filtered_running_field'][tt] += OBJ['amean_filtered_running_field'] * OBJ['area']
                    TC_this['min_inst_field'][tt] = min((TC_this['min_inst_field'][tt], OBJ['min_inst_field']))
                    TC_this['min_running_field'][tt] = min((TC_this['min_running_field'][tt], OBJ['min_running_field']))
                    TC_this['min_filtered_running_field'][tt] = min((TC_this['min_filtered_running_field'][tt], OBJ['min_filtered_running_field']))
                    TC_this['max_inst_field'][tt] = max((TC_this['max_inst_field'][tt], OBJ['max_inst_field']))
                    TC_this['max_running_field'][tt] = max((TC_this['max_running_field'][tt], OBJ['max_running_field']))
                    TC_this['max_filtered_running_field'][tt] = max((TC_this['max_filtered_running_field'][tt], OBJ['max_filtered_running_field']))

                TC_this['centroid_lon'][tt] /= TC_this['area'][tt]
                TC_this['centroid_lat'][tt] /= TC_this['area'][tt]

                TC_this['amean_inst_field'][tt] /= TC_this['area'][tt]
                TC_this['amean_running_field'][tt] /= TC_this['area'][tt]
                TC_this['amean_filtered_running_field'][tt] /= TC_this['area'][tt]


            ## Least squares linear fit for propagation speed.
            Pzonal = np.polyfit(TC_this['timestamp'],TC_this['centroid_lon'],1)
            TC_this['zonal_propagation_speed'] = Pzonal[0] * 111000.0  # deg / s --> m / s

            Pmeridional = np.polyfit(TC_this['timestamp'],TC_this['centroid_lat'],1)
            TC_this['meridional_propagation_speed'] = Pmeridional[0] * 111000.0  # deg / s --> m / s

            TC_all.append(TC_this)

    return TC_all



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

def plot_lpt_groups_time_lon_text(ax, LPT, BRANCHES, options, text_color='k'):

    objdir = options['objdir']
    dt_min = dt.datetime(1970,1,1,0,0,0) + dt.timedelta(seconds=int(np.min(LPT[:,0])))
    dt_max = dt.datetime(1970,1,1,0,0,0) + dt.timedelta(seconds=int(np.max(LPT[:,0])))

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

        plt.text(lon, dt1, (str(LPT[ii,2]) + ": " + branches_binary_str4(BRANCHES[ii]))
                  , color=this_text_color, zorder=this_zorder, fontsize=6, clip_on=True)

    ax.set_xlim([0.0, 360.0])
    ax.set_ylim([dt_min, dt_max + dt.timedelta(hours=3)])



def float_lpt_id(group, branch):
    """
    Branch is a decimal tacked on to the group ID.
    group 7, branch #1 is 7.01
    group 20, branch 10 is 20.10
    Branch > 100 will give an error message and return np.nan.
    """
    if branch > 99:
        print('ERROR! Branch number > 99.')
        float_id = np.nan
    else:
        float_id = group + branch / 100.0

    return float_id


def plot_timeclusters_time_lon(ax, TIMECLUSTERS, linewidth=2.0):

    for ii in range(len(TIMECLUSTERS)):
        x = TIMECLUSTERS[ii]['centroid_lon']
        y = TIMECLUSTERS[ii]['datetime']
        ax.plot(x, y, 'k', linewidth=linewidth)

        plt.text(x[0], y[0], str(int(ii)), fontweight='bold', color='red', clip_on=True)
        plt.text(x[-1], y[-1], str(int(ii)), fontweight='bold', color='red', clip_on=True)

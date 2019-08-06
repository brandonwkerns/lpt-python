import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import datetime as dt
from context import lpt
import os
import sys

plt.close('all')
verbose=False
#verbose=True

if len(sys.argv) < 2:
    year10 = 2011
else:
    year10 = int(sys.argv[1])

if len(sys.argv) < 3:
    year11 = 2011
else:
    year11 = int(sys.argv[2])

if len(sys.argv) < 4:
    prod = 'tmpa'
else:
    prod = sys.argv[3]


################################################################################

def feature_spread(array_in, npoints):

    ## Use convolution to expand the mask "array_in" a radius of np points.
    ## For this purpose, it takes a 3-D array with the first entry being time.

    array_out = array_in.copy()
    s = array_in.shape

    if type(npoints) is list:
        npx = npoints[0]
        npy = npoints[1]
    else:
        npx = 1*npoints
        npy = 1*npoints

    [circle_array_x, circle_array_y] = np.meshgrid(np.arange(-1*npx,npx+1), np.arange(-1*npy,npy+1))
    circle_array_dist = np.sqrt(np.power(circle_array_x,2) + np.power(circle_array_y * (npx/npy),2))
    circle_array_mask = (circle_array_dist < (npx + 0.1)).astype(np.double)
    circle_array_mask = circle_array_mask / np.sum(circle_array_mask)

    #plt.pcolormesh(circle_array_mask)
    #plt.colorbar()
    #plt.savefig('test.png')


    ## Loop over the times.
    ## For each time, use the convolution to "spread out" the effect of each time's field.
    for tt in range(s[0]):
        array_2d = array_in[tt,:,:]
        array_2d_new = array_2d.copy()
        unique_values = np.unique(array_2d)
        unique_values = unique_values[unique_values > 0]  #take out zero -- it is not a feature.
        for this_value in unique_values:
            starting_mask = (array_2d == this_value).astype(np.double)
            starting_mask_spread = scipy.ndimage.convolve(starting_mask, circle_array_mask, mode='constant')
            array_2d_new[starting_mask_spread > 0.001] = this_value

        array_out[tt,:,:] = array_2d_new

    return array_out

################################################################################


for year1 in range(year10, year11+1):

    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    print(prod + ': ' + str(year1))
    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

    ## For TMPA

    if prod == 'tmpa':
        filter='g20_72h'
        thresh='thresh12'
        interval_hours=3
        accumulation_hours=72
        filter_stdev = 10   # Points in filter width.
        YMDH1_YMDH2='{0:d}060100_{1:d}063021'.format(year1, year1+1)

    ## For CFSR

    elif prod == 'cfsr':
        filter='g10_72h'
        thresh='thresh12'
        interval_hours=6
        accumulation_hours=72
        filter_stdev = 10   # Points in filter width.
        YMDH1_YMDH2='{0:d}060100_{1:d}063018'.format(year1, year1+1)


    ## For ERA5
    elif prod == 'era5':
        filter='g18_72h'
        thresh='thresh12'
        interval_hours=12
        accumulation_hours=72
        filter_stdev = 10   # Points in filter width.
        YMDH1_YMDH2='{0:d}060100_{1:d}063012'.format(year1, year1+1)


    ## For MERRA2
    elif prod == 'merra2':
        filter='g10x8_72h'
        thresh='thresh12'
        interval_hours=6
        accumulation_hours=72
        filter_stdev = [10,8]   # Points in filter width.
        YMDH1_YMDH2='{0:d}060100_{1:d}063018'.format(year1, year1+1)

    else:

        print('Product ' + prod + ' is not valid!')
        break


    lpt_group_file = ('/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/'+prod+'/'+filter+'/'+thresh+'/systems/lpt_systems_tmpa_rt_'+YMDH1_YMDH2+'.group_array.txt')
    lpt_objects_dir = ('/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/'+prod+'/'+filter+'/'+thresh+'/objects')

    dt1 = dt.datetime(year1,6,1,0,0,0)
    dt2 = dt.datetime(year1+1,7,1,0,0,0) - dt.timedelta(hours=interval_hours)
    total_hours = (dt2 - dt1).total_seconds()/3600.0
    mask_times = [dt1 + dt.timedelta(hours=x) for x in np.arange(0,total_hours+1,interval_hours)]

    mask_arrays={} #Start with empty dictionary

    dt_idx = -1
    for this_dt in mask_times:
        dt_idx += 1

        print(this_dt.strftime('%Y%m%d%H'), flush=True)

        fn = (lpt_objects_dir + '/' + this_dt.strftime('%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc'))
        try:
            DS=Dataset(fn)
        except:
            continue


        ## Initialize the mask arrays dictionary if this is the first LP object.
        ## First, I need the grid information. Get this from the first LP object.
        if len(mask_arrays) < 1:
            lon = DS['grid_lon'][:]
            lat = DS['grid_lat'][:]
            mask_arrays['lon'] = DS['grid_lon'][:]
            mask_arrays['lat'] = DS['grid_lat'][:]
            mask_arrays_shape = [len(mask_times), len(lat), len(lon)]
            mask_arrays['mask_at_end_time'] = np.zeros(mask_arrays_shape)
            mask_arrays['mask_with_filter_at_end_time'] = np.zeros(mask_arrays_shape)
            mask_arrays['mask_with_accumulation'] = np.zeros(mask_arrays_shape)
            mask_arrays['mask_with_filter_and_accumulation'] = np.zeros(mask_arrays_shape)

        ##
        ## Get LP Object pixel information.
        ##
        #iii = DS['pixels_x'][nnnn,:].compressed()
        #jjj = DS['pixels_y'][nnnn,:].compressed()
        try:
            iii = DS['pixels_x'][:].compressed()
            jjj = DS['pixels_y'][:].compressed()
            DS.close()

            ##
            ## Fill in the mask information.
            ##

            ## For mask_at_end_time, just use the mask from the objects file.
            mask_arrays['mask_at_end_time'][dt_idx, jjj, iii] = 1

            ## For the mask with accumulation, go backwards and fill in ones.
            n_back = int(accumulation_hours/interval_hours)
            for ttt in range(dt_idx - n_back, dt_idx+1):
                mask_arrays['mask_with_accumulation'][ttt, jjj, iii] = 1

        except:
            pass

    ##
    ## Do filter width spreading.
    ##

    print('Filter width spreading...this may take awhile.', flush=True)
    mask_arrays['mask_with_filter_at_end_time'] = feature_spread(mask_arrays['mask_at_end_time'], filter_stdev)
    mask_arrays['mask_with_filter_and_accumulation'] = feature_spread(mask_arrays['mask_with_accumulation'], filter_stdev)

    ##
    ## Output.
    ##
    fn_out = ('/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/'+prod+'/'+filter+'/'+thresh+'/systems/'+YMDH1_YMDH2+'/lp_objects_mask_'+YMDH1_YMDH2+'.nc')
    os.makedirs('/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/'+prod+'/'+filter+'/'+thresh+'/systems/'+YMDH1_YMDH2, exist_ok=True)

    print('Writing to: ' + fn_out, flush=True)
    DSnew = Dataset(fn_out, 'w')
    DSnew.createDimension('time',len(mask_times))
    DSnew.createDimension('lon',len(lon))
    DSnew.createDimension('lat',len(lat))

    DSnew.createVariable('time','d',('time',))
    DSnew.createVariable('lon','d',('lon',))
    DSnew.createVariable('lat','d',('lat',))

    DSnew['time'][:] = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds()/3600.0 for x in mask_times]
    DSnew['time'].setncattr('units','hours since 1970-1-1 0:0:0')
    DSnew['lon'][:] = lon
    DSnew['lon'].setncattr('units','degrees_east')
    DSnew['lat'][:] = lat
    DSnew['lat'].setncattr('units','degrees_north')

    for mask_var in ['mask_at_end_time','mask_with_filter_at_end_time','mask_with_accumulation','mask_with_filter_and_accumulation']:
        DSnew.createVariable(mask_var,'i',('time','lat','lon'))
        DSnew[mask_var][:] = mask_arrays[mask_var]
        DSnew[mask_var].setncattr('units','1')

    DSnew.close()

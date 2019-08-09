import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import datetime as dt
from context import lpt
import os
import sys
import glob


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

data_dir = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'

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

    ## For WRF

    elif prod == 'wrf':
        filter='g14_72h'
        thresh='thresh14'
        interval_hours=1
        accumulation_hours=72
        filter_stdev = 14   # Points in filter width.
        YMDH1_YMDH2='{0:d}060100_{1:d}063018'.format(year1, year1+1)
        label='ind_20111122_ecmwf_d'


    else:

        print('Product ' + prod + ' is not valid!')
        break




    if prod == 'wrf':
        lpt_systems_file = glob.glob(data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/' + label + '/lpt_systems_'+prod+'_*.nc')[0]
        lpt_group_file = glob.glob(data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/' + label + '/lpt_systems_'+prod+'_*.group_array.txt')[0]
        lpt_objects_dir = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/objects/' + label)
    else:
        lpt_systems_file = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/lpt_systems_'+prod+'_'+YMDH1_YMDH2+'.nc')
        lpt_group_file = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/lpt_systems_'+prod+'_'+YMDH1_YMDH2+'.group_array.txt')
        lpt_objects_dir = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/objects')


    MISSING = -999.0
    FILL_VALUE = MISSING


    ## Read Stitched data.
    DS = Dataset(lpt_systems_file)
    TC={}
    TC['lptid'] = DS['lptid'][:]
    TC['i1'] = DS['lpt_begin_index'][:]
    TC['i2'] = DS['lpt_end_index'][:]
    TC['timestamp_stitched'] = DS['timestamp_stitched'][:]
    TC['datetime'] = [dt.datetime(1970,1,1,0,0,0) + dt.timedelta(hours=int(x)) if x > 100 else None for x in TC['timestamp_stitched']]
    TC['centroid_lon'] = DS['centroid_lon_stitched'][:]
    TC['centroid_lat'] = DS['centroid_lat_stitched'][:]
    TC['area'] = DS['area_stitched'][:]
    TC['max_running_field'] = DS['max_running_field'][:]
    DS.close()


    LPT, BRANCHES = lpt.lptio.read_lpt_systems_group_array(lpt_group_file)

    F = Dataset(lpt_systems_file)
    unique_lpt_ids = np.unique(F['lptid'][:])

    for this_lpt_id in unique_lpt_ids:
        print('Calculating LPT system mask for lptid = ' + str(this_lpt_id) + ' of time period ' + YMDH1_YMDH2 + '.')

        this_group = np.floor(this_lpt_id)
        this_group_lptid_list = sorted([x for x in unique_lpt_ids if np.floor(x) == this_group])

        if np.round( 100.0 * (this_group_lptid_list[0] - this_group)) > 0:
            this_branch = int(2**(np.round( 100.0 * (this_lpt_id - this_group)) - 1))
        else:
            this_branch = int(2**(np.round( 1000.0 * (this_lpt_id - this_group)) - 1))

        print((this_lpt_id, this_branch))
        this_branch_idx = [x for x in range(len(BRANCHES)) if LPT[x,2]==this_group and (BRANCHES[x] & this_branch) > 0] # bitwise and
        #lp_object_id_list = LPT[LPT[:,2]==this_group,1]
        lp_object_id_list = LPT[this_branch_idx,1]


        dt0 = dt.datetime.strptime(str(int(np.min(lp_object_id_list)))[0:10],'%Y%m%d%H') - dt.timedelta(hours=accumulation_hours)
        dt1 = dt.datetime.strptime(str(int(np.max(lp_object_id_list)))[0:10],'%Y%m%d%H')
        duration_hours = int((dt1 - dt0).total_seconds()/3600)
        mask_times = [dt0 + dt.timedelta(hours=x) for x in range(0,duration_hours+1,interval_hours)]
        mask_arrays={} #Start with empty dictionary


        ## Include some basic LPT info for user friendliness.
        lptidx = [ii for ii in range(len(TC['lptid'])) if this_lpt_id == TC['lptid'][ii]][0]
        for var in ['centroid_lon','centroid_lat','area','max_running_field']:
            mask_arrays[var] = MISSING * np.ones(len(mask_times))

        for ttt in range(TC['i1'][lptidx],TC['i2'][lptidx]+1):
            this_time_indx = [ii for ii in range(len(mask_times)) if TC['datetime'][ttt] == mask_times[ii]]
            if len(this_time_indx) > 0:
                for var in ['centroid_lon','centroid_lat','area','max_running_field']:
                    mask_arrays[var][this_time_indx] = TC[var][ttt]


        for lp_object_id in lp_object_id_list:

            nnnn = int(str(int(lp_object_id))[-4:])
            try:
                dt_this = lpt.helpers.get_objid_datetime(lp_object_id) # dt.datetime.strptime(str(int(lp_object_id))[0:10],'%Y%m%d%H')
                dt_idx = [tt for tt in range(len(mask_times)) if dt_this == mask_times[tt]]
            except:
                continue


            if len(dt_idx) < 0:
                print('This time not found in mask time list. Skipping LP object id: ' + str(int(lp_object_id)))
                continue
            elif len(dt_idx) > 1:
                print('Found more than one mask time for this LP object. This should not happen! Skipping it.')
                continue
            else:
                dt_idx = dt_idx[0]

            if prod == 'wrf':
                fn = (lpt_objects_dir + '/' + dt_this.strftime('objects_%Y%m%d%H.nc'))
            else:
                fn = (lpt_objects_dir + '/' + dt_this.strftime('%Y/%m/%Y%m%d/objects_%Y%m%d%H.nc'))
            if verbose:
                print(fn)
            DS=Dataset(fn)

            ## Initialize the mask arrays dictionary if this is the first LP object.
            ## First, I need the grid information. Get this from the first LP object.
            if not 'lon' in mask_arrays:
                lon = DS['grid_lon'][:]
                lat = DS['grid_lat'][:]
                mask_arrays['lon'] = DS['grid_lon'][:]
                mask_arrays['lat'] = DS['grid_lat'][:]
                if prod == 'wrf':
                    ny, nx = mask_arrays['lon'].shape
                    mask_arrays_shape = [len(mask_times), ny, nx]
                else:
                    mask_arrays_shape = [len(mask_times), len(lat), len(lon)]
                mask_arrays['mask_at_end_time'] = np.zeros(mask_arrays_shape)
                mask_arrays['mask_with_filter_at_end_time'] = np.zeros(mask_arrays_shape)
                mask_arrays['mask_with_accumulation'] = np.zeros(mask_arrays_shape)
                mask_arrays['mask_with_filter_and_accumulation'] = np.zeros(mask_arrays_shape)

            ##
            ## Get LP Object pixel information.
            ##
            try:
                iii = DS['pixels_x'][nnnn,:].compressed()
                jjj = DS['pixels_y'][nnnn,:].compressed()

            except:
                DS.close()
                continue

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

        ##
        ## Do filter width spreading.
        ##

        print('Filter width spreading...this may take awhile.', flush=True)
        mask_arrays['mask_with_filter_at_end_time'] = feature_spread(mask_arrays['mask_at_end_time'], filter_stdev)
        mask_arrays['mask_with_filter_and_accumulation'] = feature_spread(mask_arrays['mask_with_accumulation'], filter_stdev)

        ##
        ## Output.
        ##
        if prod == 'wrf':
            fn_out = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/'+label+'/lpt_system_mask_'+label+'.lptid{0:010.4f}.nc'.format(this_lpt_id))
        else:
            fn_out = (data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/'+YMDH1_YMDH2+'/lpt_system_mask_'+YMDH1_YMDH2+'.lptid{0:010.4f}.nc'.format(this_lpt_id))
        os.makedirs(data_dir + '/'+prod+'/'+filter+'/'+thresh+'/systems/'+YMDH1_YMDH2, exist_ok=True)

        print('Writing to: ' + fn_out, flush=True)
        DSnew = Dataset(fn_out, 'w', data_model='NETCDF4', clobber=True)
        DSnew.createDimension('time',len(mask_times))
        if prod == 'wrf':
            ny,nx = lon.shape
            DSnew.createDimension('x',nx)
            DSnew.createDimension('y',ny)
            DSnew.createVariable('lon','f4',('y','x'))
            DSnew.createVariable('lat','f4',('y','x'))
        else:
            DSnew.createDimension('lon',len(lon))
            DSnew.createDimension('lat',len(lat))
            DSnew.createVariable('lon','f4',('lon',))
            DSnew.createVariable('lat','f4',('lat',))


        DSnew.createVariable('time','d',('time',)) # I would like to use u4, but ncview complains about dimension variable being unknown type.

        DSnew.createVariable('centroid_lon','f4',('time',),fill_value=FILL_VALUE)
        DSnew.createVariable('centroid_lat','f4',('time',),fill_value=FILL_VALUE)
        DSnew.createVariable('area','d',('time',),fill_value=FILL_VALUE)
        DSnew.createVariable('max_running_field','f4',('time',),fill_value=FILL_VALUE)

        ts = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds()/3600.0 for x in mask_times]
        print(ts[0:4])
        DSnew['time'][:] = ts
        DSnew['time'].setncattr('units','hours since 1970-1-1 0:0:0')
        DSnew['lon'][:] = lon
        DSnew['lon'].setncattr('units','degrees_east')
        DSnew['lat'][:] = lat
        DSnew['lat'].setncattr('units','degrees_north')


        for mask_var in ['centroid_lon','centroid_lat','area','max_running_field']:
            DSnew[mask_var][:] = mask_arrays[mask_var]

        DSnew['centroid_lon'].setncatts({'units':'degrees_east','long_name':'centroid longitude (0-360)','standard_name':'longitude','note':'Time is end of running mean time.'})
        DSnew['centroid_lat'].setncatts({'units':'degrees_east','long_name':'centroid latitude (-90-00)','standard_name':'latitude','note':'Time is end of running mean time.'})
        DSnew['area'].setncatts({'units':'km2','long_name':'LPT System enclosed area','note':'Time is end of running mean time.'})
        DSnew['max_running_field'].setncatts({'units':'mm day-1','long_name':'LP object running mean rain rate (at end of accum time).','note':'Time is end of running mean time.'})

        for mask_var in ['mask_at_end_time','mask_with_filter_at_end_time','mask_with_accumulation','mask_with_filter_and_accumulation']:
            if prod == 'wrf':
                DSnew.createVariable(mask_var,'i',('time','y','x'), zlib=True, complevel=4)
                DSnew[mask_var][:] = mask_arrays[mask_var]
                DSnew[mask_var].setncattr('units','1')

            else:
                DSnew.createVariable(mask_var,'i',('time','lat','lon'), zlib=True, complevel=4)
                DSnew[mask_var][:] = mask_arrays[mask_var]
                DSnew[mask_var].setncattr('units','1')

        DSnew.close()

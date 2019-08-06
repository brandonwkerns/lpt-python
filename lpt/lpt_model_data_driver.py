import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage

## This driver script is for historical analysis between two times specified on command line.

def filter_str(stdev):
    if type(stdev) == int:
        strout = 'g' + str(int(stdev))
    elif type(stdev) == list:
        strout = 'g' + str(int(stdev[0])) + 'x' + str(int(stdev[1]))
    else:
        print('Warning: Wrong data type!')
        strout = None
    return strout


def get_wrfout_timestamp_list(model_output_dir, domain = 1):

    DOM = str(domain).zfill(2)
    file_list = sorted(os.listdir(model_output_dir))
    time_stamp_list = [x[11:] for x in file_list if x[0:10] == ('wrfout_d' + DOM) and not 'spinup' in x]

    return time_stamp_list


def lpt_model_data_driver(dataset,plotting,output,lpo_options,lpt_options,merge_split_options,argv):

    ## Get begin and end time from command line.
    ## Give warning message if it has not been specified.

    time_stamp_list = get_wrfout_timestamp_list(dataset['raw_data_parent_dir'])
    fmt = '%Y-%m-%d_%H:00:00'
    begin_time = dt.datetime.strptime(time_stamp_list[0], fmt)
    end_time = dt.datetime.strptime(time_stamp_list[-1], fmt)
    time_list = [dt.datetime.strptime(x, fmt) for x in time_stamp_list]



    if plotting['do_plotting']:
        fig1 = plt.figure(1, figsize = (8.5,4))
        fig2 = plt.figure(2, figsize = (8.5,11))

    if lpo_options['do_lpo_calc']:

        for end_of_accumulation_time0 in time_list:

            #try:

            YMDH = end_of_accumulation_time0.strftime('%Y%m%d%H')
            YMDH_fancy = end_of_accumulation_time0.strftime('%Y-%m-%d %H:00 UTC')

            ## Get running mean rain rate.

            hours_since_init = (end_of_accumulation_time0 - begin_time).total_seconds()/3600
            if hours_since_init < 24.0:
                beginning_of_accumulation_time = begin_time
                end_of_accumulation_time = beginning_of_accumulation_time + dt.timedelta(hours=24)
            elif hours_since_init >= 24.0 and hours_since_init <= lpo_options['accumulation_hours']:
                beginning_of_accumulation_time = begin_time
                end_of_accumulation_time = end_of_accumulation_time0
            else:
                end_of_accumulation_time = end_of_accumulation_time0
                beginning_of_accumulation_time = end_of_accumulation_time0 - dt.timedelta(hours=lpo_options['accumulation_hours'])

            print((beginning_of_accumulation_time, end_of_accumulation_time))
            hours_to_divide = (end_of_accumulation_time - beginning_of_accumulation_time).total_seconds()/3600.0


            DATA_RAW1 = dataset['read_function'](end_of_accumulation_time, raw_data_parent_dir = dataset['raw_data_parent_dir'], verbose=dataset['verbose'])
            DATA_RAW0 = dataset['read_function'](beginning_of_accumulation_time, raw_data_parent_dir = dataset['raw_data_parent_dir'], verbose=dataset['verbose'])

            DATA_RUNNING = ((DATA_RAW1['precip'] - DATA_RAW0['precip']) / hours_to_divide) * 24.0 # Get the mean in mm/day.
            print('Running mean done.',flush=True)

            ## Filter the data
            DATA_FILTERED = scipy.ndimage.gaussian_filter(DATA_RUNNING, lpo_options['filter_stdev']
                , order=0, output=None, mode='reflect', cval=0.0, truncate=3.0)
            print('filter done.',flush=True)

            ## Get LP objects.
            label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, lpo_options['thresh'], verbose=dataset['verbose'])


            OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW1['lon'], DATA_RAW1['lat']
                        , DATA_RAW1['precip'], DATA_RUNNING, DATA_FILTERED, label_im, 0
                        , end_of_accumulation_time, verbose=True)
            OBJ['units_inst'] = 'mm h-1'
            OBJ['units_running'] = 'mm day-1'
            OBJ['units_filtered'] = 'mm day-1'

            print('objects properties.',flush=True)

            """
            Object Output files
            """
            objects_dir = (output['data_dir'] + '/' + dataset['label']
                            + '/' + filter_str(lpo_options['filter_stdev'])
                            + '_' + str(int(lpo_options['accumulation_hours'])) + 'h'
                            + '/thresh' + str(int(lpo_options['thresh']))
                            + '/objects/'
                            + end_of_accumulation_time.strftime(output['sub_directory_format']))

            os.makedirs(objects_dir, exist_ok = True)
            objects_fn = (objects_dir + '/objects_' + YMDH)
            lpt.lptio.lp_objects_output_ascii(objects_fn, OBJ)
            lpt.lptio.lp_objects_output_netcdf(objects_fn + '.nc', OBJ)

            """
            Object Plot
            """
            if plotting['do_plotting']:
                plt.figure(1)
                fig1.clf()
                ax1 = fig1.add_subplot(111)
                lpt.plotting.plot_rain_map_with_filtered_contour(ax1
                        , DATA_RUNNING, OBJ
                        , plot_area = plotting['plot_area'])
                ax1.set_title((dataset['label'].upper()
                                + str(lpo_options['accumulation_hours'])
                                + '-h Rain Rate and LP Objects\nEnding ' + YMDH_fancy))

                img_dir1 = (output['img_dir'] + '/' + dataset['label']
                                + '/' + filter_str(lpo_options['filter_stdev'])
                                + '_' + str(int(lpo_options['accumulation_hours'])) + 'h'
                                + '/thresh' + str(int(lpo_options['thresh']))
                                + '/objects/'
                                + end_of_accumulation_time.strftime(output['sub_directory_format']))

                os.makedirs(img_dir1, exist_ok = True)
                file_out_base = (img_dir1 + '/lp_objects_' + dataset['label'] + '_' + YMDH)
                lpt.plotting.print_and_save(file_out_base)
                fig1.clf()


    """
    LPT Tracking Calculations (i.e., connect LP Objects in time)
    """
    options = lpt_options
    options['objdir'] = (output['data_dir'] + '/' + dataset['label']
                    + '/' + filter_str(lpo_options['filter_stdev'])
                    + '_' + str(int(lpo_options['accumulation_hours'])) + 'h'
                    + '/thresh' + str(int(lpo_options['thresh'])) + '/objects')
    options['outdir'] = (output['data_dir'] + '/' + dataset['label']
                    + '/' + filter_str(lpo_options['filter_stdev'])
                    + '_' + str(int(lpo_options['accumulation_hours'])) + 'h'
                    + '/thresh' + str(int(lpo_options['thresh'])) + '/systems')
    options['imgdir'] = (output['img_dir'] + '/' + dataset['label']
                    + '/' + filter_str(lpo_options['filter_stdev'])
                    + '_' + str(int(lpo_options['accumulation_hours'])) + 'h'
                    + '/thresh' + str(int(lpo_options['thresh'])) + '/systems')
    print(options)
    print(lpo_options)

    if options['do_lpt_calc']:

        begin_tracking_time = begin_time
        latest_lp_object_time = end_time

        YMDHb = begin_time.strftime('%Y%m%d%H')
        YMDHb_fancy = begin_time.strftime('%Y-%m-%d %H:00 UTC')

        YMDH = end_time.strftime('%Y%m%d%H')
        YMDH_fancy = end_time.strftime('%Y-%m-%d %H:00 UTC')


        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(('Doing LPT tracking for: '
                + begin_tracking_time.strftime('%Y%m%d%H')
                + ' to ' + latest_lp_object_time.strftime('%Y%m%d%H')))
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

        dt_list = time_list # [begin_tracking_time + dt.timedelta(hours=x) for x in range(0, 24*lpt_options['lpt_history_days']+1, dataset['data_time_interval'])]

        ## Initialize LPT
        LPT0, BRANCHES0 = lpt.helpers.init_lpt_group_array(dt_list, options['objdir'], fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")

        ## Remove small LPOs
        LPT0, BRANCHES0 = lpt.helpers.lpt_group_array_remove_small_objects(LPT0, BRANCHES0, options, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")

        ## Connect forward, then backwards
        print('Forward...')
        LPTf, BRANCHESf = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")
        print('Backward...')
        LPTb, BRANCHESb = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, reversed=True, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")

        ## Overlap forward and backward connections.
        LPTfb, BRANCHESfb = lpt.helpers.overlap_forward_backward(LPTf, LPTb, BRANCHESf, BRANCHESb, options, verbose=True)

        ## Allow center jumps.
        print(('Allow center jumps up to ' + str(options['center_jump_max_hours']) + ' hours.'))
        LPT_center_jumps = lpt.helpers.lpt_group_array_allow_center_jumps(LPTfb, options, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")

        ## Eliminate short duration systems.
        print(('Remove LPT shorter than ' + str(options['min_lpt_duration_hours']) + ' hours.'))
        LPT_remove_short, BRANCHES_remove_short = lpt.helpers.remove_short_lived_systems(LPT_center_jumps, BRANCHESfb, options['min_lpt_duration_hours']
                                , latest_datetime = latest_lp_object_time)




        ## Handle splitting and merging, if specified.
        if merge_split_options['allow_merge_split']:
            LPT, BRANCHES = lpt.helpers.lpt_group_id_separate_branches(LPT_remove_short, BRANCHES_remove_short, options, verbose=True, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")
            LPT, BRANCHES = lpt.helpers.lpt_split_and_merge(LPT, BRANCHES, merge_split_options)
        else:
            LPT = LPT_remove_short.copy()
            BRANCHES = BRANCHES_remove_short.copy()

        BRANCHES = lpt.helpers.reorder_LPT_branches(LPT, BRANCHES)

        ## Get "timeclusters" tracks.
        print('Calculating LPT properties.')
        if merge_split_options['allow_merge_split']:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPT, BRANCHES, options, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")
        else:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties(LPT, options, fmt = "/"+output['sub_directory_format'] + "/objects_%Y%m%d%H.nc")

        fn_tc_base = (options['outdir'] + '/' + end_time.strftime(output['sub_directory_format'])
                         + '/lpt_systems_' + dataset['label'] + '_' + YMDHb + '_' + YMDH)
        lpt.lptio.lpt_system_tracks_output_ascii(fn_tc_base + '.txt', TIMECLUSTERS)
        lpt.lptio.lpt_systems_group_array_output_ascii(fn_tc_base + '.group_array.txt', LPT, BRANCHES)
        lpt.lptio.lpt_system_tracks_output_netcdf(fn_tc_base + '.nc', TIMECLUSTERS)


        """
        LPT Plotting
        """

        if plotting['do_plotting']:
            plt.figure(2)
            fig2.clf()
            ax2 = fig2.add_subplot(111)

            timelon_rain = []
            for this_dt in dt_list:
                if (this_dt - begin_time).total_seconds()/3600 > 1.01:

                    DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'], raw_data_parent_dir = dataset['raw_data_parent_dir'])
                    DATA_RAW0 = dataset['read_function'](this_dt-dt.timedelta(hours=1), verbose=dataset['verbose'], raw_data_parent_dir = dataset['raw_data_parent_dir'])
                    DATA_RAW['precip'] -= DATA_RAW0['precip']

                else:

                    DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'], raw_data_parent_dir = dataset['raw_data_parent_dir'])

                #lat_idx, = np.where(np.logical_and(DATA_RAW['lat'] > -15.0, DATA_RAW['lat'] < 15.0))
                #timelon_rain.append(np.mean(np.array(DATA_RAW['precip'][lat_idx,:]), axis=0))
                this_row = []
                timelon_lon = []
                for ii in range(len(DATA_RAW['lat'][0,:])):
                    this_lat = DATA_RAW['lat'][:,ii]
                    this_lon = DATA_RAW['lon'][:,ii]
                    this_precip = DATA_RAW['precip'][:,ii]
                    this_row.append(np.nanmean(this_precip[np.abs(this_lat) < 15.0]))
                    timelon_lon.append(np.nanmean(this_lon[np.abs(this_lat) < 15.0]))
                timelon_rain.append(this_row)

            timelon_lon = np.array(timelon_lon)
            timelon_rain = np.array(timelon_rain)

            lpt.plotting.plot_timelon_with_lpt(ax2, dt_list, timelon_lon
                    , timelon_rain, TIMECLUSTERS, plotting['time_lon_range']
                    , accum_time_hours = lpo_options['accumulation_hours'])

            ax2.set_title((dataset['label'].upper()
                            + ' Rain Rate (15$\degree$S-15$\degree$N) and LPTs\n' + YMDHb_fancy + ' to ' + YMDH_fancy))

            ax2.text(0.87,1.02,'(<15$\degree$S, >15$\degree$N Dashed)', transform=ax2.transAxes)

            img_dir2 = (options['imgdir'] + end_time.strftime(output['sub_directory_format']))

            os.makedirs(img_dir2, exist_ok = True)
            file_out_base = (img_dir2 + '/lpt_time_lon_' + dataset['label'] + '_' + YMDHb + '_' + YMDH)
            lpt.plotting.print_and_save(file_out_base)
            fig2.clf()

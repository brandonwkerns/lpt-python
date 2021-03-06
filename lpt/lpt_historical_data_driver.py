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


def lpt_historical_data_driver(dataset,plotting,output,lpo_options,lpt_options,merge_split_options,argv):

    ## Get begin and end time from command line.
    ## Give warning message if it has not been specified.

    if len(argv) < 3:
        print('Specify begin and end time on command line, format YYYYMMDDHH.')
        return

    begin_time = dt.datetime.strptime(str(argv[1]), '%Y%m%d%H') # command line arg #1 format: YYYYMMDDHH
    end_time = dt.datetime.strptime(str(argv[2]), '%Y%m%d%H') # command line arg #1 format: YYYYMMDDHH

    hours_list = np.arange(0.0, 0.1 +(end_time-begin_time).total_seconds()/3600.0, dataset['data_time_interval'])
    time_list = [begin_time + dt.timedelta(hours=x) for x in hours_list]



    if plotting['do_plotting']:
        fig1 = plt.figure(1, figsize = (8.5,4))
        fig2 = plt.figure(2, figsize = (8.5,11))

    if lpo_options['do_lpo_calc']:

        for end_of_accumulation_time in time_list:

            try:

                YMDH = end_of_accumulation_time.strftime('%Y%m%d%H')
                YMDH_fancy = end_of_accumulation_time.strftime('%Y-%m-%d %H:00 UTC')

                beginning_of_accumulation_time = end_of_accumulation_time - dt.timedelta(hours=lpo_options['accumulation_hours'])
                print((beginning_of_accumulation_time, end_of_accumulation_time))
                if lpo_options['accumulation_hours'] > 0.1:
                    dt_list = [beginning_of_accumulation_time
                        + dt.timedelta(hours=x) for x in np.double(np.arange(0,lpo_options['accumulation_hours']
                                                          + dataset['data_time_interval'],dataset['data_time_interval']))]
                else:
                    dt_list = [end_of_accumulation_time]

                ## Get accumulated rain.
                data_collect = []
                count = 0

                for this_dt in reversed(dt_list):
                    if 'sub_area' in dataset.keys():
                        DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'], area=dataset['sub_area'])
                    else:
                        DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'])
                    DATA_RAW['precip'][DATA_RAW['precip'] < -0.01] = 0.0
                    if count < 1:
                        data_collect = np.array(DATA_RAW['precip'])
                    else:
                        data_collect += np.array(DATA_RAW['precip'])
                    count += 1

                DATA_RUNNING = (data_collect/count) * 24.0 # Get the mean in mm/day.
                print('Running mean done.',flush=True)

                fig=plt.figure()
                plt.pcolormesh(DATA_RUNNING)
                plt.colorbar()
                plt.savefig('test0.png')
                plt.close(fig)

                ## Filter the data
                if lpo_options['filter_stdev'] > 0.001:
                    DATA_FILTERED = scipy.ndimage.gaussian_filter(DATA_RUNNING, lpo_options['filter_stdev']
                        , order=0, output=None, mode='reflect', cval=0.0, truncate=3.0)
                    print('filter done.',flush=True)
                else:
                    DATA_FILTERED = DATA_RUNNING.copy()
                    print('Filter set to 0, so not doing filter.')

                ## Get LP objects.
                object_is_gt_threshold = True
                if 'object_is_gt_threshold' in lpo_options:
                    if not lpo_options['object_is_gt_threshold']:
                        object_is_gt_threshold = False

                if not 'min_points' in lpo_options:
                    lpo_options['min_points'] = 1
                label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, lpo_options['thresh'], object_is_gt_threshold=object_is_gt_threshold, min_points=lpo_options['min_points'],verbose=dataset['verbose'])
                OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW['lon'], DATA_RAW['lat']
                            , DATA_RAW['precip'], DATA_RUNNING, DATA_FILTERED, label_im, 0
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
                #if (len(OBJ['n_points']) > 0):
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

            except FileNotFoundError:
                print('Data not yet available up to this point. Skipping.')


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
        LPT0, BRANCHES0 = lpt.helpers.init_lpt_group_array(dt_list, options['objdir'])

        ## Remove small LPOs
        LPT0, BRANCHES0 = lpt.helpers.lpt_group_array_remove_small_objects(LPT0, BRANCHES0, options)

        ## Connect forward, then backwards
        print('Forward...')
        LPTf, BRANCHESf = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True)
        print('Backward...')
        LPTb, BRANCHESb = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, reversed=True)

        ## Overlap forward and backward connections.
        LPTfb, BRANCHESfb = lpt.helpers.overlap_forward_backward(LPTf, LPTb, BRANCHESf, BRANCHESb, options, verbose=True)

        ## Allow center jumps.
        print(('Allow center jumps up to ' + str(options['center_jump_max_hours']) + ' hours.'))
        LPT_center_jumps = lpt.helpers.lpt_group_array_allow_center_jumps(LPTfb, options)

        ## Eliminate short duration systems.
        print(('Remove LPT shorter than ' + str(options['min_lpt_duration_hours']) + ' hours.'))
        LPT_remove_short, BRANCHES_remove_short = lpt.helpers.remove_short_lived_systems(LPT_center_jumps, BRANCHESfb, options['min_lpt_duration_hours']
                                , latest_datetime = latest_lp_object_time)




        ## Handle splitting and merging, if specified.
        if merge_split_options['allow_merge_split']:
            LPT, BRANCHES = lpt.helpers.lpt_group_id_separate_branches(LPT_remove_short, BRANCHES_remove_short, options, verbose=True)
            LPT, BRANCHES = lpt.helpers.lpt_split_and_merge(LPT, BRANCHES, merge_split_options, options)
        else:
            LPT = LPT_remove_short.copy()
            BRANCHES = BRANCHES_remove_short.copy()


        BRANCHES = lpt.helpers.reorder_LPT_branches(LPT, BRANCHES)

        ## Get "timeclusters" tracks.
        print('Calculating LPT properties.')
        if merge_split_options['allow_merge_split']:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPT, BRANCHES, options)
        else:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties(LPT, options)

        fn_tc_base = (options['outdir'] #+ '/' + end_time.strftime(output['sub_directory_format'])
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
                if 'sub_area' in dataset.keys():
                    DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'], area=dataset['sub_area'])
                else:
                    DATA_RAW = dataset['read_function'](this_dt, verbose=dataset['verbose'])

                lat_idx, = np.where(np.logical_and(DATA_RAW['lat'] > -15.0, DATA_RAW['lat'] < 15.0))
                timelon_rain.append(np.mean(np.array(DATA_RAW['precip'][lat_idx,:]), axis=0))


            lpt.plotting.plot_timelon_with_lpt(ax2, dt_list, DATA_RAW['lon']
                    , timelon_rain, TIMECLUSTERS, plotting['time_lon_range']
                    , accum_time_hours = lpo_options['accumulation_hours'])

            ax2.set_title((dataset['label'].upper()
                            + ' Rain Rate (15$\degree$S-15$\degree$N) and LPTs\n' + YMDHb_fancy + ' to ' + YMDH_fancy))

            ax2.text(0.87,1.02,'(<15$\degree$S, >15$\degree$N Dashed)', transform=ax2.transAxes)

            img_dir2 = (output['img_dir'] + '/' + dataset['label'] + '/systems/')
            #                + end_time.strftime(output['sub_directory_format']))

            os.makedirs(img_dir2, exist_ok = True)
            file_out_base = (img_dir2 + '/lpt_time_lon_' + dataset['label'] + '_' + YMDHb + '_' + YMDH)
            lpt.plotting.print_and_save(file_out_base)
            fig2.clf()

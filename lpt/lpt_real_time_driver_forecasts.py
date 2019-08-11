import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pyplot as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage

## This driver script is for forecasts, e.g., going forward in time.

def lpt_real_time_driver_forecasts(dataset,plotting,output,lpo_options,lpt_options,merge_split_options,argv):

    ## Use current real time, or specified time from the command line args.
    if len(argv) < 2:
        base_time0 = dt.datetime.utcnow()
        year0, month0, day0, hour0 = base_time0.timetuple()[0:4]  #Note: this is the initialization time.
        base_time = dt.datetime(year0, month0, day0, 0, 0, 0)   # 00Z of same day.

    else:
        base_time = dt.datetime.strptime(str(argv[1]), '%Y%m%d%H') # command line arg #1 format: YYYYMMDDHH

    year, month, day, hour = base_time.timetuple()[0:4]  #Note: this is the initialization time.
    hour = 0 #int(np.floor(hour/dataset['data_time_interval']) * dataset['data_time_interval'])

    if plotting['do_plotting']:
        fig1 = plt.figure(1, figsize = (8.5,4))
        fig2 = plt.figure(2, figsize = (8.5,11))

    model_init_time = dt.datetime(year, month, day, hour, 0, 0)


    try:
        if 'sub_area' in dataset.keys():
            DATA_RAW = dataset['read_function'](model_init_time, verbose=dataset['verbose'], area=dataset['sub_area'])
        else:
            DATA_RAW = dataset['read_function'](model_init_time, verbose=dataset['verbose'])

        for fcst_hours in range(lpo_options['accumulation_hours'], dataset['max_forecast_hours']+dataset['forecast_hours_interval'], dataset['forecast_hours_interval']):

            record_num = int(fcst_hours/dataset['forecast_hours_interval'])-1
            record_num0 = int((fcst_hours-lpo_options['accumulation_hours'])/dataset['forecast_hours_interval'])

            end_of_accumulation_time = base_time + dt.timedelta(hours=fcst_hours)

            YMDH = end_of_accumulation_time.strftime('%Y%m%d%H')
            YMDH_fancy = end_of_accumulation_time.strftime('%Y-%m-%d %H:00 UTC')
            YMDHi = model_init_time.strftime('%Y%m%d%H')
            YMDHi_fancy = model_init_time.strftime('%Y-%m-%d %H:00 UTC')

            beginning_of_accumulation_time = end_of_accumulation_time - dt.timedelta(hours=lpo_options['accumulation_hours'])
            print((beginning_of_accumulation_time, end_of_accumulation_time))
            dt_list = [beginning_of_accumulation_time
                + dt.timedelta(hours=x) for x in np.double(np.arange(0,lpo_options['accumulation_hours']
                                                  + dataset['forecast_hours_interval'],dataset['forecast_hours_interval']))]

            record_list = range(record_num0, record_num+1)

            ## Get accumulated rain.
            data_collect = []
            count = 0

            for this_record in record_list:
                if count < 1:
                    data_collect = np.array(DATA_RAW['precip'][this_record,:,:])
                else:
                    data_collect += np.array(DATA_RAW['precip'][this_record,:,:])
                count += 1

            DATA_RUNNING = (data_collect/count) * 24.0 # Get the mean in mm/day.
            print('Accum done.',flush=True)

            ## Filter the data
            DATA_FILTERED = scipy.ndimage.gaussian_filter(DATA_RUNNING, lpo_options['filter_stdev']
                , order=0, output=None, mode='reflect', cval=0.0, truncate=3.0)
            print('filter done.',flush=True)

            ## Get LP objects.
            label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, lpo_options['thresh'], verbose=dataset['verbose'])
            OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW['lon'], DATA_RAW['lat']
                        , DATA_RAW['precip'], DATA_RUNNING, DATA_FILTERED, label_im, 0
                        , end_of_accumulation_time, verbose=True)
            print('objects properties.',flush=True)

            """
            Object Output files
            """
            objects_dir = (output['data_dir'] + '/' + dataset['label'] + '/objects/'
                            + model_init_time.strftime(output['sub_directory_format']))
            os.makedirs(objects_dir, exist_ok = True)
            objects_fn = (objects_dir + '/objects_' + YMDH)
            lpt.lptio.lp_objects_output_ascii(objects_fn, OBJ)
            if (len(OBJ['n_points']) > 0):
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
                ax1.set_title((dataset['label'].upper() + ' FCST '
                                + str(lpo_options['accumulation_hours'])
                                + '-h Rain Rate and LP Objects'
                                + '\n             Init: ' + YMDHi_fancy
                                + '\nValid: ending ' + YMDH_fancy))

                img_dir1 = (output['img_dir'] + '/' + dataset['label'] + '/objects/'
                                + model_init_time.strftime(output['sub_directory_format']))

                os.makedirs(img_dir1, exist_ok = True)
                file_out_base = (img_dir1 + '/lp_objects_' + dataset['label'] + '_rt_' + YMDH)
                lpt.plotting.print_and_save(file_out_base)
                fig1.clf()

    except FileNotFoundError:
        print('Data not yet available up to this point. Skipping.')



    """
    LPT Tracking Calculations (i.e., connect LP Objects in time)
    """
    options = lpt_options
    options['objdir'] = (output['data_dir'] + '/' + dataset['label'] + '/objects/'
        + model_init_time.strftime(output['sub_directory_format']))
    options['outdir'] = (output['data_dir'] + '/' + dataset['label'] + '/systems/'
        + model_init_time.strftime(output['sub_directory_format']))

    if options['do_lpt_calc']:

        latest_lp_object_time = lpt.helpers.get_latest_lp_object_time(options['objdir'], level=0)
        begin_tracking_time = model_init_time

        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(('Doing LPT tracking for: '
                + begin_tracking_time.strftime('%Y%m%d%H')
                + ' to ' + latest_lp_object_time.strftime('%Y%m%d%H')))
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

        #dt_list = [begin_tracking_time + dt.timedelta(hours=x) for x in range(0, 24*lpt_options['lpt_history_days']+1, dataset['data_time_interval'])]
        dt_list = [begin_tracking_time + dt.timedelta(hours=x) for x in range(lpo_options['accumulation_hours'], dataset['max_forecast_hours']+dataset['forecast_hours_interval'], dataset['forecast_hours_interval'])]
        print(dt_list)

        ## Initialize LPT
        #fmt = model_init_time.strftime("/%Y/%m/%Y%m%d") + "/objects_%Y%m%d%H.nc"
        fmt = "/objects_%Y%m%d%H.nc"
        LPT0, BRANCHES0 = lpt.helpers.init_lpt_group_array(dt_list, options['objdir'], fmt)
        ## Remove small LPOs
        LPT0, BRANCHES0 = lpt.helpers.lpt_group_array_remove_small_objects(LPT0, BRANCHES0, options, fmt=fmt)

        ## Connect forward, then backwards
        print('Forward...')
        LPTf, BRANCHESf = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, fmt=fmt)
        print('Backward...')
        LPTb, BRANCHESb = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, reversed=True, fmt=fmt)

        ## Overlap forward and backward connections.
        LPTfb, BRANCHESfb = lpt.helpers.overlap_forward_backward(LPTf, LPTb, BRANCHESf, BRANCHESb, options, verbose=True)

        ## Allow center jumps.
        print(('Allow center jumps up to ' + str(options['center_jump_max_hours']) + ' hours.'))
        LPT_center_jumps = lpt.helpers.lpt_group_array_allow_center_jumps(LPTfb, options, fmt=fmt)

        ## Eliminate short duration systems.
        print(('Remove LPT shorter than ' + str(options['min_lpt_duration_hours']) + ' hours.'))
        LPT_remove_short, BRANCHES_remove_short = lpt.helpers.remove_short_lived_systems(LPT_center_jumps, BRANCHESfb, options['min_lpt_duration_hours']
                                , latest_datetime = latest_lp_object_time - dt.timedelta(hours=options['min_lpt_duration_hours']))



        ## Handle splitting and merging, if specified.
        if merge_split_options['allow_merge_split']:
            LPT, BRANCHES = lpt.helpers.lpt_group_id_separate_branches(LPT_remove_short, BRANCHES_remove_short, options, verbose=True, fmt=fmt)
            LPT, BRANCHES = lpt.helpers.lpt_split_and_merge(LPT, BRANCHES, merge_split_options, options)
        else:
            LPT = LPT_remove_short.copy()
            BRANCHES = BRANCHES_remove_short.copy()

        ## Get "timeclusters" tracks.
        print('Calculating LPT properties.')
        if merge_split_options['allow_merge_split']:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPT, BRANCHES, options, fmt)
        else:
            TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties(LPT, options, fmt)


        ## Get "timeclusters" tracks.
        #print('Calculating LPT properties.')
        #TIMECLUSTERS = lpt.helpers.calc_lpt_system_group_properties(LPT_remove_short, options, fmt=fmt)

        fn_tc_base = (options['outdir'] #+ '/' + end_of_accumulation_time.strftime(output['sub_directory_format'])
                         + '/lpt_systems_' + dataset['label'] + '_init' + YMDHi + '__' + str(options['lpt_history_days']) + 'days')
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
        if 'sub_area' in dataset.keys():
            DATA_RAW = dataset['read_function'](model_init_time, verbose=dataset['verbose'], area=dataset['sub_area'])
        else:
            DATA_RAW = dataset['read_function'](model_init_time, verbose=dataset['verbose'])

        lat_avg_list = [x for x in range(len(DATA_RAW['lat'])) if abs(DATA_RAW['lat'][x]) < 15.0]
        timelon_rain = np.mean(DATA_RAW['precip'][:,lat_avg_list,:], axis=1)
        all_model_time = [begin_tracking_time + dt.timedelta(hours=x) for x in range(dataset['forecast_hours_interval'], dataset['max_forecast_hours']+dataset['forecast_hours_interval'], dataset['forecast_hours_interval'])]

        lpt.plotting.plot_timelon_with_lpt(ax2, all_model_time, DATA_RAW['lon']
                , timelon_rain, TIMECLUSTERS, plotting['time_lon_range']
                , accum_time_hours = lpo_options['accumulation_hours'], label_font_size=plotting['label_font_size'])

        ax2.set_title((dataset['label'].upper() + ' ' + str(lpt_options['lpt_history_days']) + ' Day'
                        + ' FCST 15S-15N Rain Rate and LPTs\nInitialized: ' + YMDHi_fancy), fontsize=plotting['label_font_size'])

        img_dir2 = (output['img_dir'] + '/' + dataset['label'] + '/systems/'
                        + model_init_time.strftime(output['sub_directory_format']))

        os.makedirs(img_dir2, exist_ok = True)
        file_out_base = (img_dir2 + '/lpt_time_lon_' + dataset['label'] + '_rt_init' + YMDHi
                            + '__' + str(options['lpt_history_days']) + 'days')
        lpt.plotting.print_and_save(file_out_base)
        fig2.clf()

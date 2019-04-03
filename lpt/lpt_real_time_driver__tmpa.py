import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage

plt.close('all')
plt.ioff()

"""
Main settings for lpt
"""
THRESH=12.0
accumulation_hours = 72
data_time_interval = 3 #Time resolution of the data in hours.
filter_stdev = 20 # in terms of number of grid points.

plot_area = [50, 200, -30, 30]
img_dir = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images'
data_dir = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'

################################################################################

## Use current real time, or specified time from the command line args.
if len(sys.argv) < 2:
    base_time = dt.datetime.utcnow()
else:
    base_time = dt.datetime.strptime(str(sys.argv[1]), '%Y%m%d%H') # command line arg #1 format: YYYYMMDDHH

if len(sys.argv) < 3:
    hours_to_go_back = 24
else:
    hours_to_go_back = int(sys.argv[2]) # command line arg #2: hours to go back.

year, month, day, hour = base_time.timetuple()[0:4]
hour = int(np.floor(hour/3) * 3)

current_end_of_accumulation_time = dt.datetime(year,month,day,hour,0,0)

fig = plt.figure(figsize=(8.5,4))

## Check back 24 h from current time.
for hours_back in range(0, hours_to_go_back+1, data_time_interval):

    try:
        end_of_accumulation_time = current_end_of_accumulation_time - dt.timedelta(hours=hours_back)

        YMDH=end_of_accumulation_time.strftime('%Y%m%d%H')
        YMDH_fancy=end_of_accumulation_time.strftime('%Y-%m-%d %H:00 UTC')

        beginning_of_accumulation_time = end_of_accumulation_time - dt.timedelta(hours=accumulation_hours)
        print((beginning_of_accumulation_time, end_of_accumulation_time))
        dt_list = [beginning_of_accumulation_time
            + dt.timedelta(hours=x) for x in np.double(np.arange(0,accumulation_hours
                                              + data_time_interval,data_time_interval))]

        ## Get accumulated rain.
        data_collect = []
        count = 0
        for this_dt in reversed(dt_list):
            DATA_RAW=lpt.readdata.read_tmpa_at_datetime(this_dt, verbose=True)
            if count < 1:
                data_collect = np.array(DATA_RAW['precip'])
            else:
                data_collect += np.array(DATA_RAW['precip'])
            count += 1

        DATA_ACCUM = (data_collect/count) * 24.0 # Get the mean in mm/day.
        print('Accum done.',flush=True)

        ## Filter the data
        DATA_FILTERED = scipy.ndimage.gaussian_filter(DATA_ACCUM, filter_stdev
            , order=0, output=None, mode='reflect', cval=0.0, truncate=3.0)
        print('filter done.',flush=True)

        ## Get LP objects.
        label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, THRESH, verbose=True)
        OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW['lon'], DATA_RAW['lat']
                    , DATA_RAW['precip'], DATA_ACCUM, label_im, 0
                    , end_of_accumulation_time, verbose=True)
        print('objects properties.',flush=True)

        ## Output files
        objects_dir = (data_dir + '/tmpa/objects/' + str(end_of_accumulation_time.year)
         + '/' + str(end_of_accumulation_time.month).zfill(2)
         + '/' + end_of_accumulation_time.strftime('%Y%m%d'))
        os.makedirs(objects_dir, exist_ok = True)
        objects_fn = (objects_dir + '/objects_' + str(end_of_accumulation_time.year)
         + str(end_of_accumulation_time.month).zfill(2)
         + str(end_of_accumulation_time.day).zfill(2)
         + str(end_of_accumulation_time.hour).zfill(2))
        lpt.lptio.lp_objects_output_ascii(objects_fn, OBJ)
        if (len(OBJ['n_points']) > 0):
            lpt.lptio.lp_objects_output_netcdf(objects_fn + '.nc', OBJ)

        ## Plot
        plt.clf()
        ax1 = fig.add_subplot(111)
        map1=lpt.plotting.plot_map_background(plot_area)
        cmap = lpt.plotting.cmap_map(lambda x: x/2 + 0.5, plt.cm.jet)
        cmap.set_under(color='white')
        H1 = map1.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'],DATA_ACCUM, cmap=cmap, vmin=1, vmax=50)
        H2 = plt.contour(DATA_RAW['lon'], DATA_RAW['lat'],DATA_FILTERED, [THRESH,], colors='k', linewidths=1.0)

        map1.plot(OBJ['lon'], OBJ['lat'], 'kx', markersize=7)
        plt.colorbar(H1)

        ax1.set_title('TMPA RT 3-Day Rain Rate and LP Objects\n' + YMDH_fancy)

        img_dir2 = (img_dir + '/tmpa/objects/' + str(end_of_accumulation_time.year)
                    + '/' + str(end_of_accumulation_time.month).zfill(2)
                    + '/' + end_of_accumulation_time.strftime('%Y%m%d'))
        os.makedirs(img_dir2, exist_ok = True)
        file_out_base = (img_dir2 + '/lp_objects_tmpa_rt_' + YMDH)

        lpt.plotting.print_and_save(file_out_base)
        plt.clf()

    except FileNotFoundError:
        print('Data not yet available up to this point. Skipping.')

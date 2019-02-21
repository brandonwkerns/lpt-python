import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import os
import matplotlib.colors as colors
import scipy.ndimage

plt.close('all')

"""
Main settings for lpt
"""
THRESH=14.0
accumulation_hours = 72
data_time_interval = 1
filter_stdev = 70 # in terms of number of grid points.

plot_area = [50, 200, -30, 30]
img_dir = '/home/orca/bkerns/public_html/realtime_mjo_tracking/images/lpt'
data_dir = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'

################################################################################

def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return colors.LinearSegmentedColormap('colormap',cdict,1024)






year, month, day, hour = dt.datetime.utcnow().timetuple()[0:4]
hour = int(np.floor(hour/3) * 3)

current_end_of_accumulation_time = dt.datetime(year,month,day,hour,0,0)

## Check back 24 h from current time.
for hours_back in range(0,13,data_time_interval):

    try:
        end_of_accumulation_time = current_end_of_accumulation_time - dt.timedelta(hours=hours_back)

        YMDH=end_of_accumulation_time.strftime('%Y%m%d%H')
        YMDH_fancy=end_of_accumulation_time.strftime('%Y-%m-%d %H:00 UTC')

        dt_list = [end_of_accumulation_time - dt.timedelta(hours=accumulation_hours)
                  + dt.timedelta(hours=x) for x in np.double(np.arange(0,accumulation_hours
                                              + 1, 1))]

        ## Get accumulated rain.
        data_collect = []
        for this_dt in reversed(dt_list):
            DATA_RAW=lpt.readdata.read_cmorph_at_datetime(this_dt, area=[40,210,-40,40], verbose=True)
            data_collect.append(DATA_RAW['precip'][0,:,:]+DATA_RAW['precip'][1,:,:])
        data_collect3d = np.array(data_collect)
        DATA_ACCUM = np.nanmean(data_collect3d,axis=0) * 24.0
        print('Accum done.',flush=True)
        ## Filter the data
        #DATA_FILTERED = lpt.helpers.gauss_smooth(DATA_ACCUM, filter_stdev)
        DATA_FILTERED = scipy.ndimage.gaussian_filter(DATA_ACCUM, filter_stdev
            , order=0, output=None, mode='reflect', cval=0.0, truncate=3.0)
        print('filter done.',flush=True)
        label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, THRESH, verbose=True)
        print('objects id.',flush=True)
        OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW['lon'], DATA_RAW['lat']
                    , DATA_RAW['precip'], DATA_ACCUM, label_im, verbose=True)
        print('objects properties.',flush=True)

        fig = plt.figure(figsize=(8.5,4))

        ax1 = fig.add_subplot(111)
        map1=lpt.helpers.plot_map_background(plot_area)
        cmap = cmap_map(lambda x: x/2 + 0.5, plt.cm.jet)
        cmap.set_under(color='white')
        H1 = map1.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'],DATA_ACCUM, cmap=cmap
                            , vmin=1, vmax=50)
        H2 = plt.contour(DATA_RAW['lon'], DATA_RAW['lat'],DATA_FILTERED, [THRESH,], colors='k', linewidths=1.0)

        map1.plot(OBJ['lon'], OBJ['lat'], 'kx', markersize=7)
        plt.colorbar(H1)

        ax1.set_title('CMORPH RT 3-Day Rain Rate and LP Objects\n' + YMDH_fancy)

        file_out_base = (img_dir + '/lp_objects_cmorph_rt_' + YMDH)

        lpt.helpers.print_and_save(file_out_base)

    except FileNotFoundError:
        print('Data not yet available up to this point. Skipping.')

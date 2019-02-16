import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import os

plt.close('all')

kernel = lpt.helpers.gauss_smooth_kernel(121,121,20,20)

THRESH=12.0
accumulation_hours = 72
data_time_interval = 3
end_of_accumulation_time = dt.datetime(2011,11,24,0,0,0)
#end_of_accumulation_time = dt.datetime(2019,2,15,0,0,0)
filter_stdev = 20 # in terms of number of grid points.

dt_list = [end_of_accumulation_time - dt.timedelta(hours=accumulation_hours)
          + dt.timedelta(hours=x) for x in np.double(np.arange(0,accumulation_hours
                                      + data_time_interval,data_time_interval))]
print(dt_list)

## Get accumulated rain.
data_collect = []
for this_dt in dt_list:
    DATA_RAW=lpt.readdata.read_tmpa_at_datetime(this_dt, verbose=True)
    data_collect.append(DATA_RAW['precip'])
data_collect3d = np.array(data_collect)
DATA_ACCUM = np.nanmean(data_collect3d,axis=0) * 24.0

## Filter the data
DATA_FILTERED = lpt.helpers.gauss_smooth(DATA_ACCUM, filter_stdev)

## Show plot
"""
fig = plt.figure(figsize=(10,6))

cmap = 'jet'
ax1 = fig.add_subplot(221)
ax1.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'], 24.0*data_collect3d[0,:,:], cmap=cmap, vmin=0, vmax=20)

ax2 = fig.add_subplot(222)
ax2.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'], 24.0*data_collect3d[-1,:,:], cmap=cmap, vmin=0, vmax=20)

ax3 = fig.add_subplot(223)
ax3.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'], DATA_ACCUM, cmap=cmap, vmin=0, vmax=20)

ax4 = fig.add_subplot(224)
ax4.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'], DATA_FILTERED, cmap=cmap, vmin=0, vmax=20)
ax4.contour(DATA_RAW['lon'], DATA_RAW['lat'], DATA_FILTERED, [12,], colors='k')
"""

label_im = lpt.helpers.identify_lp_objects(DATA_FILTERED, THRESH, verbose=True)

OBJ = lpt.helpers.calculate_lp_object_properties(DATA_RAW['lon'], DATA_RAW['lat']
            , DATA_RAW['precip'], DATA_ACCUM, label_im, verbose=True)


fig2 = plt.figure(figsize=(6,4))
H=plt.pcolormesh(DATA_RAW['lon'], DATA_RAW['lat'],label_im, cmap='jet')
plt.plot(OBJ['lon'], OBJ['lat'], 'kx')
plt.colorbar(H)
plt.show()

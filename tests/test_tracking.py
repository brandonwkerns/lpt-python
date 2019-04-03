import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt

plt.close('all')


options={}
options['min_overlap_points'] = 4000
options['min_overlap_frac'] = 0.5
options['min_lp_objects_points'] = 100


## For TMPA
options['objdir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/tmpa/objects'
#dt_list = [dt.datetime(2019,2,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,61*24.0+1,3.0)]
dt_list = [dt.datetime(2019,3,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,15*24.0+1,3.0)]


## For CMORPH
#options['objdir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/cmorph/objects'
#dt_list = [dt.datetime(2019,3,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,31*24.0+1,1.0)]


print('Initialization.')
print(dt.datetime.now(), flush=True)
LPT = lpt.helpers.init_lpt_group_array(dt_list, options['objdir'])
print(dt.datetime.now(), flush=True)

print('Remove small LPOs.')
print(dt.datetime.now(), flush=True)
LPT0 = lpt.helpers.lpt_group_array_remove_small_objects(LPT.copy(), options)
print(dt.datetime.now(), flush=True)

print('Looping forward')
print(dt.datetime.now(), flush=True)
LPTf = lpt.helpers.calc_lpt_group_array(LPT0.copy(), options, verbose=True)
print(dt.datetime.now(), flush=True)

fig0 = plt.figure(figsize=(11.0,8.5))
ax0 = fig0.add_subplot(131)
#lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax0, LPTf, options)
ax0.set_title('Forward')


print('Looping backwards')
print(dt.datetime.now(), flush=True)
LPTb = lpt.helpers.calc_lpt_group_array(LPT0.copy(), options, verbose=True, reversed=True)
print(dt.datetime.now(), flush=True)

#fig1 = plt.figure(figsize=(8.5,11))
ax1 = fig0.add_subplot(132)
#lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax1, LPTb, options)
ax1.set_title('Backwards')

print('Merging forward and backwards')
print(dt.datetime.now(), flush=True)
LPTfb = lpt.helpers.overlap_forward_backward(LPTf.copy(), LPTb.copy(), options, verbose=True)
print(dt.datetime.now(), flush=True)


print('Calc LP system group bulk properties.')
print(dt.datetime.now(), flush=True)
TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties(LPTfb, options)
print(dt.datetime.now(), flush=True)

#print('Calc LP system group bulk properties with branches.')
#LPTbranches = lpt.helpers.separate_lpt_system_branches(LPTfb, LPTf, LPTb, options)

#fig = plt.figure(figsize=(8.5,11))
ax = fig0.add_subplot(133)
lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax, LPTfb, options, text_color='darkgrey')
ax.set_title('Forward and Backwards, No Branches')

plt.savefig('test_tracking.png')

fn_tc = 'TEST.txt'
lpt.lptio.lpt_system_tracks_output_ascii(fn_tc, TIMECLUSTERS0)

plt.show()

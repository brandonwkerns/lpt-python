import numpy as np
from context import lpt
import matplotlib.pylab as plt
from matplotlib import dates
import datetime as dt

plt.close('all')


options={}
options['min_overlap_points'] = 4000
options['min_overlap_frac'] = 0.5
options['min_lp_objects_points'] = 100
options['center_jump_max_hours'] = 3*24   # How long to allow center jumps
options['min_lpt_duration_hours'] = 7*24  # Minumum duration to keep it as an LPT


## For TMPA
options['objdir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/tmpa/objects'
#dt_list = [dt.datetime(2019,2,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,61*24.0+1,3.0)]
dt_list = [dt.datetime(2019,3,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,16*24.0+1,3.0)]
#dt_list = [dt.datetime(2019,3,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,4*24.0+1,3.0)]


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

fig0 = plt.figure(figsize=(50.0,30))
ax0 = fig0.add_subplot(131)
lpt.helpers.plot_lpt_groups_time_lon_text(ax0, LPTf, options)
ax0.set_title('Forward')


print('Looping backwards')
print(dt.datetime.now(), flush=True)
LPTb = lpt.helpers.calc_lpt_group_array(LPT0.copy(), options, verbose=True, reversed=True)
print(dt.datetime.now(), flush=True)

ax1 = fig0.add_subplot(132)
lpt.helpers.plot_lpt_groups_time_lon_text(ax1, LPTb, options)
ax1.set_title('Backwards')

print('Merging forward and backwards')
LPTfb0 = lpt.helpers.overlap_forward_backward(LPTf.copy(), LPTb.copy(), options, verbose=True)

## Allow center jumps.
print(('Allow center jumps up to ' + str(options['center_jump_max_hours']) + ' hours.'))
LPTfb1 = lpt.helpers.lpt_group_array_allow_center_jumps(LPTfb0, options)

print(('Remove LPT shorter than ' + str(options['min_lpt_duration_hours']) + ' hours.'))
LPTfb2 = lpt.helpers.remove_short_lived_systems(LPTfb1, options['min_lpt_duration_hours'])

LPTfb = lpt.helpers.lpt_group_id_separate_branches(LPTfb2, options, verbose=True)

print('Calc LP system group bulk properties.')
TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPTfb, options)


## Re-combine when splits/mergers are < 3 days.
lpt_diff = lpt.helpers.lpt_branches_difference(LPTfb, 1, 1, 2)
print('Diff')
print(str(lpt_diff))
print(str(LPTfb[lpt_diff,1].astype('int')))

lpt_intersection = lpt.helpers.lpt_branches_intersection(LPTfb, 1, 1, 2)
print('Intersection')
print(str(lpt_intersection))
print(str(LPTfb[lpt_intersection,1].astype('int')))

intersection_timestamp_begin =  LPTfb[lpt_intersection[0],0]
intersection_timestamp_end =  LPTfb[lpt_intersection[-1],0]

#fig = plt.figure(figsize=(8.5,11))
ax = fig0.add_subplot(133)
lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax, LPTfb, options, text_color='k')
ax.set_title('Forward and Backwards, No Branches')
ax.set_yticks(dt_list[::8])
ax.grid()
ax.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))

plt.savefig('test_tracking.png')

fn_tc = './TEST.txt'
lpt.lptio.lpt_system_tracks_output_ascii(fn_tc, TIMECLUSTERS0)



#plt.show()

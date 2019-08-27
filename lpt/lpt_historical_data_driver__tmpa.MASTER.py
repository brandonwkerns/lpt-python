import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage
from lpt_historical_data_driver import *

plt.close('all')
plt.ioff()

"""
**** LPT driver script for historical TRMM TMPA "3B42" data. ****


Example usage:
   python lpt_historical_data_driver__tmpa.py 2011060100 2012063021


*** Note: raw data directory on your system is set in readdata.py ***
"""

"""
Dataset Case Settings
"""
dataset={}
dataset['label'] = 'tmpa'
dataset['raw_data_parent_dir'] = '/home/orca/data/satellite/trmm_global_rainfall'
dataset['sub_directory_format'] = '%Y'
dataset['data_time_interval'] = 3           # Time resolution of the data in hours.
dataset['read_function'] = lpt.readdata.read_tmpa_at_datetime
dataset['verbose'] = True

"""
Main settings for lpt
"""
## Plot settings.
plotting = {}
plotting['do_plotting'] = True               # True or False -- Should I make plots?
plotting['plot_area'] = [0, 360, -50, 50]   # Plotting area for maps.
plotting['time_lon_range'] = [40, 200]       # Longitude Range for time-longitude plots.

## High level output directories. Images and data will go in here.
output={}
output['img_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images'
output['data_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'
output['sub_directory_format'] = '%Y/%m/%Y%m%d'

## LP Object settings
lpo_options={}
#lpo_options['do_lpo_calc'] = True
lpo_options['do_lpo_calc'] = False
lpo_options['thresh'] = 12.0                 # LP Objects threshold
lpo_options['accumulation_hours'] = 72       # Accumulation period for LP objects.
lpo_options['filter_stdev'] = 20             # Gaussian filter width, in terms of grid points.

## LPT Settings
lpt_options={}
#lpt_options['do_lpt_calc'] = False
lpt_options['do_lpt_calc'] = True
lpt_options['min_overlap_points'] = 1600      # LP object connectivity is based on points
lpt_options['min_overlap_frac'] = 0.5         # -- OR fraction of either LP object.
lpt_options['min_lp_objects_points'] = 400    # Disregard LP objects smaller than this.
lpt_options['min_lpt_duration_hours'] = 7*24  # Minumum duration to keep it as an LPT
lpt_options['center_jump_max_hours'] = 3*24   # How long to allow center jumps

## Merging/Splitting settings
merge_split_options={}
merge_split_options['allow_merge_split'] = True
#merge_split_options['allow_merge_split'] = False
merge_split_options['split_merger_min_hours'] = 72     # Min duration of a split/merging track to separate it.

"""
Call the real time driver function.
"""

lpt_historical_data_driver(dataset,plotting,output,lpo_options,lpt_options, merge_split_options, sys.argv)

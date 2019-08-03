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
**** Historical Data LPT driver script for MERRA2 Rainfall data. ****

You need to provied the beginning and ending date as command line args.

Example usage:
   python lpt_historical_data_driver__merra2.py 2011110100 2011113018

*** Note: raw data directory on your system is set in readdata.py ***
"""

"""
Dataset Case Settings
"""
dataset={}
dataset['label'] = 'merra2'
dataset['raw_data_parent_dir'] = '/home/orca/asavarin/LPT/MERRA2'
dataset['sub_directory_format'] = '%Y'
dataset['data_time_interval'] = 6           # Time resolution of the data in hours.
dataset['read_function'] = lpt.readdata.get_merra2_6h_rain
dataset['verbose'] = False

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
lpo_options['do_lpo_calc'] = False
lpo_options['thresh'] = 12.0                 # LP Objects threshold
lpo_options['accumulation_hours'] = 72       # Accumulation period for LP objects.
lpo_options['filter_stdev'] = [10,8]         # Gaussian filter width, in terms of grid points.
                                             # If it's a list, specify [x,y] standard deviation.

## LPT Settings
lpt_options={}
lpt_options['do_lpt_calc'] = True
#lpt_options['do_lpt_calc'] = False
lpt_options['lpt_history_days'] = 60          # How many days to go back for LPT tracking and time-lon plot.
lpt_options['min_overlap_points'] = 12000     # LP object connectivity is based on points
lpt_options['min_overlap_frac'] = 0.5         # -- OR fraction of either LP object.
lpt_options['min_lp_objects_points'] = 300    # Disregard LP objects smaller than this.
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

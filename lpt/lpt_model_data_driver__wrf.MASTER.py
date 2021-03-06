import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage
from lpt_model_data_driver import *

plt.close('all')
plt.ioff()

"""
**** LPT Driver Script for WRF model output. ****

No command line args needed. Time range will be detected from wrfout_d01 files.

Example usage:
   python lpt_model_data_driver__wrf.py

*** Note: the path of wrfout data is set as dataset['raw_data_parent_dir'] below. ***
"""

"""
Dataset Case Settings
"""
dataset={}
dataset['label'] = 'wrf'
dataset['raw_data_parent_dir'] = '/home/orca/asavarin/umcm/output/ind_20111122_ecmwf_p'
dataset['sub_directory_format'] = ''
dataset['data_time_interval'] = 1           # Time resolution of the data in hours.
dataset['read_function'] = lpt.readdata.get_wrfout_rain
dataset['verbose'] = True

"""
Main settings for lpt
"""
## Plot settings.
plotting = {}
plotting['do_plotting'] = True               # True or False -- Should I make plots?
plotting['plot_area'] = [30, 180, -30, 30]   # Plotting area for maps.
plotting['time_lon_range'] = [40, 200]       # Longitude Range for time-longitude plots.

## High level output directories. Images and data will go in here.
output={}
output['img_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images'
output['data_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'
output['sub_directory_format'] = 'ind_20111122_ecmwf_p'

## LP Object settings
lpo_options={}
#lpo_options['do_lpo_calc'] = True
lpo_options['do_lpo_calc'] = False
lpo_options['thresh'] = 15.0                 # LP Objects threshold
lpo_options['accumulation_hours'] = 72       # Accumulation period for LP objects.
lpo_options['filter_stdev'] = 14             # Gaussian filter width, in terms of grid points.

## LPT Settings
lpt_options={}
#lpt_options['do_lpt_calc'] = False
lpt_options['do_lpt_calc'] = True
lpt_options['min_overlap_points'] = 830      # LP object connectivity is based on points
lpt_options['min_overlap_frac'] = 0.5         # -- OR fraction of either LP object.
lpt_options['min_lp_objects_points'] = 200    # Disregard LP objects smaller than this.
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

lpt_model_data_driver(dataset,plotting,output,lpo_options,lpt_options, merge_split_options, sys.argv)

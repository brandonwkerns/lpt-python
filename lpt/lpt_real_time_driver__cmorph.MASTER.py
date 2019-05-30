import matplotlib; matplotlib.use('agg')
import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt
import sys
import os
import matplotlib.colors as colors
import scipy.ndimage
from lpt_real_time_driver import *

plt.close('all')
plt.ioff()

"""
**** Real-time LPT driver script for CMORPH RT data. ****

If no command line arguements, it uses the latest data closest to real time.
The date and how far back to calculate can be specified on the command line
which can be useful for filling in previous/missing data.

Example usage:
   python lpt_real_time_driver__cmorph.py
   python lpt_real_time_driver__cmorph.py 2019040100 72


*** Note: raw data directory on your system is set in readdata.py ***
*** Note: the "72" back filling time is for updating LP Objects.  ***
***    LPT goes back "history_length" specified below.            ***
"""

"""
Dataset Case Settings
"""
dataset={}
dataset['label'] = 'cmorph'
dataset['raw_data_parent_dir'] = '/home/orca/data/satellite/cmorph/'
dataset['sub_directory_format'] = '%Y/%m/%Y%m%d'
dataset['data_time_interval'] = 1           # Time resolution of the data in hours.
dataset['read_function'] = lpt.readdata.read_cmorph_at_datetime
dataset['verbose'] = True
dataset['sub_area'] = [40,210,-40,40] # Use for CMORPH only. Very slow if you use full global data.

"""
Main settings for lpt
"""
## Plot settings.
plotting = {}
plotting['do_plotting'] = True               # True or False -- Should I make plots?
plotting['plot_area'] = [50, 200, -30, 30]   # Plotting area for maps.
plotting['time_lon_range'] = [40, 200]       # Longitude Range for time-longitude plots.

## High level output directories. Images and data will go in here.
output={}
output['img_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images'
output['data_dir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data'
output['sub_directory_format'] = '%Y/%m/%Y%m%d'

## LP Object settings
lpo_options={}
lpo_options['thresh'] = 12.0                 # LP Objects threshold
lpo_options['accumulation_hours'] = 72       # Accumulation period for LP objects.
lpo_options['filter_stdev'] = 70             # Gaussian filter width, in terms of grid points.

## LPT Settings
lpt_options={}
lpt_options['do_lpt_calc'] = True
#lpt_options['do_lpt_calc'] = False
lpt_options['lpt_history_days'] = 60          # How many days to go back for LPT tracking and time-lon plot.
lpt_options['min_overlap_points'] = 48000      # LP object connectivity is based on points
lpt_options['min_overlap_frac'] = 0.5         # -- OR fraction of either LP object.
lpt_options['min_lp_objects_points'] = 1200    # Disregard LP objects smaller than this.
lpt_options['min_lpt_duration_hours'] = 7*24  # Minumum duration to keep it as an LPT
lpt_options['center_jump_max_hours'] = 3*24   # How long to allow center jumps


"""
Call the real time driver function.
"""

lpt_real_time_driver(dataset,plotting,output,lpo_options,lpt_options, None, sys.argv)

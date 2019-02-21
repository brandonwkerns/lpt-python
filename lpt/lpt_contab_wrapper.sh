#!/bin/bash


## Get in to script directory.
cd /home/orca/bkerns/lib/lpt/lpt-python/lpt

## Activate the Anaconda Python module with all the dependencies.
source /home/disk/atmos/bkerns/anaconda3/bin/activate meteo
#PATH=/home/users/user_name/miniconda2/envs/my_env/bin
#PATH=/home/disk/atmos/bkerns/anaconda3/envs/meteo/bin

## Call the Python driver scripts.
echo Updating LPT.
python lpt_real_time_driver__tmpa.py >& log.rt.tmpa
python lpt_real_time_driver__cmorph.py >& log.rt.cmorph

##
## Update the animations.
##
echo Updating Animations.
cd /home/orca/bkerns/public_html/realtime_mjo_tracking/images/lpt

# Last 3 days
file_list=`ls -1 lp_objects_tmpa_rt_*.png | tail -25 | tr '\n' ' '`
/usr/bin/convert -delay 20 $file_list lp_objects_tmpa_rt_LAST3DAYS.gif

file_list=`ls -1 lp_objects_cmorph_rt_*.png | tail -73 | tr '\n' ' '`
/usr/bin/convert -delay 10 $file_list lp_objects_cmorph_rt_LAST3DAYS.gif

# Last 10 days
file_list=`ls -1 lp_objects_tmpa_rt_*.png | tail -81 | tr '\n' ' '`
/usr/bin/convert -delay 20 $file_list lp_objects_tmpa_rt_LAST10DAYS.gif

file_list=`ls -1 lp_objects_cmorph_rt_*.png | tail -241 | tr '\n' ' '`
/usr/bin/convert -delay 10 $file_list lp_objects_cmorph_rt_LAST10DAYS.gif




echo Done.

exit 0

#!/bin/bash


## Get in to script directory.
cd /home/orca/bkerns/lib/lpt/lpt-python/lpt

## Activate the Anaconda Python module with all the dependencies.
source /home/disk/atmos/bkerns/anaconda3/bin/activate meteo

## Call the Python driver scripts.
echo Updating LPT.
python lpt_real_time_driver__tmpa.py >& log.rt.tmpa
python lpt_real_time_driver__cmorph.py >& log.rt.cmorph

##
## Update the animations.
##
echo Updating Animations.
cd /home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images
rm *.png

# Last 3 days
#file_list=`find tmpa/objects/*/*/* | grep .png | tail -25 | tr '\n' ' '`
ln -s `find tmpa/objects/*/*/* | grep .png | tail -25 ` .
/usr/bin/convert -delay 20 *.png lp_objects_tmpa_rt_LAST3DAYS.gif
rm *.png

ln -s `find cmorph/objects/*/*/* | grep .png | tail -73 ` .
/usr/bin/convert -delay 20 *.png lp_objects_cmorph_rt_LAST3DAYS.gif
rm *.png


#file_list=`ls -1 lp_objects_cmorph_rt_*.png | tail -73 | tr '\n' ' '`
#/usr/bin/convert -delay 10 $file_list lp_objects_cmorph_rt_LAST3DAYS.gif

# Last 10 days

#file_list=`ls -1 lp_objects_tmpa_rt_*.png | tail -81 | tr '\n' ' '`
#/usr/bin/convert -delay 20 $file_list lp_objects_tmpa_rt_LAST10DAYS.gif

#file_list=`ls -1 lp_objects_cmorph_rt_*.png | tail -241 | tr '\n' ' '`
#/usr/bin/convert -delay 10 $file_list lp_objects_cmorph_rt_LAST10DAYS.gif


ln -s `find tmpa/objects/*/*/* | grep .png | tail -81 ` .
/usr/bin/convert -delay 20 *.png lp_objects_tmpa_rt_LAST10DAYS.gif
rm *.png

ln -s `find cmorph/objects/*/*/* | grep .png | tail -241 ` .
/usr/bin/convert -delay 20 *.png lp_objects_cmorph_rt_LAST10DAYS.gif
rm *.png

## Make pause at the beginning and end of animation repeats.
for ff in *.gif
do
  convert $ff \( +clone -set delay 100 \) +swap +delete $ff
done




echo Done.

exit 0

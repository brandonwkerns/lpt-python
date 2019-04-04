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

# Last 10 days

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

## Get latest time-longitude plots.
ln -sf `find tmpa/systems/*/*/* | grep .png | tail -1 ` lpt_time_lon_tmpa_LATEST.png
ln -sf `find cmorph/systems/*/*/* | grep .png | tail -1 ` lpt_time_lon_cmorph_LATEST.png 



echo Done.

exit 0

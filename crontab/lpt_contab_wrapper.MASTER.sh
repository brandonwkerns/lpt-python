#!/bin/bash


## Get in to script directory.
cd /home/orca/bkerns/lib/lpt/lpt-python/lpt

## Activate the Anaconda Python module with all the dependencies.
source /home/disk/atmos/bkerns/anaconda3/bin/activate meteo

## Call the Python driver scripts.
echo Updating LPT...
echo TMPA
python lpt_real_time_driver__tmpa.py >& log.rt.tmpa
echo CMORPH
python lpt_real_time_driver__cmorph.py >& log.rt.cmorph
echo CFS Forecast
python lpt_real_time_driver__cfs_forecast.py >& log.rt.cfs

##
## Update the animations.
##
echo Updating Animations.
cd /home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/images
rm -f *.png *.gif

ln -s `find cfs/objects/*/*/* | grep .png | tail -168 ` .
/usr/bin/convert -delay 15 *.png lp_objects_cfs_FCST45DAYS.gif
rm *.png

ln -s `find cmorph/objects/*/*/* \( -name "*00.png" -or -name "*06.png" -or -name "*12.png" -or -name "*18.png" \)  | tail -241` .
/usr/bin/convert -delay 15 *.png lp_objects_cmorph_rt_LAST45DAYS.gif
rm *.png

ln -s `find tmpa/objects/*/*/* \( -name "*.png" \)  | tail -481` .
/usr/bin/convert -delay 15 *.png lp_objects_tmpa_rt_LAST45DAYS.gif
rm *.png

## Make pause at the beginning and end of animation repeats.
for ff in *.gif
do
  convert $ff \( +clone -set delay 100 \) +swap +delete $ff
done

## Get latest time-longitude plots.
ln -sf `find tmpa/systems/*/*/* | grep .png | tail -1 ` lpt_time_lon_tmpa_LATEST.png
ln -sf `find cmorph/systems/*/*/* | grep .png | tail -1 ` lpt_time_lon_cmorph_LATEST.png
ln -sf `find cfs/systems/*/*/* | grep .png | tail -1 ` lpt_time_lon_cfs_LATEST.png



echo Done.

exit 0

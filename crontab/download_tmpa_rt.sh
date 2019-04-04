#!/bin/sh -fx

#######  Download 3B42 real-time data from NASA.

ftpsite=ftp://trmmopen.gsfc.nasa.gov/pub/merged/mergeIRMicro

hours=" 00 03 06 09 12 15 18 21 "

workdir=/home/orca/data/satellite/trmm_global_rainfall

VER=7

#######################################################################


cd $workdir

today=`date -u +%Y%m%d`
yyyy=`date -u +%Y`
mm=`date -u +%m`

      
for hh in $hours
  do
	
  filewanted=3B42RT.${today}${hh}.$VER.bin.gz
  /usr/bin/wget $ftpsite/$yyyy/$filewanted      
  
  if [ -e $filewanted ]
      then
      
      mkdir -p rt/$yyyy/$mm/$today
      mv $filewanted rt/$yyyy/$mm/$today
      /bin/gunzip -f  rt/$yyyy/$mm/$today/$filewanted
      
  fi
    
done


yesterday=`date --date=${today}-1day  +%Y%m%d`
yyyy=`date --date=${today}-1day +%Y`
mm=`date --date=${today}-1day +%m`

      
for hh in $hours
  do
	
  filewanted=3B42RT.${yesterday}${hh}.$VER.bin.gz
  /usr/bin/wget $ftpsite/$yyyy/$filewanted      
  
  if [ -e $filewanted ]
      then
      
      mkdir -p rt/$yyyy/$mm/$yesterday
      mv $filewanted rt/$yyyy/$mm/$yesterday
      /bin/gunzip -f  rt/$yyyy/$mm/$yesterday/$filewanted
      
  fi
    
done

echo Done.

exit 0


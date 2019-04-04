#!/bin/sh -fx

#######  Download 3B42 real-time data from NASA.

ftpsite=ftp://ftp.cpc.ncep.noaa.gov/precip/CMORPH_RT/GLOBE/data
workdir=/home/orca/data/satellite/cmorph

#######################################################################


cd $workdir

today=`date -u +%Y%m%d`
yyyy=`date -u +%Y`
mm=`date -u +%m`

      
for hh in {00..23}
  do
	
  filewanted=CMORPH_V0.x_RT_8km-30min_$today$hh.gz
  /usr/bin/wget $ftpsite/$yyyy/$yyyy$mm/$filewanted      
  
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

      
for hh in {00..23}
  do
	
  filewanted=CMORPH_V0.x_RT_8km-30min_$yesterday$hh.gz
  /usr/bin/wget $ftpsite/$yyyy/$yyyy$mm/$filewanted      

  
  if [ -e $filewanted ]
      then
      
      mkdir -p rt/$yyyy/$mm/$yesterday
      mv $filewanted rt/$yyyy/$mm/$yesterday
      /bin/gunzip -f  rt/$yyyy/$mm/$yesterday/$filewanted
      
  fi
    
done

echo Done.

exit 0


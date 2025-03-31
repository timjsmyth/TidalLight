#!/usr/bin/bash

# Datum fraction (e.g. 0.5 = 50%)
DP=0.25

# Time increment in hours (e.g. 0.25 is 15 minutes)
TI=0.5

# Additional depth (in metres)
#DEPTH=-4.0

# Start date
#START="2020-01-01"
# End date
#END="2020-01-05"

# Start date
START="1985-01-01"
# End date
END="1985-01-06"

TPXO_dir="/data/proteus1/scratch/gle/TPXO9_atlas_v5_nc/"

FILES=("Requests/Ruben/"*Ruben_*.csv)
for f in "${FILES[@]}"
do
  echo "Creating symbolic link to: ${f}"
  ln -sf "${f}" .
  xbase=${f##*/}
  
  # Run the code here
  xpref=${xbase%.*}
  location=`echo $xpref | sed -e 's/Kd_Falchi_Output_//g'`
  echo $location
  location_number=`echo $location | sed -e 's/Ruben_//g'`
  echo $location_number
  DEPTH=`cat Ruben_station_depths.txt | awk -F ',' -v station="$location_number $location" '$1 ~ station {print $2}'`
  echo $DEPTH
  echo " Running code for $location"
  /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -d $DEPTH -o -loc $location -p -start $START -end $END -TPXO_d $TPXO_dir
  echo "  Removing symbolic link: $xbase"
  /bin/rm $xbase
  echo "============================"
  exit 0
done

exit 0


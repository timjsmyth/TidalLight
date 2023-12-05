#!/usr/bin/bash

count=1
while read latlon; do
  #echo "$latlon"
  lat=`echo $latlon | awk -F "," '{print $1}'`
  lon=`echo $latlon | awk -F "," '{print $2}'`
  #echo $lat $lon
  echo "./Falchi_Kd_Position.py -lat $lat -lon $lon -ofile Requests/Kd_Falchi_Output_Ruben_${count}.csv"
  ./Falchi_Kd_Position.py -lat $lat -lon $lon -ofile Requests/Kd_Falchi_Output_Ruben_${count}.csv
  count=$((count+1))
  echo $count
done < Requests/Sites_dedup.csv

exit 0

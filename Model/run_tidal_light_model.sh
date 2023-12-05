#!/usr/bin/bash
## declare an array variable holding the locations

#declare -a locations=("Plymouth_Dockyard" "Tokyo" "Eilat" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")
declare -a locations=("Tokyo")
declare -a locations=("Ruben_1")

# Datum fraction (e.g. 0.5 = 50%)
DP=0.25

# Time increment in hours (e.g. 0.25 is 15 minutes)
TI=0.25

# Start date
START="2020-01-01"
# End date
END="2020-01-31"

## MULTIPLE YEAR RUN ##
## loop through the locations array
#for location in "${locations[@]}"
#do
#   echo "$location"
#   year="2001"
#   while [ $year -lt 2020 ] 
#   do 
#      echo $year
#      START=$year"-01-01"
#      END=$[$year+1]"-01-01"
#      echo $START
#      echo $END
      ## This version without plots (long runs)
#      /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -start $START -end $END #-TC
#      year=$[$year+1]
#   done
#done

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   ## Run the TidalLight_Model.py
   ## This version if you want to plot (i.e. short duration runs)
   #/bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -p -start $START -end $END #-TC
   /bin/python3 TidalLight_Model.py -s -A 2 -l -T $TI -dp $DP -o -loc $location -p -start $START -end $END #-TC

   ## This version without plots (long runs)
   #/bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -start $START -end $END #-TC

done

exit 0


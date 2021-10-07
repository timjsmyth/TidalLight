#!/usr/bin/bash
## declare an array variable holding the locations

declare -a locations=("Tokyo" "Eilat" "Plymouth_Dockyard" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")
declare -a locations=("BuenosAires")

# Datum fraction (e.g. 0.5 = 50%)
DP=0.25

# Time increment in hours (e.g. 0.25 is 15 minutes)
TI=0.25

# Start date
START="2021-09-28"
# End date
END="2021-10-10"

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   # Run the TidalLight_Model.py
   /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -p -start $START -end $END #-TC
done

exit 0


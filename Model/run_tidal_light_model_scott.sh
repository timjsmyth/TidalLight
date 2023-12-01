#!/usr/bin/bash
## declare an array variable holding the locations
declare -a locations=("Plymouth_Dockyard")

# Datum fraction (e.g. 0.5 = 50%)
DP=0.25

# Time increment in hours (e.g. 0.25 is 15 minutes)
TI=12.0

# Start date
START="1990-01-01"
# End date
END="2029-12-31"

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   ## This version without plots (long runs)
   /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -start $START -end $END
done

exit 0


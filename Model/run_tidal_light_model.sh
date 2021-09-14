#!/usr/bin/bash
## declare an array variable holding the locations
declare -a locations=("Tokyo" "Eilat" "Plymouth_L4" "Plymouth_Dockyard" "Lagos" "NewYork" "LosAngeles" "BuenosAires")
declare -a locations=("BuenosAires")

# Datum fraction (0.5 = 50%)
DP=0.25

# Time increment (0.25 is 15 minutes)
TI=0.25

# Start date
START="2021-07-20"
# End date
END="2021-07-25"

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   # Run the TidalLight_Model.py
   /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -p -start $START -end $END

done


exit 0


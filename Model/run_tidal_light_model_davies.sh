#!/usr/bin/bash
## declare an array variable holding the locations
#Site_85: 32.167, 130.03 # Japan
#Site_78: -5.00, 119.00  # Indonesia
#Site_183: 29.50, 34.92  # Eilat
#Site_80: -5.06, 119.329 # Indonesia

#declare -a locations=("Plymouth_Dockyard" "Tokyo" "Eilat" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")
declare -a locations=("Site_85" "Site_78" "Site_183" "Site_80")

# Datum fraction (e.g. 0.5 = 50%)
DP=0.25

# Time increment in hours (e.g. 0.25 is 15 minutes)
TI=0.08333

# Start date
START="2020-03-01"
# End date
END="2020-05-31"

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   ## This version without plots (long runs)
   /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -start $START -end $END
done


# Start date
START="2020-09-01"
# End date
END="2020-11-30"

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   ## Run the TidalLight_Model.py
   ## This version if you want to plot (i.e. short duration runs)
   /bin/python3 TidalLight_Model.py -s -A 2 -l -t -T $TI -dp $DP -o -loc $location -start $START -end $END

done

exit 0


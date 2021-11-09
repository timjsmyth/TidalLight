#!/usr/bin/bash
## declare an array variable holding the locations

declare -a locations=("Plymouth_Dockyard" "Tokyo" "Eilat" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")

for location in "${locations[@]}"
do
   echo "$location"
   analyse_dosage.py --loc $location
done

exit 0

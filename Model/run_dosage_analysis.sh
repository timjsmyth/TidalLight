#!/usr/bin/bash
## declare an array variable holding the locations

#declare -a locations=("Plymouth_Dockyard" "Tokyo" "Eilat" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")
declare -a locations=("Plymouth_Dockyard" "Tokyo" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")

#declare -a locations=("Shanghai")

for location in "${locations[@]}"
do
   echo "$location"
   analyse_dosage.py --loc $location -p
done

exit 0

#for location in "${locations[@]}"
#do
#   echo "$location"
#   analyse_dosage.py --loc $location -p -m
#done

#analyse_dosage.py -p -pos intertidal
analyse_dosage.py -p -pos surface
analyse_dosage.py -p -pos intertidal

exit 0

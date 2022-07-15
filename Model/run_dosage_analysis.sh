#!/usr/bin/bash
## declare an array variable holding the locations

declare -a locations=("Plymouth_Dockyard" "Tokyo" "Eilat" "Lagos" "NewYork" "LosAngeles" "BuenosAires" "Shanghai" "Mumbai")

#declare -a locations=("Tokyo")

#for location in "${locations[@]}"
#do
#   echo "$location"
#   analyse_dosage.py --loc $location -p
#done

#for location in "${locations[@]}"
#do
#   echo "$location"
#   analyse_dosage.py --loc $location -p -m
#done

#analyse_dosage.py -p -pos intertidal
analyse_dosage.py -p -pos surface

exit 0

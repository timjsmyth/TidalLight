#!/bin/bash
# FUNCTION:
# script to reduce all tidal data files to 'datetime, tide (m)' format. 
# Adam wright 16/09/20
readfile=TideRead_LOCATION
filepath=$(nawk '/</ {} NR==4' $readfile) # retrieve variable from row 4 of the text file 
filepath=${filepath#?} # remove "<" marker
echo 	"FILEPATH:"	$filepath # print variable
filename=$(nawk '/</ {} NR==2' $readfile)
filename=${filename#?}
echo	"FILENAME:"	$filename
fname_new="${filename%????}_new.csv" # rename variable by removing last 4 characters
echo	"NEW FILENAME:"	$fname_new
pattern=$(nawk '/</ {} NR==6' $readfile)
pattern=${patttern#?}
echo 	"PATTERN:"		$pattern
tidal_col=$(nawk '/</ {} NR==8' $readfile); tidal_col=${tidal_col#?} 
date_col=$(nawk '/</ {} NR==10' $readfile); date_col=${date_col#?}
start_line=$(grep -n -m 1 $pattern $filename) # Find first instance of the year/pattern match
echo "HEADER LINE" $start_line
start_row=${start_line:0:2} # Retrieve the row number of the first data recording for the year/pattern CAREFUL OF SINGLmE (0:1) AND DOUBLE DIGITS (0:2) 
echo "Original file:"
head -5 $filename
sed -e 1,${start_row}d $filename > new-file # rewmove all lines prior to data and save to a tmp file (new-file)
sed -i "s/.*)//" new-file #remove everything prior to the first )
sed -i "s/   \+/, /g" new-file # Convert any sequence of 3 or more blank spaces to a ,
sed -i "s/[^0-9 .,:/ -]*//g" new-file # remove all characteers except those inlcuded in the []
cut -d ',' -f $date_col,$tidal_col new-file > $fname_new
rm new-file

echo "New file:"
head $fname_new

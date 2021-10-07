#!/bin/bash
# FUNCTION:
# script to reduce all tidal data files to 'datetime, tide (m)' format. 
# Adam Wright 16/09/2020
# Tim Smyth 13/09/2021
readfile=TideReadShanghai
filepath=$(nawk '/</ {} NR==4' $readfile) # retrieve variable from row 4 of the text file 
filepath=${filepath#?} # remove "<" marker
echo 	"FILEPATH:"	$filepath # print variable
filename=$(nawk '/</ {} NR==2' $readfile)
filename=${filename#?}
echo	"FILENAME:"	$filename
fstem="${filename%???.*}" # remove extension from filename and _L0
fname_new="${fstem}_L1.csv" # rename to L1
echo	"NEW FILENAME:"	$fname_new
pattern=$(nawk '/</ {} NR==6' $readfile)
pattern=${pattern#?}
echo 	"PATTERN:"		$pattern
tidal_col=$(nawk '/</ {} NR==8' $readfile); tidal_col=${tidal_col#?} 
date_col=$(nawk '/</ {} NR==10' $readfile); date_col=${date_col#?}
#start_line=$(grep -n $pattern $filename | head -1) # Find first instance of the year/pattern match
start_line=$(grep -n -m 1 $pattern $filename)
echo "HEADER LINE" $start_line
start_row=$(grep -n -m 1 $pattern $filename |sed 's/\([0-9]*\).*/\1/') # Retrieve the row number of the first data recording for the year/pattern CAREFUL OF SINGLE (0:1) AND DOUBLE DIGITS (0:2) 
echo $start_row
echo "Original file:"
head -5 $filename
sed -e 1,${start_row}d $filename > new-file # remove all lines prior to data and save to a tmp file (new-file)
sed '/-32767/d' new-file > tmp.txt
mv tmp.txt new-file

cat new-file | awk 'BEGIN {FS=","};{printf("%04d/%02d/%02d %02d:00:00, %5.3f\n", $1,$2,$3,$4,$5/1000.)}' > $fname_new
rm new-file

echo "New file:"
head $fname_new

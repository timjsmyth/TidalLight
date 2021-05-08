#!/bin/bash
#ead2) script to reduce all tidal data files to 'datetime, tide (m)' format. 
readfile=TideReadEilatGulf
filepath=$(nawk '/</ {} NR==4' $readfile) # retrieve variable from row 4 of the text file 
filepath=${filepath#?} # remove "<" marker
echo 	"FILEPATH:"	$filepath # print variable
filename=$(nawk '/</ {} NR==2' $readfile)
filename=${filename#?}
echo	"FILENAME:"	$filename
fname_new="${filename%????}_new.csv" # rename variable by removing last 4 characters
echo	"NEW FILENAME:"	$fname_new
year=$(nawk '/</ {} NR==6' $readfile)
year=${year#?}
echo 	"YEAR:"		$year
tidal_col=$(nawk '/</ {} NR==8' $readfile); tidal_col=${tidal_col#?} 
date_col=$(nawk '/</ {} NR==10' $readfile); date_col=${date_col#?}
start_line=$(grep -n -m 1 $year $filename) # Find first instance of the year
echo $start_line
start_row=${start_line:0:1} # Retrieve the row number of the first data recording for the year of interest 
#sed -i "s/${year}/${year},/g" $filename
echo "Original file:"
head $filename
sed -e 1,${start_row}d $filename > $fname_new   # remove all daa + headers prior to first instance and save to a new filem
echo "first new-file"
head $fname_new
#!/bin/bash
cut -d ',' -f 1-2 $fname_new | sed -i 's/,/\ /' $fname_new # replace the first comma with a whitespace 
echo "New file:"
head $fname_new

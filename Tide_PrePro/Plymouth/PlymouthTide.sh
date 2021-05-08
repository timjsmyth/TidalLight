#!/bin/bash
#ead2) script to reduce all tidal data files to 'datetime, tide (m)' format. 
readfile=TideReadPlymouth
filepath=$(nawk '/</ {} NR==4' $readfile) # retrieve variable from row 4 of the text file 
filepath=${filepath#?} # remove "<" marker
echo 	"FILEPATH:"	$filepath # print variable
filename=$(nawk '/</ {} NR==2' $readfile)
filename=${filename#?}
echo	"FILENAME:"	$filename
fname_new="${filename%????}_new.csv" # rename variable by removing last 4 characters
echo	"NEW FILENAME:"	$fname_new
pattern=$(nawk '/</ {} NR==6' $readfile)
pattern=${pattern#?}
echo 	"PATTERN:"		$pattern
tidal_col=$(nawk '/</ {} NR==8' $readfile); tidal_col=${tidal_col#?} 
date_col=$(nawk '/</ {} NR==10' $readfile); date_col=${date_col#?}
start_line=$(grep -n -m 1 $pattern $filename) # Find first instance of the year
echo $start_line
start_row=${start_line:0:2} # Retrieve the row number of the first data recording for the year of interest 
#sed -i "s/${year}/${year},/g" $filename
echo "Original file:"
head -5 $filename
sed -e 1,${start_row}d $filename > new-file
sed -i "s/.*) //" new-file 
sed -i "s/   \+/, /g" new-file 
sed -i "s/[^0-9 .,:/ -]*//g" new-file 
cut -d ',' -f $date_col,$tidal_col new-file > $fname_new
rm new-file
# rewmove all daa + headers prior to first instance and save to a new file
#!/bin/bash
#sed -i 's/.*)//' $fname_new #remove everything prior to the first )
#sed -i 's/[^0-9 .,:/ -]*/ /g' $fname_new
#sed -i 's/   \+/, /g' $fname_new # Convert any sequence of 3 or more blank spaces to a , 
#sed -i 's/[^0-9 .,:/ -]*//g' $fname_new # remove all characteers except those inlcuded in the []
#echo "AMENDED FILE"
#head $fname_new
#cut -d ',' -f 1,2 $fname_new > $fname_new 
echo "New file:"
head $fname_new

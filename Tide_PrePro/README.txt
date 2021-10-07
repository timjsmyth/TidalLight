Make a copy of both "LOCATION_Tide.sh" and "TideRead_LOCATION"

Replace _LOCATION with a relevant Geotag
Create a new folder with the same name as the Geotag to place the amended copies into. 
Place the data file of the Geoteg locations tide gauge into the same folder. 

open Linux terminal and edit TideRead<Geotag> with and editor of your choice.
e.g. vim TideReadPlymouth

Edit the file such that all lines beginnig with '<' have the appropriate information. 

Edit the <Geotag>Tide.sh file such that the output is two columns: 
column 1 = date, 
column 2 = sealevel 

NOTE if sealevel is not in 'm' this must be changed withing the python script

Common things to change: 
- start_row=${start_line:0:n} - n must be an integer, this will be 1 if the header line </= 9, 2 if </99, 3 if </= 999 
- sed statements this will depend on the source of the dataset

Overall goal is to remove all but the required data, LOCATION_Tide.sh serves as a basic guide to assist with this. 

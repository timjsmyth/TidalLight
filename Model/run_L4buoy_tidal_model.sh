#!/usr/bin/bash
## declare an array variable holding the locations
declare -a locations=("Plymouth_L4")

# Issues with the Profiler only operating on the Devil's time
# The following assumes that the UK and France shift to Daylight Savings simultaneously
dst_start_end_file=Required/daylight_saving_start_end.txt
YYYY=`date --utc +%Y`
JJJ=`date --utc +%j | sed 's/^0//'`
DST_start=`cat $dst_start_end_file | grep $YYYY | awk '{print $2}' | sed 's/^0//'`
DST_end=`cat $dst_start_end_file | grep $YYYY | awk '{print $3}' | sed 's/^0//'`

if [ "$JJJ" -ge "$DST_start" -a "$JJJ" -le "$DST_end" ]; then 
   ZONE="CEST"
else 
   ZONE="CET"
fi

TIME=$(date +%s --utc -d "12:00:00 $ZONE")
UTC_TIME=$(date +%s --utc -d "12:00:00")
((DIFF=UTC_TIME-TIME))
FROFF=`echo - | awk -v SECS=$DIFF '{printf "%d",SECS/(60*60)}'`

# Make sure that the current time has not been monkeyed with
echo "Current time"
date

# Run model for a week
# Start date
START=`date --date='-1 day' +%Y-%m-%d`
# End date
END=`date --date='7 day' +%Y-%m-%d`

## loop through the locations array
for location in "${locations[@]}"
do
   echo "$location"
   ## Run the L4buoy_tides.py model
   /bin/python3 L4buoy_tides.py -t -o -start $START -end $END #-TC

done

cat Output/L4buoy_tides.txt | sed -e 's/\"//g' > /tmp/L4buoy_tides.txt

# TODAY
echo "======================="
echo "Today's times of Hi/Lo water @Devonport"
TODAY=`date +%Y-%m-%d`
grep $TODAY /tmp/L4buoy_tides.txt
grep $TODAY /tmp/L4buoy_tides.txt > /tmp/L4buoy_tides_today.txt
echo "======================="

echo "Offset in French time from God's own time(UTC): $FROFF"
while read entry
do
   DATETIME=`echo "$entry" | awk '{print $1,$2}'`
   FRENCHTIME=`date -d"$DATETIME $FROFF hour" +"%Y-%m-%d %H:%M:%S"`
   echo $FRENCHTIME
done < /tmp/L4buoy_tides_today.txt

exit 0


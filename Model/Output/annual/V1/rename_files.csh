#!/bin/csh -f

foreach file(ALAN*.csv)
   set nfile = `echo $file | sed -e 's/DATA\_//g' | sed -e 's/1\-366/2020\-01\-01\_2020\-12\-31/g'`
   /bin/mv $file $nfile   
   #echo "/bin/mv $file $nfile"   
end

foreach file(Lunar*.csv)
   set nfile = `echo $file | sed -e 's/DATA\_//g' | sed -e 's/Day//g' | sed -E "s/\(//g" | sed -E "s/\)//g" | sed -e 's/1\-366/2020\-01\-01\_2020\-12\-31/g'`
   /bin/mv "${file}" $nfile   
   #echo "/bin/mv "${file}" $nfile"   
end

foreach file(Solar*.csv)
   set nfile = `echo $file | sed -e 's/DATA\_//g' | sed -e 's/Day//g' | sed -E "s/\(//g" | sed -E "s/\)//g" | sed -e 's/1\-366/2020\-01\-01\_2020\-12\-31/g'`
   /bin/mv "${file}" $nfile   
   #echo "/bin/mv "${file}" $nfile"   
end

exit(0)


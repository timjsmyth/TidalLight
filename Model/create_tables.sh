# Surface
cat Dosage_*2020.csv | awk 'BEGIN{FS=","} /Tot/{printf("%s,Solar,%.2e,%.2e,%.2e\n%s,Lunar,%.2e,%.2e,%.2e\n%s,ALAN,%.2e,%.2e,%.2e\n%s,Twilight,%.2e,%.2e,%.2e\n", $3,$8,$9,$10,$3,$11,$12,$13,$3,$14,$15,$16,$3,$17,$18,$19)}'

# Intertidal
cat Dosage_*2020.csv | awk 'BEGIN{FS=","} /Tot/{printf("%s,Solar,%.2e,%.2e,%.2e\n%s,Lunar,%.2e,%.2e,%.2e\n%s,ALAN,%.2e,%.2e,%.2e\n%s,Twilight,%.2e,%.2e,%.2e\n", $3,$24,$25,$26,$3,$27,$28,$29,$3,$30,$31,$32,$3,$33,$34,$35)}'



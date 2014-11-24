#!/bin/sh

cd files

for i in ccsrnies cgcm1 cgcm2 csiromk2b echam4 gfdlr15 gfdlr30 hadcm2 hadcm3 ncarpcm; do
	../load-csv.pl /home/websites/www.cics.uvic.ca/htdocs/scenarios/data/gcminfo/125.csv | grep $i | grep "nam" | ../scenarios2netcdf ../$i\_nam.dat
	../load-csv.pl /home/websites/www.cics.uvic.ca/htdocs/scenarios/data/gcminfo/125.csv | grep $i | grep "cdn" | ../scenarios2netcdf ../$i\_cdn.dat
done

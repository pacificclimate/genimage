#!/usr/bin/perl -w

use strict;

use NetCDF;

my($fullfile) = "./sresa1b/rsds/cnrm_cm3/run1/rsds_A1.nc";
chomp($fullfile);
my($scen, $var, $model, $run, $filename) = (split(/\//, $fullfile))[1..4];
my($startyear, $endyear);

my($calendar, @time);
my($ncfile) = NetCDF::open($fullfile, NetCDF::NOWRITE);

if(NetCDF::inq_dimid($ncfile, "time") == 0) {
    print "Time is present!\n";
}

NetCDF::close($ncfile);

#!/usr/bin/perl -w

use strict;

use NetCDF;

my($listfile) = "listing.txt";
my(%months) = ( "Jan" => 0, "Feb" => 1, "Mar" => 2, "Apr" => 3, "May" => 4, "Jun" => 5, "Jul" => 6, "Aug" => 7, "Sep" => 8, "Oct" => 9, "Nov" => 10, "Dec" => 11 );

open($fd, $listfile) or die("Couldn't open " . $merfile . "\n");
while(<$fd>) {
    my($fullfile) = $_;
    chomp($fullfile);
    my($scen, $var, $model, $run, $filename) = split(/\//, $fullfile)[1..4];
    my($startyear, $endyear);

    my($calendar, @time);
    my($ncfile) = NetCDF::open($fullfile);



	print STDERR join(" , ", $file, $time, $esttime, $startyear, $endyear, $oldstartyear, $oldendyear) . "\n";
    }

    print join(",", $file, $bits[3], $bits[1], $bits[4], join("-", $startyear, $startmonth), join("-", $endyear, $endmonth)) . "\n";
}

close($fd);

    if(0) {
	my(@tmp);
	if($file =~ /([0-9]{4})_([A-Za-z]{3})_to_([0-9]{4})_([A-Za-z]{3})/) {
	    # Year_StringMonth_to_Year_StringMonth
	    @tmp = ($1, $months{$2}, $3, $months{$4});
	} elsif($file =~ /([0-9]{4})-([0-9]{2})(.*?)([0-9]{4})-([0-9]{2})/) {
	    # Year_Month*Year_Month
	    @tmp = ($1, $2 - 1, $4, $5 - 1);
	} elsif($file =~ /([0-9]{4})(.*?)([0-9]{4})/) {
	    # Year*Year
	    @tmp = ($1, $startmonth, $3, $endmonth);
	}
	my($oldstartyear, $oldendyear) = ($startyear, $endyear);
	if($#tmp != -1 && $tmp[0] >= 1800 && $tmp[0] <= 3000 && $tmp[2] >= 1800 && $tmp[2] <= 3000) {
	    ($startyear, $startmonth, $endyear, $endmonth) = @tmp;
	}
    my($esttime) = 12 * ($endyear - $startyear) + ($endmonth - $startmonth) + 1;
    if(defined($time) && $time != $esttime) {
	if(0) {
	if($oldstartyear > $startyear) { 
	    print STDERR join(",", $file, $startyear, $tmp[0]) . "\n";
	}
	if($oldendyear < $endyear) {
	    print STDERR join(",", $file, $endyear, $tmp[2]) . "\n";
	}
}	

    }


#!/usr/bin/perl -w

use strict;

use NetCDF;

my($merfile) = "modelexptrun.txt";
my($listfile) = "listing.txt";

my(%mer);

my($fd);

my(%months) = ( "Jan" => 0, "Feb" => 1, "Mar" => 2, "Apr" => 3, "May" => 4, "Jun" => 5, "Jul" => 6, "Aug" => 7, "Sep" => 8, "Oct" => 9, "Nov" => 10, "Dec" => 11 );

open($fd, $merfile) or die("Couldn't open " . $merfile . "\n");
while(<$fd>) {
    chomp;
    my(@bits) = split(/,/);
    $mer{join("", @bits[0..2])} = [ @bits[3..4] ];
}
close($fd);

open($fd, $listfile) or die("Couldn't open " . $merfile . "\n");
while(<$fd>) {
    chomp;
    my($file) = $_;
    my(@bits) = split(/\//, $file);
    my($key) = join("", $bits[3], $bits[1], $bits[4]);
    my($startyear, $endyear) = @{$mer{$key}};
    my($startmonth, $endmonth) = (0, 11);

    if($startyear =~ /([0-9]+)-([0-9]+)/) {
	$startyear = $1;
	$startmonth = $2;
    }
    if($endyear =~ /([0-9]+)-([0-9]+)/) {
	$endyear = $1;
	$endmonth = $2;
    }

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

    my($calendar, $time);
    NetCDF::open();
    open($cmd, "ncdump -h ".$file."|"); 
    while(<$cmd>) { 
	if(/time = UNLIMITED ; \/\/ \(([0-9]+) currently\)/) {
	    $time = $1;
	} elsif(/
    }
    close($cmd);

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

	print STDERR join(" , ", $file, $time, $esttime, $startyear, $endyear, $oldstartyear, $oldendyear) . "\n";
    }

    print join(",", $file, $bits[3], $bits[1], $bits[4], join("-", $startyear, $startmonth), join("-", $endyear, $endmonth)) . "\n";
}

close($fd);

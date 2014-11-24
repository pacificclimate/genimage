#!/usr/bin/perl -w

use strict;

my($merfile) = "modelexptrun.txt";
my($listfile) = "modelexptrunvar.txt";

my(%mer);

my($fd);

open($fd, $merfile) or die("Couldn't open " . $merfile . "\n");
while(<$fd>) {
    chomp;
    my(@bits) = split(/,/);
    $mer{join("", @bits[0..2])} = [ @bits[3..4] ];
}
close($fd);

open($fd, $listfile) or die("Couldn't open " . $merfile . "\n");
print "Model,Expt/Scenario,Run,Variable,Start_year,End_year,";
while(<$fd>) {
    chomp;
    my($line) = $_;
    my(@bits) = split(/,/, $line);
    my($key) = join("", $bits[2], $bits[0], $bits[3]);

    print join(",", $bits[2], $bits[0], $bits[3], $bits[1], @{$mer{$key}}) . ",\n";
}
close($fd);

#!/usr/bin/perl -w

use strict;

while(<>) {
    chomp;
    my(@bits) = split(/,/);

    print join(",", $bits[0], $bits[4], $bits[5]) . "\n";
    system('ncdump -v time ' . $bits[0] . ' | egrep "(currently|time:units|time:calendar|time = )"');
}

#!/usr/bin/perl -w

use strict;

sub is_available {
  my($token) = @_;

  return (!defined($token) || !($token eq "o" || $token eq "0"));
}

my(%mhash);
my($oldmodel) = "";;
my($fd);
my($i);
my($fn) = $ARGV[0];
open($fd, $fn) or die("Couldn't open file!\n");

my(@dictionary) = split(/,/, <$fd>);
my(@timeslices) = ("1961_1990", @dictionary[7..9]);
my(@vars) = @dictionary[10..30];
my(@regions) = ( "cdn", "nam" );
while(<$fd>) {
  my(@line) = split(/,/);
  my($scenario) = $line[5];
  my($model) = lc($line[0]);
  my(@mytimeslices) = ();
  my(@myvars) = ();
  push(@mytimeslices, $timeslices[0]);
  for($i = 0; $i <= 2; $i++) {
    if(is_available($line[$i + 6])) {
      push(@mytimeslices, $timeslices[$i + 1]);
    }
  }
  for($i = 0; $i <= 20; $i++) {
    if(is_available($line[$i + 9])) {
      push(@myvars, $vars[$i]);
    }
  }

  foreach(@regions) {
    my($region) = $_;
    foreach(@mytimeslices) {
      my($ts) = $_;
      foreach(@myvars) {
	my($var) = $_;
	print join("_", $model, $scenario, $ts, $var, $region) . ".dat\n";
      }
    }
  }
  if(!defined($mhash{$model})) {
    foreach(@regions) {
      my($region) = $_;
      print join("_", $model, "slmask", $region) . ".dat\n";
      print join("_", $model, "lats", $region) . ".dat\n";
      print join("_", $model, "longs", $region) . ".dat\n";
    }
  }
  $mhash{$model} = 1;
}

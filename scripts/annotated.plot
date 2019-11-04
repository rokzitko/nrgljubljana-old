#!/usr/bin/env perl

use warnings;

my $cnt = 0;

while(<>) {
  !/^#/ or next;
  chomp;
  if ($_ eq "" ) {
     $cnt++;
     next;
  }
  ($e) = split;
  print "$cnt $e\n";
}

#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin);

exec $^X, "$RealBin/compare.pl", '--ignore-signs', @ARGV;
die "Can't execute $RealBin/compare.pl: $!\n";

#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin);

exec $^X, "$RealBin/mycomp.pl", '--ignore-signs', @ARGV;
die "Can't execute $RealBin/mycomp.pl: $!\n";

#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($RealBin);
use lib $RealBin;
use PhysicalOutput qw(is_physical_output);

@ARGV == 0 or die "Usage: clean-physical-output.pl\n";
opendir(my $directory, '.') or die "Can't read current directory: $!\n";
my @entries = grep { $_ ne '.' && $_ ne '..' } readdir($directory);
closedir($directory) or die "Can't finish reading current directory: $!\n";

for my $name (@entries) {
    next unless $name eq 'DONE' || is_physical_output($name);
    -d $name and die "Refusing to remove output-like directory $name\n";
    unlink($name) or die "Can't remove stale result $name: $!\n";
}

#!/usr/bin/env perl
# Compare all files in the directory passed as the command line parameter
# with the files in the current directory. Use mycompl.pl as diff tool.
# Rok Zitko, Nov 2019

use warnings;
use File::Basename;
use File::Spec;
my $scriptdir = File::Spec->rel2abs(dirname(__FILE__));

$refdir = shift;
if (!defined($refdir)) {
    die "Usage: compare.pl ref_dir";
}
-d $refdir or die "Directory $refdir does not exist.\n";

chomp(my $here = `pwd`);

chdir $refdir;
foreach (<*>) {
  $fn = $_;
  -f $fn or next;
  -f "$here/$fn" or die "Result $fn does not exist.";
  print "Comparing $fn\n";
  system("$scriptdir/mycomp-ignoresigns.pl", "$fn", "$here/$fn") == 0 or die "diff failed: $?";
}
print "OK!\n";

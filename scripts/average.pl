#!/usr/bin/env perl
# Spectral function averaging
# Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Oct 2008
# $Id: average.pl,v 1.1 2008/10/08 08:11:40 rzitko Exp rzitko $

use strict;
use warnings;

$|=1;

sub avg {
    my $Nz = shift; # number of z-values
    my $suffix = shift; # suffix of the spectral files (prefixed by 1..Nz)
    my $output = shift; # filename of the averaged spectrum
    print "$suffix -> $output (Nz=$Nz)\n";
    
    # Algorithm: just sum both columns and divide the sums by Nz at the
    # end. No attempt is made to check for the sanity of input files.
    my @sumx; 
    my @sumy;
    
    for (my $i = 1; $i <= $Nz; $i++) {
        my $filename = "$i$suffix";
	print "$i ";
	open(F, "<$filename") or die "Can't open $filename: $!\n";
	
	my $cnt = 0; # line counter
	while (<F>) {
	    chomp;
	    if (/^#/) { # skip comment lines
		next;
	    }
	
	    my $line = $_;
	    $line =~ s/^\s+(.*)/$1/; # strip leading white space
	    my @vals = split(/\s+/, $line); # split on whitespace
	    my $x = $vals[0];
	    my $y = $vals[1];
		
	    if ($i == 1) { 
		push(@sumx, $x); # push back in the first iteration
		push(@sumy, $y);
	    } else {
		$sumx[$cnt] += $x; # add in later iterations
		$sumy[$cnt] += $y;
	    }
	    $cnt++;
	}
	close(F);
    }
    
    my $len = @sumx;

    open(OUT, ">$output") or die "Can't open $output for writing: $!\n";
    
    for (my $l = 0; $l < $len; $l++) {
	$sumx[$l] /= $Nz;
	$sumy[$l] /= $Nz;
	print OUT "$sumx[$l] $sumy[$l]\n";
    }
    
    close(OUT);
    
    print "\n";
}

1;

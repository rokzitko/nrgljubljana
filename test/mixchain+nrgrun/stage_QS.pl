#!/usr/bin/env perl
# Stage a mixchain chain for the SIAM/QS template: rename the per-element files into the coefficient sets that the
# template's instantiate and matrix read.
#
#   stage_qs.pl [DIRECTORY]      (default: the working directory)
#
# For one channel the chain is a 1x1 matrix per site, so the mapping is a copy:
#
#   T11.dat -> xi1.dat      the hopping between shells, sites 0..Nmax, one value per line
#   E11.dat -> zeta1.dat    the on-site energy, the same layout
#   V11.dat -> theta1.dat   one number, squared: matrix reads theta and uses gammaPolCh = sqrt(theta/pi),
#                           while mixchain writes V = theta^(1/2), in the normalization of its input
#
# The layout is the one nrgchain writes, so for one channel the renaming is the whole of it. The U1 and SPSU2 cases
# map several elements each; see stage_U1.pl and stage_SPSU2.pl.

use strict;
use warnings;

my $directory = shift // ".";
chdir $directory or die "stage_qs: cannot enter $directory: $!\n";

# One column per line, as mixchain writes a real chain. A complex chain has the pair "Re Im" instead, which the QS
# template cannot take; that is caught here rather than silently truncated.
sub read_values {
    my $filename = shift;
    open(my $in, "<", $filename)
      or die "stage_qs: cannot read $filename: $!\nRun mixchain with discretization_files=true.\n";
    my @values;
    # An explicit loop variable, so that nothing is written into $_, which the caller may have aliased.
    while (defined(my $line = <$in>)) {
        next if $line =~ /^\s*(#|$)/;
        my @fields = split ' ', $line;
        die "stage_qs: $filename line $.: expected one real value, found '$line'\n" if @fields != 1;
        push @values, $fields[0];
    }
    close $in;
    die "stage_qs: $filename is empty.\n" unless @values;
    return @values;
}

sub write_values {
    my ($filename, @values) = @_;
    open(my $out, ">", $filename) or die "stage_qs: cannot write $filename: $!\n";
    printf $out "%.17g\n", $_ for @values;
    close $out or die "stage_qs: error writing $filename: $!\n";
}

my @xi   = read_values("T11.dat");
my @zeta = read_values("E11.dat");
my @v    = read_values("V11.dat");
die "stage_qs: V11.dat holds " . scalar(@v) . " values instead of one.\n" if @v != 1;
die "stage_qs: T11.dat and E11.dat hold different numbers of sites.\n" if @xi != @zeta;

write_values("xi1.dat",    @xi);
write_values("zeta1.dat",  @zeta);
write_values("theta1.dat", $v[0] * $v[0]);

printf "stage_qs: %d sites, theta1=%.17g\n", scalar(@xi), $v[0] * $v[0];

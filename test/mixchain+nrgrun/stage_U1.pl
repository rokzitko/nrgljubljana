#!/usr/bin/env perl
# Stage a 2x2 mixchain chain for the SIAM/U1 template with pol2x2=true, where the z block holds four coefficient
# sets per channel.
#
#   stage_u1.pl [DIRECTORY]      (default: the working directory)
#
# Channel 1 of Gamma is up and channel 2 is down. The four sets are UP, DOWN, UPDO and DOUP, in the order
# c++/coef.hpp maps them, and the terms they multiply are the macros of c++/symmetry/sym-U1-impl.hpp:
#
#   xi1 <- T11 (up,up)      xi2 <- T22 (down,down)      xi3 <- T12 (up,down)      xi4 <- T21 (down,up)
#   zeta1 <- E11            zeta2 <- E22                zeta3 <- E12              zeta4 <- E21
#   V{i}{j}1 <- V(i,j)      the four files matrix -s reads as coefV[i,j]
#
# nrg builds the Hamiltonian from zeta1, zeta2 and zeta4 alone: there is no DIAG_UPDO macro, so zeta3 enters only
# through the seed Hamiltonian of the template, which uses coefzeta[3,0]. For real data the two are the same number,
# and this script checks that rather than assuming it. A complex chain cannot go through this route at all, since
# matrix stores doubles.

use strict;
use warnings;

my $directory = shift // ".";
chdir $directory or die "stage_u1: cannot enter $directory: $!\n";

sub read_values {
    my $filename = shift;
    open(my $in, "<", $filename)
      or die "stage_u1: cannot read $filename: $!\nRun mixchain with discretization_files=true.\n";
    my @values;
    # An explicit loop variable: with while (<$in>) the line would go into $_, which is aliased to a read-only
    # literal when this is called from map.
    while (defined(my $line = <$in>)) {
        next if $line =~ /^\s*(#|$)/;
        my @fields = split ' ', $line;
        die "stage_u1: $filename line $.: expected one real value, found '$line'\n" if @fields != 1;
        push @values, $fields[0];
    }
    close $in;
    die "stage_u1: $filename is empty.\n" unless @values;
    return @values;
}

sub write_values {
    my ($filename, @values) = @_;
    open(my $out, ">", $filename) or die "stage_u1: cannot write $filename: $!\n";
    printf $out "%.17g\n", $_ for @values;
    close $out or die "stage_u1: error writing $filename: $!\n";
}

# T and E, one file per matrix element, sites 0..Nmax.
my %t = map { $_ => [ read_values("T$_.dat") ] } qw(11 22 12 21);
my %e = map { $_ => [ read_values("E$_.dat") ] } qw(11 22 12 21);
my %v = map { $_ => [ read_values("V$_.dat") ] } qw(11 22 12 21);

my $sites = scalar @{ $t{11} };
for my $element (qw(11 22 12 21)) {
    die "stage_u1: T$element.dat holds " . scalar(@{ $t{$element} }) . " sites instead of $sites.\n"
      if @{ $t{$element} } != $sites;
    die "stage_u1: E$element.dat holds " . scalar(@{ $e{$element} }) . " sites instead of $sites.\n"
      if @{ $e{$element} } != $sites;
    die "stage_u1: V$element.dat holds " . scalar(@{ $v{$element} }) . " values instead of one.\n"
      if @{ $v{$element} } != 1;
}

# E must be symmetric for the template route: zeta3 and zeta4 are one number each, and nrg and the seed Hamiltonian
# read different ones of them.
#
# Measured against the largest element of the whole chain, since the entries compared here vanish for a particle-hole
# symmetric band and a ratio of two roundings means nothing.
my $scale = 0;
for my $element (qw(11 22 12 21)) {
    for my $n (0 .. $sites - 1) {
        $scale = abs($e{$element}[$n]) if abs($e{$element}[$n]) > $scale;
        $scale = abs($t{$element}[$n]) if abs($t{$element}[$n]) > $scale;
    }
}
die "stage_u1: every coefficient of the chain is zero.\n" unless $scale > 0;

my $worst = 0;
for my $n (0 .. $sites - 1) {
    my $deviation = abs($e{12}[$n] - $e{21}[$n]) / $scale;
    $worst = $deviation if $deviation > $worst;
}
die sprintf("stage_u1: E12 and E21 differ by %.2e relative; the U1 template cannot express that.\n", $worst)
  if $worst > 1e-12;

write_values("xi1.dat", @{ $t{11} });
write_values("xi2.dat", @{ $t{22} });
write_values("xi3.dat", @{ $t{12} });
write_values("xi4.dat", @{ $t{21} });

write_values("zeta1.dat", @{ $e{11} });
write_values("zeta2.dat", @{ $e{22} });
write_values("zeta3.dat", @{ $e{12} });
write_values("zeta4.dat", @{ $e{21} });

# coefV[i,j] for the one channel: V{i}{j}{ch}.dat.
write_values("V${_}1.dat", $v{$_}[0]) for qw(11 22 12 21);

printf "stage_u1: %d sites, V = [[%.17g, %.17g], [%.17g, %.17g]], E off-diagonal asymmetry %.2e\n",
  $sites, $v{11}[0], $v{12}[0], $v{21}[0], $v{22}[0], $worst;

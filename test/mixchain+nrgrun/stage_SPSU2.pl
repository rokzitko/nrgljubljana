#!/usr/bin/env perl
# Stage a Nambu-gauge mixchain chain for the SIAM/SPSU2 template.
#
#   stage_spsu2.pl [--flip-scdelta] [DIRECTORY]      (default: the working directory)
#
# Channel 1 of Gamma is the particle component and channel 2 the hole component, psi = (c_up, c_down^dagger). From
# the macros of c++/symmetry/sym-SPSU2-impl.hpp, SPSU2 keeps one number per site from each block of the Nambu
# structure and reconstructs the rest:
#
#   xi1      <- T11   normal hopping, the diagonal of T
#   sckappa1 <- T12   anomalous hopping, the off-diagonal of T
#   zeta1    <- E11   on-site energy, the diagonal of E
#   scdelta1 <- E12   on-site pairing, the off-diagonal of E
#   V{i}{j}1 <- V(i,j)  the seed Hamiltonian of the template reads coefV[1,1] and coefV[1,2]
#
# That reconstruction only holds in the Nambu gauge, so the chain must come from mixchain with chain_gauge=nambu:
# E(2,2) = -E(1,1) and T(2,2) = -conj(T(1,1)). This script checks both rather than trusting the gauge.
#
# The mapping is structural: every entry goes through as it stands. What the sign of scdelta then is follows from the
# input, since the anomalous part of Gamma carries a Nambu basis phase -- psi = (c_up, c_down^dagger) against
# (c_up, -c_down^dagger) turns it over, and E12 with it.
#
# Which of the two phases matches nrg's scdelta is what this case measured: the one make_input.sh writes, where the
# anomalous part is negative above the gap. With the other one, <pair_d> comes out equal and opposite to the
# reference. That is a statement about the input, not about this script.
#
# --flip-scdelta writes -E12 instead. The test does not use it; it is how the measurement was made, by running both
# signs against a bath whose physics is known, and it stays so that it can be repeated.

use strict;
use warnings;

my $flip = 0;
my @rest;
for my $argument (@ARGV) {
    if ($argument eq "--flip-scdelta") { $flip = 1; next; }
    push @rest, $argument;
}
my $directory = $rest[0] // ".";
chdir $directory or die "stage_spsu2: cannot enter $directory: $!\n";

sub read_values {
    my $filename = shift;
    open(my $in, "<", $filename)
      or die "stage_spsu2: cannot read $filename: $!\nRun mixchain with discretization_files=true.\n";
    my @values;
    while (defined(my $line = <$in>)) {
        next if $line =~ /^\s*(#|$)/;
        my @fields = split ' ', $line;
        die "stage_spsu2: $filename line $.: expected one real value, found '$line'\n" if @fields != 1;
        push @values, $fields[0];
    }
    close $in;
    die "stage_spsu2: $filename is empty.\n" unless @values;
    return @values;
}

sub write_values {
    my ($filename, @values) = @_;
    open(my $out, ">", $filename) or die "stage_spsu2: cannot write $filename: $!\n";
    printf $out "%.17g\n", $_ for @values;
    close $out or die "stage_spsu2: error writing $filename: $!\n";
}

my %t = map { $_ => [ read_values("T$_.dat") ] } qw(11 22 12 21);
my %e = map { $_ => [ read_values("E$_.dat") ] } qw(11 22 12 21);
my %v = map { $_ => [ read_values("V$_.dat") ] } qw(11 22 12 21);
my $sites = scalar @{ $t{11} };

# The Nambu structure the template's reconstruction assumes. Real data here, so the conjugation is the identity.
#
# The deviations are measured against the largest element of the whole chain, not against the site they occur on. At
# particle-hole symmetry the on-site energies vanish, so E(1,1) and E(2,2) are rounding, and their ratio says nothing
# at all: two numbers of order 1e-200 differ by a factor of two as easily as not.
my $scale = 0;
for my $element (qw(11 22 12 21)) {
    for my $n (0 .. $sites - 1) {
        $scale = abs($e{$element}[$n]) if abs($e{$element}[$n]) > $scale;
        $scale = abs($t{$element}[$n]) if abs($t{$element}[$n]) > $scale;
    }
}
die "stage_spsu2: every coefficient of the chain is zero.\n" unless $scale > 0;

my ($worst_e, $worst_t) = (0, 0);
for my $n (0 .. $sites - 1) {
    my $deviation_e = abs($e{22}[$n] + $e{11}[$n]) / $scale;
    my $deviation_t = abs($t{22}[$n] + $t{11}[$n]) / $scale;
    $worst_e = $deviation_e if $deviation_e > $worst_e;
    $worst_t = $deviation_t if $deviation_t > $worst_t;
}
die sprintf("stage_spsu2: E(2,2) departs from -E(1,1) by %.2e relative; is chain_gauge=nambu set?\n", $worst_e)
  if $worst_e > 1e-8;
die sprintf("stage_spsu2: T(2,2) departs from -T(1,1) by %.2e relative; is chain_gauge=nambu set?\n", $worst_t)
  if $worst_t > 1e-8;

write_values("xi1.dat",      @{ $t{11} });
write_values("sckappa1.dat", @{ $t{12} });
write_values("zeta1.dat",    @{ $e{11} });
write_values("scdelta1.dat", map { $flip ? -$_ : $_ } @{ $e{12} });
write_values("V${_}1.dat", $v{$_}[0]) for qw(11 22 12 21);

printf "stage_spsu2: %d sites%s, V11=%.17g V12=%.17g, scdelta[0]=%.17g, sckappa[0]=%.17g, Nambu deviation %.1e/%.1e\n",
  $sites, ($flip ? " (scdelta written as -E12, the phase the physics rejects)" : ""), $v{11}[0], $v{12}[0],
  ($flip ? -$e{12}[0] : $e{12}[0]), $t{12}[0], $worst_e, $worst_t;

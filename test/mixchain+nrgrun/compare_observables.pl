#!/usr/bin/env perl
# The impurity observables of several runs side by side, at the lowest temperature of each.
#
#   compare_observables.pl [--tol=0.05] [--ops=a,b,c] REFERENCE_DIR OTHER_DIR...
#
# Each directory holds the custom file of one run. The first is the reference; the others are reported relative to
# it. The exit status is 1 if any of them differs from the reference by more than the tolerance, or if the sign of
# pair_d disagrees.
#
# Runs that describe different baths are not expected to agree to rounding, so the tolerance is a statement about how
# far apart the physics of two routes may be, not about arithmetic.

use strict;
use warnings;

my $tolerance = 0.05;
my @observables = qw(pair_d n_d n_d^2 Himp Hhyb);
@ARGV = grep { !/^--tol=(.*)$/ or ($tolerance = $1, 0) } @ARGV;
@ARGV = grep { !/^--ops=(.*)$/ or (@observables = split(/,/, $1), 0) } @ARGV;
my @directories = @ARGV;
die "usage: compare_observables.pl [--tol=X] [--ops=a,b] REFERENCE_DIR OTHER_DIR...\n" if @directories < 2;

sub read_last_row {
    my $directory = shift;
    my $path = "$directory/custom";
    open(my $in, "<", $path) or die "compare_observables: cannot read $path: $!\n";
    my (@names, @last);
    while (defined(my $line = <$in>)) {
        chomp $line;
        if ($line =~ /^#/) {
            my @fields = split ' ', $line;
            shift @fields;
            @names = @fields if @fields && $fields[0] eq "T";
            next;
        }
        next unless $line =~ /\S/;
        @last = split ' ', $line;
    }
    close $in;
    die "compare_observables: $path has no column names.\n" unless @names;
    die "compare_observables: $path has no data rows.\n" unless @last;
    my %value;
    $value{ $names[$_] } = $last[$_] for 0 .. $#names;
    return \%value;
}

my %value = map { $_ => read_last_row($_) } @directories;
my @labels = map { my $l = $_; $l =~ s{/$}{}; $l =~ s{.*/}{}; $l } @directories;

printf "%-10s", "observable";
printf " %20s", $_ for @labels;
print "\n";
for my $observable (@observables) {
    printf "%-10s", $observable;
    printf " %20s", (defined $value{$_}{$observable} ? $value{$_}{$observable} : "-") for @directories;
    print "\n";
}
printf "%-10s", "T";
printf " %20s", $value{$_}{T} for @directories;
print "\n\n";

my $reference = $directories[0];
my $status = 0;
for my $directory (@directories[1 .. $#directories]) {
    for my $observable (@observables) {
        my ($a, $b) = ($value{$reference}{$observable}, $value{$directory}{$observable});
        next unless defined $a && defined $b;
        my $scale = abs($a) > abs($b) ? abs($a) : abs($b);
        next unless $scale > 1e-12;
        my $deviation = abs($a - $b) / $scale;
        if ($observable eq "pair_d" && $a * $b < 0) {
            printf "SIGN     %-12s %-10s %.10g against %.10g\n", $directory, $observable, $b, $a;
            $status = 1;
            next;
        }
        if ($deviation > $tolerance) {
            printf "ABOVE    %-12s %-10s %.3f%% from the reference\n", $directory, $observable, 100 * $deviation;
            $status = 1;
        } else {
            printf "ok       %-12s %-10s %.3f%%\n", $directory, $observable, 100 * $deviation;
        }
    }
}
exit $status;

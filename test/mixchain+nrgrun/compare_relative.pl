#!/usr/bin/env perl
# Compare the numbers in the same files of two runs, relatively.
#
#   compare.pl [--tol=1e-10] [--floor=0] DIR_A DIR_B FILE...
#
# Each file is read from both directories and compared token by token. Numeric tokens are compared relatively,
# |a-b| / max(|a|,|b|), since chain coefficients fall like Lambda^(-n/2) and an absolute tolerance says nothing
# about the late sites. Non-numeric tokens (the headers of td and custom) must match exactly. The largest relative
# difference of each file is printed, and the exit status is 1 if any file is missing, has a different shape, or
# exceeds the tolerance.
#
# --floor is the magnitude below which a value counts as zero: where both sides are under it, they are equal
# whatever their ratio. A quantity that vanishes by symmetry -- zeta at particle-hole symmetry, <S_x> without a
# transverse field -- comes out as rounding on one side and as an exact zero on the other, and comparing those two
# relatively says nothing. It is off by default, so a case has to name the scale on which its numbers are zero.

use strict;
use warnings;

my $tolerance = 1e-10;
my $floor     = 0;
@ARGV = grep { !/^--tol=(.*)$/   or ($tolerance = $1, 0) } @ARGV;
@ARGV = grep { !/^--floor=(.*)$/ or ($floor     = $1, 0) } @ARGV;
my ($dir_a, $dir_b, @files) = @ARGV;
die "usage: compare.pl [--tol=X] DIR_A DIR_B FILE...\n" unless defined $dir_b && @files;

sub slurp {
    my $path = shift;
    open(my $in, "<", $path) or return undef;
    my @lines = <$in>;
    close $in;
    chomp @lines;
    return \@lines;
}

# A token is numeric when it parses as a number in full; anything else is a label.
sub numeric { return $_[0] =~ /^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$/; }

my $status = 0;
for my $file (@files) {
    my $a = slurp("$dir_a/$file");
    my $b = slurp("$dir_b/$file");
    if (!defined $a || !defined $b) {
        printf "MISSING  %-12s %s\n", $file, !defined $a ? "$dir_a/$file" : "$dir_b/$file";
        $status = 1;
        next;
    }
    if (@$a != @$b) {
        printf "DIFFERS  %-12s %d lines vs %d\n", $file, scalar @$a, scalar @$b;
        $status = 1;
        next;
    }

    my ($worst, $where, $failed) = (0, "", 0);
    for my $line (0 .. $#$a) {
        my @ta = split ' ', $a->[$line];
        my @tb = split ' ', $b->[$line];
        if (@ta != @tb) {
            printf "DIFFERS  %-12s line %d: %d columns vs %d\n", $file, $line + 1, scalar @ta, scalar @tb;
            $failed = 1;
            last;
        }
        for my $column (0 .. $#ta) {
            my ($x, $y) = ($ta[$column], $tb[$column]);
            if (!numeric($x) || !numeric($y)) {
                next if $x eq $y;
                printf "DIFFERS  %-12s line %d column %d: '%s' vs '%s'\n", $file, $line + 1, $column + 1, $x, $y;
                $failed = 1;
                next;
            }
            next if $x == $y;
            next if abs($x) <= $floor && abs($y) <= $floor;
            my $scale = abs($x) > abs($y) ? abs($x) : abs($y);
            my $deviation = $scale > 0 ? abs($x - $y) / $scale : 0;
            if ($deviation > $worst) {
                $worst = $deviation;
                $where = sprintf "line %d column %d: %.17g vs %.17g", $line + 1, $column + 1, $x, $y;
            }
        }
    }

    if ($failed) { $status = 1; next; }
    if ($worst > $tolerance) {
        printf "ABOVE    %-12s %.2e at %s\n", $file, $worst, $where;
        $status = 1;
    } else {
        printf "ok       %-12s largest relative difference %.2e\n", $file, $worst;
    }
}

exit $status;

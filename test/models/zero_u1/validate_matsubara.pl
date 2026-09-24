#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw(isfinite);

@ARGV == 2 or die "Usage: $0 REFERENCE_DIR ACTUAL_DIR\n";
my ($reference, $actual) = @ARGV;
my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;

sub read_table {
    my ($path) = @_;
    open my $in, '<', $path or die "$path: $!\n";
    my @rows;
    while (<$in>) {
        next if /^\s*$/;
        my @fields = split ' ';
        die "$path: expected three numeric columns\n" unless @fields == 3;
        for my $field (@fields) {
            die "$path: nonfinite or invalid number '$field'\n"
                unless $field =~ /\A$number\z/ && isfinite(0 + $field);
            my $mantissa = $field;
            $mantissa =~ s/[eE].*//;
            die "$path: unrepresentable number '$field'\n" if 0 + $field == 0 && $mantissa =~ /[1-9]/;
        }
        push @rows, [map { 0 + $_ } @fields];
    }
    close $in or die "$path: $!\n";
    die "$path: expected five Matsubara records\n" unless @rows == 5;
    return \@rows;
}

sub matches {
    my ($left, $right, $atol, $rtol) = @_;
    my $scale = abs($left) > abs($right) ? abs($left) : abs($right);
    return abs($left - $right) <= $atol + $rtol * $scale;
}

# Independent flat Z-discretization at Lambda=2, z=1. Each uncoupled
# channel/spin is a three-site chain; keep=5000 retains all 4096 states.
my $t0_squared = 1 / (7 * log(2)**2);
my $t1_squared = 18 / (217 * log(2)**2);
my $pi = 4 * atan2(1, 1);
for my $spin (qw(u d)) {
    my $name = "spec_FDMmats_dens_A_f_$spin-A_f_$spin.dat";
    my $expected = read_table("$reference/$name");
    my $observed = read_table("$actual/$name");
    for my $n (0 .. 4) {
        my ($frequency, $real, $imaginary) = @{$observed->[$n]};
        # Same matched-record rules as compare.pl --strict for spectral output.
        die "$name: record $n frequency differs from reference\n"
            unless matches($frequency, $expected->[$n][0], 1e-12, 1e-5);
        die "$name: record $n imaginary part differs from reference\n"
            unless matches($imaginary, $expected->[$n][2], 1e-12, 2e-2);
        die "$name: record $n real part exceeds the analytic-zero bound 1e-12\n" if abs($real) > 1e-12;

        my $omega = (2 * $n + 1) * $pi * 0.01;
        my $analytic_imaginary = -($omega**2 + $t1_squared) / ($omega * ($omega**2 + $t0_squared + $t1_squared));
        die "$name: record $n is not on the analytic Matsubara grid\n"
            if abs($frequency - $omega) > 5e-10 * $omega;
        die "$name: record $n imaginary part differs from the analytic three-site resolvent\n"
            if abs($imaginary - $analytic_imaginary) > 5e-10 * abs($analytic_imaginary);
    }
    print "Validated $name (analytic zero and three-site resolvent; matched reference grid/imaginary part)\n";
}

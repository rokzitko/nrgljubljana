#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Test::More;

my $work = tempdir(CLEANUP => 1);
make_path("$work/ref", "$work/actual");
# Neutral static data for this known free model, not either production backend.
my @baseline = map { [split ' '] } split /\n/, <<'DATA';
0.03141592654 0 -11.73521686
0.09424777961 0 -4.022186552
0.1570796327 0 -2.539502447
0.2199114858 0 -1.938837019
0.2827433388 0 -1.624498748
DATA

sub write_table {
    my ($directory, $spin, $rows) = @_;
    my $path = "$work/$directory/spec_FDMmats_dens_A_f_$spin-A_f_$spin.dat";
    open my $out, '>', $path or die "$path: $!";
    print {$out} join(' ', @$_), "\n" for @$rows;
    close $out or die "$path: $!";
}

sub validate {
    my ($label, $spin, $directory, $mutate, $diagnostic) = @_;
    for my $dir (qw(ref actual)) {
        write_table($dir, $_, \@baseline) for qw(u d);
    }
    my @rows = map { [@$_] } @baseline;
    $mutate->(\@rows);
    write_table($directory, $spin, \@rows);
    my $pid = fork();
    die "fork: $!" unless defined $pid;
    if (!$pid) {
        open STDOUT, '>', "$work/validator.log" or die $!;
        open STDERR, '>&', \*STDOUT or die $!;
        exec $^X, "$RealBin/validate_matsubara.pl", "$work/ref", "$work/actual";
        die "exec: $!";
    }
    waitpid($pid, 0);
    my $status = $?;
    open my $log, '<', "$work/validator.log" or die $!;
    my $output = do { local $/; <$log> };
    if (defined $diagnostic) {
        isnt($status, 0, "$spin rejects $label");
        like($output, $diagnostic, "$spin diagnoses $label");
    } else {
        is($status, 0, "$spin accepts $label") or diag $output;
    }
}

for my $spin (qw(u d)) {
    validate('analytic static data', $spin, 'actual', sub {}, undef);
    validate('unchanged real absolute boundary', $spin, 'actual', sub { $_[0][0][1] = '1e-12'; $_[0][1][1] = '-1e-12' }, undef);
    validate('historical reference real roundoff', $spin, 'ref', sub { $_[0][0][1] = '-9.594024077e-13' }, undef);
    for my $real ('1.000001e-12', '-1.000001e-12') {
        validate("nonzero real $real", $spin, 'actual', sub { $_[0][0][1] = $real }, qr/analytic-zero bound/);
    }
    validate('wrong grid', $spin, 'actual', sub { $_[0][0][0] *= 1.1 }, qr/frequency differs from reference/);
    validate('small grid error', $spin, 'actual', sub { $_[0][0][0] *= 1 + 1e-7 }, qr/analytic Matsubara grid/);
    validate('wrong imaginary part', $spin, 'actual', sub { $_[0][0][2] *= 1.1 }, qr/imaginary part differs from reference/);
    validate('small imaginary error', $spin, 'actual', sub { $_[0][0][2] *= 1 + 1e-7 }, qr/analytic three-site resolvent/);
    validate('wrong reference grid', $spin, 'ref', sub { $_[0][0][0] *= 1.1 }, qr/frequency differs from reference/);
    validate('wrong reference imaginary part', $spin, 'ref', sub { $_[0][0][2] *= 1.1 }, qr/imaginary part differs from reference/);
    validate('missing row', $spin, 'actual', sub { pop @{$_[0]} }, qr/expected five/);
    validate('extra row', $spin, 'actual', sub { push @{$_[0]}, [@{$_[0][0]}] }, qr/expected five/);
    validate('reordered rows', $spin, 'actual', sub { @{$_[0]}[0, 1] = @{$_[0]}[1, 0] }, qr/frequency differs from reference/);
    validate('wrong column count', $spin, 'actual', sub { push @{$_[0][0]}, 0 }, qr/expected three/);
    for my $column (0 .. 2) {
        for my $invalid ('NaN', 'Inf', '1e9999') {
            validate("nonfinite column $column ($invalid)", $spin, 'actual', sub { $_[0][0][$column] = $invalid }, qr/nonfinite or invalid/);
        }
    }
    validate('underflowed nonzero value', $spin, 'actual', sub { $_[0][0][1] = '1e-9999' }, qr/unrepresentable/);
}
done_testing();

#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Temp qw(tempdir);
use POSIX qw(isfinite);
use Scalar::Util qw(looks_like_number);
my $tools = abs_path(shift @ARGV // die "Usage: check.pl BUILD_TOOLS\n");
my $backend = $ENV{CHAIN_TEST_BACKEND} // die "Missing backend\n";
my $method = $backend eq 'legacy' ? 'lanczos' : $backend eq 'rkpw' ? 'rkpw' : die "Invalid backend\n";
my $work = tempdir('star-ranges-XXXXXX', DIR => '.');
chdir $work or die $!;
sub write_text {
    open my $out, '>', $_[0] or die "$!: $_[0]";
    print {$out} $_[1]; close $out or die $!;
}
sub text { open my $in, '<', $_[0] or die "$!: $_[0]"; local $/; return <$in>; }
sub parameters {
    my %p = (tridiag_method => $method, Lambda => 2, z => 0.5, Nmax => 4, mMAX => 4,
             band => 'flat', adapt => 'false', hardgap => 'true', boundary => 0.25,
             preccpp => 4096, xmax => 6, outputstep => 0.125, @_);
    write_text('param', "[param]\n" . join('', map { "$_=$p{$_}\n" } sort keys %p));
}
sub run {
    my ($label, $tool, $error, @args) = @_;
    my $pid = fork(); die $! unless defined $pid;
    if (!$pid) {
        open STDOUT, '>', "$label.out" or die $!;
        open STDERR, '>', "$label.err" or die $!;
        exec {"$tools/$tool"} "$tools/$tool", @args; die $!;
    }
    waitpid($pid, 0);
    my $status = $?;
    die "$label failed ($status): " . text("$label.err") if ($status & 127) || ($error ? $status == 0 : $status != 0);
    if ($error) { die "$label: wrong error: " . text("$label.err") unless text("$label.err") =~ $error; }
}
sub read_values {
    my @v = split ' ', text($_[0]);
    die "Empty/nonfinite $_[0]" if !@v || grep { !looks_like_number($_) || !isfinite(0+$_) } @v;
    return \@v;
}
sub near {
    my ($label, $a, $b, $tol) = @_;
    die "$label: $a != $b" if abs($a-$b) > $tol;
}
sub compare_chains {
    my ($xi, $zeta) = @_;
    my $x = read_values('xi.dat'); my $e = read_values('zeta.dat');
    die 'Wrong chain length' unless @$x == @$xi && @$e == @$zeta;
    for my $n (0..$#$xi) {
        near("xi[$n]", $x->[$n]/$xi->[$n], 1, 2e-12);
        near("zeta[$n]", ($e->[$n]-$zeta->[$n])/$xi->[$n], 0, 2e-12);
    }
}
sub flat_energy {
    my ($x) = @_;
    my $h = $x < 2 ? 2-$x+(1-2**(1-$x))/log(2) : 2**(2-$x)/(2*log(2));
    return 0.25+0.75*$h;
}

for my $interpolation (qw(linear steffen)) {
    parameters(density_interpolation => $interpolation);
    run("flat-$interpolation", 'nrgchain', undef);
    my $xi = read_values('xi.dat'); my $zeta = read_values('zeta.dat');
    my $theta = 0.75*(1-2**(-4.5));
    near('retained theta', read_values('theta.dat')->[0], $theta, 2e-15);
    run('flat-save', 'nrgchain', undef, 's');
    my $de = read_values('de_pos.dat'); my $du = read_values('du_pos.dat');
    for my $m (0..4) {
        my $upper = $m == 0 ? 1 : 0.25+0.75*2**(0.5-$m);
        my $lower = 0.25+0.75*2**(-0.5-$m);
        near('analytic representative', $de->[$m], flat_energy(1.5+$m), 2e-15);
        near('analytic amplitude', $du->[$m], sqrt(($upper-$lower)/(2*$theta)), 2e-15);
    }
    # Actual adapt export, no GSOL, same normalized density as the analytic branch.
    write_text('Delta.dat', "-1 0.5\n-0.5 0.5\n-0.25 0.5\n0.25 0.5\n0.5 0.5\n1 0.5\n");
    parameters(band => 'adapt', density_interpolation => $interpolation);
    run('adapt-pos', 'adapt', undef, '--integral', 'P');
    run('adapt-neg', 'adapt', undef, '--integral', 'N');
    run('table-gap', 'nrgchain', undef, 's');
    run('table-gap-chain', 'nrgchain', undef);
    compare_chains($xi, $zeta);
    run('gap-load', 'nrgchain', undef, 'l');
    compare_chains($xi, $zeta);
    parameters(band => 'adapt', bandrescale => 2, nrgchain_tables_load => 'true');
    run('gap-instantiate', 'instantiate', undef, '--wilson-only');
    my $scaled = read_values('xi1.dat');
    for my $n (0..$#$xi) { near('physical hopping scale', $scaled->[$n]/$xi->[$n], 2, 4e-12); }
    near('theta not scaled twice', read_values('theta1.dat')->[0], $theta, 2e-15);
}

# Off-grid x must interpolate nodal energies, not exponentially growing f.
parameters(band => 'adapt', z => 0.37);
run('off-grid', 'nrgchain', undef, 's');
my @rows = map { [split ' '] } split /\n/, text('FSOL.dat');
my $energies = read_values('de_pos.dat');
for my $m (0..4) {
    my $x = 1.37+$m;
    my $j = 0; ++$j while $rows[$j+1][0] < $x;
    my ($x0,$f0) = @{$rows[$j]}; my ($x1,$f1) = @{$rows[$j+1]};
    my $e = $f0*2**(2-$x0) + ($x-$x0)/($x1-$x0)*($f1*2**(2-$x1)-$f0*2**(2-$x0));
    near('off-grid energy interpolation', $energies->[$m], $e, 2e-15);
}
my %before = map { $_ => text($_) } qw(theta.dat xi.dat zeta.dat de_pos.dat de_neg.dat du_pos.dat du_neg.dat);
sub preserved { for my $name (keys %before) { die "Failure modified $name" unless text($name) eq $before{$name}; } }
parameters(band => 'adapt', mMAX => 20, xmax => 999);
run('short-table', 'nrgchain', qr/does not cover/, 's'); preserved();
my $negative = text('FSOLNEG.dat');
write_text('FSOLNEG.dat', "1 0.5\n2 0.8\n");
parameters(band => 'adapt'); run('short-negative', 'nrgchain', qr/FSOLNEG.dat.*does not cover/); preserved();
write_text('FSOLNEG.dat', "2 0.8\n1 0.5\n6 10\n");
run('unordered-negative', 'nrgchain', qr/increasing abscissas/); preserved();
write_text('FSOLNEG.dat', "1 0.01\n6 0.01\n");
run('bad-gap-energy', 'nrgchain', qr/outside its hard-gap shell/); preserved();
write_text('FSOLNEG.dat', $negative);
parameters(mMAX => 200); run('collapsed-gap', 'nrgchain', qr/collapsed/); preserved();
parameters(adapt => 'true'); run('adaptive-gap', 'nrgchain', qr/require adapt=false/); preserved();
parameters(adapt => 'true', hardgap => 'false'); run('adaptive-flat', 'nrgchain', qr/require adapt=false/); preserved();
parameters(boundary => 1); run('bad-boundary', 'nrgchain', qr/boundary/); preserved();
parameters(band => 'adapt', max_abs => 0.6, xmax => 10);
my $positive = text('FSOL.dat');
run('short-export', 'adapt', qr/max_abs before the requested extent/, '--integral', 'P');
die 'Failed export clobbered FSOL' unless text('FSOL.dat') eq $positive;
run('gap-ode', 'adapt', qr/integral/, 'P');

# Deep masses may underflow even though amplitudes are representable. No ODE
# approximation enters this oracle: f is the exact flat Z prefactor for x>=2.
for my $interpolation (qw(linear steffen)) {
    my ($normal_xi, $normal_zeta);
    for my $factor (1, 1e-200) {
        write_text('Delta.dat', join('', map { "$_ $factor\n" } (-1,-0.5,-0.25,0.25,0.5,1)));
        my $c = 1/(2*log(2));
        write_text('FSOL.dat', "1 0.5\n2 $c\n"); write_text('FSOLNEG.dat', text('FSOL.dat'));
        parameters(band => 'adapt', hardgap => 'false', z => 1, mMAX => 800, Nmax => 100,
                   density_interpolation => $interpolation);
        run("deep-save-$interpolation-$factor", 'nrgchain', undef, 's');
        run("deep-chain-$interpolation-$factor", 'nrgchain', undef, 'l');
        my $du = read_values('du_pos.dat');
        die 'Wrong deep star length' unless @$du == 801;
        for my $m (0..800) {
            # Normal masses retain the old cumulative-integral compatibility
            # path (1024 epsilon budget). Underflow recovery uses local means.
            my $tol = $factor*2**(-$m-1) < 2.2250738585072014e-308 ? 5e-14 : 5e-13;
            near("deep amplitude $interpolation rho=$factor m=$m", $du->[$m]/(0.5*2**(-$m/2)), 1, $tol);
        }
        near('physical total density scale', read_values('theta.dat')->[0]/(2*$factor), 1, 3e-14);
        my $xi = read_values('xi.dat');
        for my $n (0..100) {
            my $expected = $c*2**(-$n/2)*(1-2**(-$n-1))/sqrt((1-2**(-2*$n-1))*(1-2**(-2*$n-3)));
            near('deep analytic hopping', $xi->[$n]/$expected, 1, 5e-12);
        }
        if ($factor == 1) { $normal_xi = $xi; $normal_zeta = read_values('zeta.dat'); }
        else { compare_chains($normal_xi, $normal_zeta); }
    }
}
# A linear pseudogap reaches underflowed shell masses without a density floor
# anywhere in the retained mesh. The expected amplitude is analytic.
write_text('Delta.dat', "-1 1\n-0.5 0.5\n-2.2250738585072014e-308 2.2250738585072014e-308\n2.2250738585072014e-308 2.2250738585072014e-308\n0.5 0.5\n1 1\n");
my $c = sqrt(3/(8*log(2)));
write_text('FSOL.dat', "1 0.5\n2 $c\n"); write_text('FSOLNEG.dat', text('FSOL.dat'));
parameters(band => 'adapt', hardgap => 'false', z => 1, mMAX => 800, Nmax => 100);
run('pseudogap-save', 'nrgchain', undef, 's');
run('pseudogap-chain', 'nrgchain', undef, 'l');
my $du = read_values('du_pos.dat');
for my $m (0..800) { near('pseudogap normalized amplitude', $du->[$m]/(sqrt(3/8)*2**(-$m)), 1, 5e-14); }
near('pseudogap first hopping', read_values('xi.dat')->[0], $c/sqrt(1.25), 2e-14);
# Each branch retains 1.5 least-subnormal units: normalize before rounding the
# two branch masses separately. The exported total is exactly three units.
write_text('Delta.dat', join('', map { "$_ 9.8813129168249309e-324\n" } (-1,-0.5,-0.25,0.25,0.5,1)));
my $flat_c = 1/(2*log(2));
write_text('FSOL.dat', "1 0.5\n2 $flat_c\n"); write_text('FSOLNEG.dat', text('FSOL.dat'));
parameters(band => 'adapt', hardgap => 'false', z => 1, mMAX => 1, Nmax => 1);
run('subnormal-theta', 'nrgchain', undef, 's');
my $subnormal_du = read_values('du_pos.dat');
near('subnormal theta amplitude 0', $subnormal_du->[0], 1/sqrt(3), 5e-15);
near('subnormal theta amplitude 1', $subnormal_du->[1], 1/sqrt(6), 5e-15);
near('subnormal combined theta', read_values('theta.dat')->[0]/4.9406564584124654e-324, 3, 0);
run('subnormal-theta-load', 'nrgchain', undef, 'l');
print "Scalar star range tests passed for $backend (artifacts in $work).\n";

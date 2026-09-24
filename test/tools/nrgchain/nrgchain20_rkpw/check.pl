#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Copy qw(copy);
use File::Temp qw(tempdir);
use POSIX qw(isfinite);
use Scalar::Util qw(looks_like_number);

my $tools = abs_path(shift @ARGV // die "Usage: check.pl BUILD_TOOLS\n");
my %base;
open my $param, '<', 'param' or die "param: $!";
while (<$param>) {
    $base{$1} = $2 if /^([^=\s]+)=(\S+)\s*$/;
}
close $param;
my $backend = $ENV{CHAIN_TEST_BACKEND} // die "CHAIN_TEST_BACKEND is required\n";
die "Unknown backend: $backend\n" unless $backend eq 'legacy' || $backend eq 'rkpw';
my $method = $backend eq 'legacy' ? 'lanczos' : 'rkpw';
die "Staged param does not select $method\n" unless ($base{tridiag_method} // '') eq $method;
my $work = tempdir("$backend-XXXXXX", DIR => '.');
copy('Delta.dat', "$work/Delta.dat") or die "copy Delta.dat: $!";
chdir $work or die "chdir $work: $!";

sub write_text {
    my ($file, $text) = @_;
    open my $out, '>', $file or die "$file: $!";
    print {$out} $text;
    close $out or die "$file: $!";
}

sub text {
    my ($file) = @_;
    open my $in, '<', $file or die "$file: $!";
    local $/;
    return <$in>;
}

sub parameters {
    my %values = (%base, @_);
    write_text('current.param', "[param]\n" . join('', map { "$_=$values{$_}\n" } sort keys %values));
}

sub run {
    my ($label, $tool, $failure, @args) = @_;
    unless ($failure) {
        my @outputs = $tool eq 'nrgchain' ? qw(xi.dat zeta.dat)
                    : $tool eq 'instantiate' ? qw(theta1.dat xi1.dat zeta1.dat) : ();
        unlink($_) or die "unlink $_: $!" for grep { -e $_ } @outputs;
    }
    my $pid = fork();
    die "fork: $!" unless defined $pid;
    if (!$pid) {
        open STDOUT, '>', "$label.out" or die $!;
        open STDERR, '>', "$label.err" or die $!;
        exec {"$tools/$tool"} "$tools/$tool", @args;
        die "exec $tool: $!";
    }
    waitpid($pid, 0);
    my $status = $?;
    my $err = text("$label.err");
    die "$label: unexpected exit status $status\n$err\nSee $work/$label.out\n"
        if ($status & 127) || ($failure ? $status == 0 : $status != 0);
    return $err;
}

sub read_values {
    my ($file) = @_;
    my @values = split ' ', text($file);
    die "$file: empty or nonfinite output\n"
        if !@values || grep { !looks_like_number($_) || !isfinite(0 + $_) } @values;
    return \@values;
}

sub chain {
    my ($suffix, $count) = @_;
    my %chain = map { $_ => read_values("$_$suffix.dat") } qw(theta xi zeta);
    die "Wrong coefficient count\n" unless @{$chain{theta}} == 1 && @{$chain{xi}} == $count && @{$chain{zeta}} == $count;
    return \%chain;
}

sub near {
    my ($label, $expected, $actual) = @_;
    die "$label: count mismatch\n" unless @$expected == @$actual;
    for my $i (0 .. $#$expected) {
        my ($a, $b) = ($expected->[$i], $actual->[$i]);
        my $scale = abs($a) > abs($b) ? abs($a) : abs($b);
        # Standard tool-comparator tolerances; no relaxed golden-reference overrides.
        die "$label\[$i\]: $a != $b\n" if abs($a - $b) > 1e-12 + 1e-5 * $scale;
    }
}

sub compare_chain {
    my ($label, $expected, $actual, $rescale) = @_;
    near("$label theta", $expected->{theta}, $actual->{theta});
    die "$label: count mismatch\n" unless @{$expected->{xi}} == @{$actual->{xi}}
        && @{$expected->{zeta}} == @{$actual->{zeta}};
    for my $i (0 .. $#{$expected->{xi}}) {
        my $hopping = abs($expected->{xi}[$i]);
        if ($hopping == 0) {
            die "$label xi[$i]: terminal hopping is not exactly zero\n" unless $actual->{xi}[$i] == 0;
        } else {
            die "$label xi[$i]: relative error exceeds 2e-12\n"
                if abs(($actual->{xi}[$i] - $expected->{xi}[$i]) / $hopping) > 2e-12;
        }
        # Onsite errors are measured against the adjacent (un-rescalexi'd)
        # hopping, including the preceding link at an exactly terminal site.
        my $local = 0;
        for my $j ($i > 0 ? ($i - 1, $i) : ($i)) {
            my $link = abs($expected->{xi}[$j]);
            $link *= (1 - 1 / $base{Lambda}) / log($base{Lambda})
                     * $base{Lambda} ** (-$j / 2 + 1 - $base{z}) if $rescale;
            $local = $link if $link > $local;
        }
        my $difference = abs($actual->{zeta}[$i] - $expected->{zeta}[$i]);
        die "$label zeta[$i]: local-scale error exceeds 2e-12\n"
            if ($local == 0 ? $difference != 0 : $difference / $local > 2e-12);
    }
}

my $count = $base{Nmax} + 1;
for my $band (qw(flat adapt)) {
    parameters(band => $band);
    if ($band eq 'adapt') {
        run('adapt-positive', 'adapt', 0, '--integral', 'P', 'current.param');
        run('adapt-negative', 'adapt', 0, '--integral', 'N', 'current.param');
    }
    # All producers and oracles in this registration use only the selected backend.
    run("$band-save", 'nrgchain', 0, 's', 'current.param');
    my %saved = map { $_ => text("$_.dat") } qw(theta de_pos de_neg du_pos du_neg);
    my $err = run("$band-load", 'nrgchain', 0, '-v', 'l', 'current.param');
    die "Missing selected backend configuration\n" unless $err =~ /^  tridiag_method=\Q$method\E$/m;
    my $raw = chain('', $count);
    my %loaded_bytes = map { $_ => text("$_.dat") } qw(xi zeta);
    run("$band-load-repeat", 'nrgchain', 0, 'l', 'current.param');
    die "Repeated load changed the outputs\n"
        if grep { text("$_.dat") ne $loaded_bytes{$_} } qw(xi zeta);

    parameters(band => $band, $backend eq 'rkpw' ? (preccpp => 0) : ());
    run("$band-calculate", 'nrgchain', 0, 'current.param');
    compare_chain("$band calculate/save/load", $raw, chain('', $count));
    die "Calculation changed saved star or theta\n"
        if grep { text("$_.dat") ne $saved{$_} } keys %saved;
    run("$band-calculate-instantiate", 'instantiate', 0, '--wilson-only', '--param', 'current.param');
    compare_chain("$band calculated result arrays", $raw, chain('1', $count));

    for my $scale (1, 0.125, 8) {
        for my $rescale (0, 1) {
            my $label = "$band-scale$scale-rescale$rescale";
            my %settings = (band => $band, bandrescale => $scale, rescalexi => $rescale,
                             nrgchain_tables_load => 'true');
            parameters(%settings, $backend eq 'rkpw' ? (preccpp => 0) : ());
            $err = run("$label-nrgchain", 'nrgchain', 0, '-v', 'l', 'current.param');
            die "Missing selected backend configuration\n" unless $err =~ /^  tridiag_method=\Q$method\E$/m;
            die "Missing inactive RKPW precision diagnostic\n"
                if $backend eq 'rkpw' && $err !~ /gmp_precision=auto -> inactive/;
            my $scaled = chain('', $count);

            my @xi = map {
                my $factor = (1 - 1 / $base{Lambda}) / log($base{Lambda})
                             * $base{Lambda} ** (-$_ / 2 + 1 - $base{z});
                $raw->{xi}[$_] * $scale / ($rescale ? $factor : 1);
            } 0 .. $count - 1;
            near("$label hopping scaling", \@xi, $scaled->{xi});
            near("$label onsite scaling", [map { $_ * $scale } @{$raw->{zeta}}], $scaled->{zeta});
            near("$label unchanged theta", $raw->{theta}, $scaled->{theta});

            $err = run("$label-instantiate", 'instantiate', 0, '--wilson-only', '-v', '--param', 'current.param');
            die "Missing instantiate backend configuration\n" unless $err =~ /^  nrgchain.tridiag_method=\Q$method\E$/m;
            die "Missing inactive instantiate precision diagnostic\n"
                if $backend eq 'rkpw' && $err !~ /nrgchain.preccpp=auto -> inactive/;
            compare_chain("$label instantiate result arrays", $scaled, chain('1', $count), $rescale);
        }
    }
}

for my $tool (qw(nrgchain instantiate)) {
    my @args = $tool eq 'nrgchain' ? ('current.param') : ('--wilson-only', '--param', 'current.param');
    for my $method ('bogus', 'RKPW', '') {
        parameters(tridiag_method => $method);
        my $err = run("$tool-invalid-$method", $tool, 1, @args);
        die "Missing method diagnostic\n" unless $err =~ /tridiag_method.*expected lanczos or rkpw/;
    }
    for my $precision ('-1', 'no', '1.5', '2147483648') {
        parameters(preccpp => $precision);
        my $err = run("$tool-invalid-precision-$precision", $tool, 1, @args);
        die "Missing precision parsing diagnostic\n" unless $err =~ /preccpp/;
    }
    if ($backend eq 'legacy') {
        parameters(preccpp => 10);
        my $err = run("$tool-low-lanczos-precision", $tool, 1, @args);
        die "Missing legacy precision diagnostic\n" unless $err =~ /preccpp must be greater than 10/;
    }
}
parameters(tridiag_method => 'invalid');
my $err = run('save-invalid-method', 'nrgchain', 1, 's', 'current.param');
die "Save-only mode did not validate the method\n" unless $err =~ /tridiag_method/;

# Exact finite-support termination and scaled-output rejection are RKPW-specific contracts.
if ($backend eq 'legacy') {
    print "Scalar legacy integration passed (artifacts in $work).\n";
    exit 0;
}

# Eight rows reduce to two supported energies: duplicates combine, zero weights vanish.
write_text('de_pos.dat', "0.9\n0.9\n0.45\n0.2\n");
write_text('de_neg.dat', "0.7\n0.7\n0.3\n0.1\n");
write_text('du_pos.dat', "0.5\n0.5\n0\n0\n");
write_text('du_neg.dat', "0.5\n0.5\n0\n0\n");
write_text('theta.dat', "2.75\n");
parameters(mMAX => 3, Nmax => 1, nrgchain_tables_load => 'true');
run('finite-support', 'nrgchain', 0, 'l', 'current.param');
my $finite = chain('', 2);
near('finite hopping', [0.8, 0], $finite->{xi});
near('finite onsite', [0.1, 0.1], $finite->{zeta});
near('finite theta', [2.75], $finite->{theta});
die "Terminal hopping is not exactly zero\n" unless $finite->{xi}[-1] == 0;
parameters(mMAX => 3, Nmax => 1, preccpp => 0, nrgchain_tables_load => 'true');
run('finite-support-zero-precision', 'nrgchain', 0, 'l', 'current.param');
compare_chain('RKPW ignores preccpp', $finite, chain('', 2));
run('finite-support-instantiate', 'instantiate', 0, '--wilson-only', '--param', 'current.param');
compare_chain('finite support result arrays', $finite, chain('1', 2));

parameters(mMAX => 3, Nmax => 2, nrgchain_tables_load => 'true');
my %before = map { $_ => text("$_.dat") } qw(xi zeta xi1 zeta1 theta1);
for my $tool (qw(nrgchain instantiate)) {
    my @args = $tool eq 'nrgchain' ? ('l', 'current.param') : ('--wilson-only', '--param', 'current.param');
    $err = run("$tool-overlong", $tool, 1, @args);
    die "Missing finite-support diagnostic: $err\n" unless $err =~ /support/i;
    die "Failed calculation truncated coefficient outputs\n"
        if grep { text("$_.dat") ne $before{$_} } keys %before;
}

# Kernel results can be valid while bandrescale makes their output unrepresentable.
for my $case (
    ['hopping-overflow', '1.7976931348623157e308', 4, 2, 0.7071067811865476, 0.7071067811865476],
    ['onsite-overflow', '1e200', '1e200', 1, 1, 0],
    ['hopping-underflow', '4.9406564584124654e-324', 0.125, 0.125, 0.7071067811865476, 0.7071067811865476],
    ['onsite-underflow', '4.9406564584124654e-324', 1.125, 0.875, 0.7071067811865476, 0.7071067811865476],
) {
    my ($label, $scale, $ep, $em, $up, $um) = @$case;
    write_text('de_pos.dat', "$ep\n0.05\n");
    write_text('de_neg.dat', "$em\n0.05\n");
    write_text('du_pos.dat', "$up\n0\n");
    write_text('du_neg.dat', "$um\n0\n");
    parameters(mMAX => 1, Nmax => 0, bandrescale => $scale,
               nrgchain_tables_load => 'true');
    for my $tool (qw(nrgchain instantiate)) {
        my @args = $tool eq 'nrgchain' ? ('l', 'current.param') : ('--wilson-only', '--param', 'current.param');
        $err = run("$tool-$label", $tool, 1, @args);
        die "Missing scaled-coefficient diagnostic: $err\n" unless $err =~ /scaled coefficient.*nonfinite or underflowed/;
        die "Invalid scaling truncated coefficient outputs\n"
            if grep { text("$_.dat") ne $before{$_} } keys %before;
    }
}
print "Scalar RKPW integration passed (artifacts in $work).\n";

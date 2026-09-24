#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Spec;

my $root = File::Spec->rel2abs(shift @ARGV // die "Usage: chain-backend.t TEST_SOURCE\n");
my $nrgchain = shift @ARGV;
my $helper = "$root/chain-backend.pl";
my $work = tempdir(CLEANUP => 1);
mkdir "$work/source" or die $!;
mkdir "$work/run" or die $!;

sub write_text {
    my ($path, $text) = @_;
    open my $out, '>', $path or die $!;
    print {$out} $text;
    close $out or die $!;
}
sub read_text {
    open my $in, '<', $_[0] or die $!;
    local $/;
    return <$in>;
}
sub stage {
    my ($backend, $mode, $text, $succeeds) = @_;
    write_text("$work/source/param", $text);
    write_text("$work/run/param", 'sentinel');
    my $status = system($^X, $helper, '--backend', $backend, '--mode', $mode,
                        '--input', "$work/source/param", '--output', "$work/run/param");
    is($status == 0 ? 1 : 0, $succeeds, "$backend/$mode staging status");
    is(read_text("$work/source/param"), $text, 'source stays immutable');
    is(read_text("$work/run/param"), 'sentinel', 'failure preserves destination') unless $succeeds;
    return read_text("$work/run/param");
}

for my $backend (qw(legacy rkpw)) {
    my $method = $backend eq 'legacy' ? 'lanczos' : 'rkpw';
    my $tri = $backend eq 'legacy' ? 'old' : 'rkpw';
    my $text = "# fixture\n[param]\n U = !1/2\n[extra]\ntri=untouched\ntridiag_method=untouched\n";
    my $actual = stage($backend, 'initializer', $text, 1);
    like($actual, qr/tri=$tri\ntridiag_method=$method\n\[extra\]/, 'insert inside the correct block');
    like($actual, qr/ U = !1\/2\n/, 'expressions preserved');
    like($actual, qr/\[extra\]\ntri=untouched\ntridiag_method=untouched\n/, 'other blocks preserved');
    is(stage($backend, 'initializer', $actual, 1), $actual, 'idempotent');
    for my $existing (qw(old rkpw cpp none)) {
        $actual = stage($backend, 'initializer', "[param]\n tri = $existing # comment\ntridiag_method=other\n", 1);
        my $expected = $existing =~ /^(cpp|none)$/ ? $existing : $tri;
        like($actual, qr/^tri=$expected$/m, 'initializer dispatch preserved or selected');
        like($actual, qr/^tridiag_method=$method$/m, 'explicit backend');
    }
    for my $mode (qw(tool runtime)) {
        $actual = stage($backend, $mode, "[param]\ntri=manual\npreccpp=4000", 1);
        like($actual, qr/tri=manual\npreccpp=4000\ntridiag_method=$method\n/, 'manual mode and precision preserved');
    }
}
for my $text ("[extra]\ntri=old\n", "[param]\n[extra]\n[param]\n", "[param]\ntri=old\n tri = old\n",
              "[param]\ntridiag_method=lanczos\ntridiag_method=rkpw\n", "[param]\ntri=sc\n") {
    stage('rkpw', 'initializer', $text, 0);
}
stage('unknown', 'tool', "[param]\n", 0);
stage('rkpw', 'unknown', "[param]\n", 0);

my $canonical = stage('rkpw', 'tool', "  [param] # comment\n band=flat\nLambda=2\nNmax=1\nmMAX=2\n", 1);
like($canonical, qr/\A\[param\]\n/, 'canonical header understood by both producer parsers');
if (defined $nrgchain) {
    $nrgchain = File::Spec->rel2abs($nrgchain);
    my $pid = fork();
    die $! unless defined $pid;
    if (!$pid) {
        chdir "$work/run" or die $!;
        open STDOUT, '>', 'producer.out' or die $!;
        open STDERR, '>', 'producer.err' or die $!;
        exec {$nrgchain} $nrgchain, '-v', 's', 'param';
        die $!;
    }
    waitpid($pid, 0);
    is($?, 0, 'real tool consumes the staged noncanonical-header input');
    like(read_text("$work/run/producer.err"), qr/tridiag_method=rkpw/, 'producer reports explicit selected backend');
}

write_text("$work/source/param", "[param]\ntri=old\ntridiag_method=lanczos\n");
for (1 .. 2) {
    write_text("$work/run/param", "[param]\ntri=broken\ntridiag_method=broken\n");
    my $status = system($^X, $helper, '--backend', 'rkpw', '--mode', 'initializer',
        '--source', "$work/source", '--work', "$work/run", '--param', 'param', '--',
        $^X, '-e', 'open my $f,"<","param" or die $!; local $/; my $p=<$f>; die unless $p =~ /^tri=rkpw$/m && $ENV{CHAIN_TEST_BACKEND} eq "rkpw";');
    is($status, 0, 'wrapper restages pristine selected inputs every invocation');
}
my $status = system($^X, $helper, '--backend', 'rkpw', '--mode', 'tool', '--source', "$work/source",
                    '--work', "$work/run", '--', $^X, '-e', 'exit 7');
is($status >> 8, 7, 'producer failure is propagated');
done_testing();

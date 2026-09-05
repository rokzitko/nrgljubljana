#!/usr/bin/env perl

use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir tempfile);
use POSIX ();
use Test::More;

my ($test_dir) = @ARGV;
defined($test_dir) or die "Usage: comparators.t TEST_DIR\n";
$test_dir = File::Spec->rel2abs($test_dir);

my $mycomp = File::Spec->catfile($test_dir, 'mycomp.pl');
my $mycomp_signs = File::Spec->catfile($test_dir, 'mycomp-ignoresigns.pl');
my $compare = File::Spec->catfile($test_dir, 'compare.pl');
my $compare_signs = File::Spec->catfile($test_dir, 'compare-ignoresigns.pl');
my $clean = File::Spec->catfile($test_dir, 'clean-physical-output.pl');

sub write_file {
    my ($path, $content) = @_;
    open(my $handle, '>', $path) or die "Can't write $path: $!\n";
    print {$handle} $content or die "Can't write $path: $!\n";
    close($handle) or die "Can't finish $path: $!\n";
}

sub read_file {
    my ($path) = @_;
    open(my $handle, '<', $path) or die "Can't read $path: $!\n";
    my $content = do { local $/; <$handle> };
    close($handle) or die "Can't finish $path: $!\n";
    return $content;
}

sub run_command {
    my ($directory, @command) = @_;
    my ($output, $output_path) = tempfile(UNLINK => 1);
    my ($error, $error_path) = tempfile(UNLINK => 1);
    my $pid = fork();
    defined($pid) or die "Can't fork: $!\n";
    if ($pid == 0) {
        chdir($directory) or POSIX::_exit(126);
        open(STDOUT, '>&', $output) or POSIX::_exit(126);
        open(STDERR, '>&', $error) or POSIX::_exit(126);
        exec @command or POSIX::_exit(127);
    }
    waitpid($pid, 0);
    my $child_status = $?;
    seek($output, 0, 0) or die "Can't rewind $output_path: $!\n";
    seek($error, 0, 0) or die "Can't rewind $error_path: $!\n";
    my $stdout = do { local $/; <$output> } // '';
    my $stderr = do { local $/; <$error> } // '';
    close($output) or die "Can't close $output_path: $!\n";
    close($error) or die "Can't close $error_path: $!\n";
    my $signal = $child_status & 127;
    my $status = $signal ? 128 + $signal : $child_status >> 8;
    return ($status, $stdout, $stderr);
}

sub compare_text {
    my ($directory, $expected, $actual, @options) = @_;
    write_file("$directory/expected", $expected);
    write_file("$directory/actual", $actual);
    return run_command($directory, $^X, $mycomp, @options, 'expected', 'actual');
}

subtest 'finite numeric and structural comparison' => sub {
    my $dir = tempdir(CLEANUP => 1);
    my ($status) = compare_text(
        $dir,
        "# T value\nI=0 E=1 <n_d>=0.5 z=(1,2)\n",
        "# T value\nI=0 E=1.000001 <n_d>=0.500001 z=(1.000001,1.999999)\n",
        '--check-comments');
    is($status, 0, 'close labeled real and complex values pass');

    ($status) = compare_text($dir, "I=0 E=1\n", "J=0 E=1\n");
    is($status, 1, 'label changes fail');
    ($status) = compare_text($dir, "# x y\n1 2\n", "# y x\n1 2\n", '--check-comments');
    is($status, 1, 'checked header changes fail');
    ($status) = compare_text($dir, "1 2\n", "1 2 3\n");
    is($status, 1, 'field-count changes fail');
    ($status) = compare_text($dir, "1\n2\n", "1\n");
    is($status, 1, 'record-count changes fail');
    ($status) = compare_text($dir, "4.60495394416e-11\n", "0\n");
    is($status, 1, 'a small physical tail cannot collapse to zero');
    ($status) = compare_text($dir, "1000\n", "1000.009\n");
    is($status, 0, 'relative tolerance scales for large finite values');
    ($status) = compare_text($dir, "0\n", "0.000001\n", '--atol', '0.000001', '--rtol', '0');
    is($status, 0, 'explicit absolute tolerance includes its boundary');
    ($status) = compare_text($dir, "0\n", "0.0000011\n", '--atol', '0.000001', '--rtol', '0');
    is($status, 1, 'explicit absolute tolerance rejects values outside its boundary');
};

subtest 'invalid values never compare equal' => sub {
    my $dir = tempdir(CLEANUP => 1);
    for my $invalid ('nan', '-NAN', 'inf', '-Infinity', 'nan(payload)', '1e9999', '1e-9999', '1e', '(1,nan)') {
        my ($status) = compare_text($dir, "$invalid\n", "$invalid\n");
        is($status, 2, "identical invalid value '$invalid' is rejected as input");
    }

    my ($status) = compare_text($dir, "1\n", "garbage\n");
    is($status, 1, 'numeric-to-text replacement fails');
    ($status) = compare_text($dir, "garbage\n", "1\n");
    is($status, 1, 'text-to-numeric replacement fails');
    ($status) = compare_text($dir, "foo\n", "bar\n");
    is($status, 1, 'unequal text fails');

    write_file("$dir/empty", '');
    write_file("$dir/data", "1\n");
    ($status) = run_command($dir, $^X, $mycomp, 'empty', 'data');
    is($status, 2, 'empty files are rejected');
    write_file("$dir/comments", "# only metadata\n\n");
    ($status) = run_command($dir, $^X, $mycomp, 'comments', 'comments');
    is($status, 2, 'comments-only files are rejected');
    ($status) = run_command($dir, $^X, $mycomp, 'missing', 'data');
    is($status, 2, 'missing files are rejected');

    my $many_expected = "foo\n" x 256;
    my $many_actual = "bar\n" x 256;
    ($status) = compare_text($dir, $many_expected, $many_actual);
    is($status, 1, '256 differences cannot wrap to a successful exit status');
};

subtest 'sign relaxation is explicit and field-aware' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/expected", "x=1 z=(1,2)\n");
    write_file("$dir/actual", "x=-1 z=(-1,-2)\n");
    my ($status) = run_command($dir, $^X, $mycomp_signs, 'expected', 'actual');
    is($status, 0, 'scalar and whole-complex sign flips pass in sign mode');
    write_file("$dir/actual", "x=-1 z=(1,-2)\n");
    ($status) = run_command($dir, $^X, $mycomp_signs, 'expected', 'actual');
    is($status, 1, 'independent complex-component sign changes fail');
    write_file("$dir/actual", "x=-1 z=(0,2.2360679775)\n");
    ($status) = run_command($dir, $^X, $mycomp_signs, 'expected', 'actual');
    is($status, 1, 'equal complex magnitude with a different phase fails');

    make_path("$dir/ref", "$dir/results");
    write_file("$dir/ref/data", "1\n");
    write_file("$dir/results/data", "-1\n");
    ($status) = run_command($dir, $^X, $compare_signs, '--actual', "$dir/results", "$dir/ref");
    is($status, 0, 'directory sign wrapper forwards comparator options');
};

subtest 'per-reference tolerance overrides' => sub {
    my $dir = tempdir(CLEANUP => 1);
    my $ref = "$dir/ref";
    my $actual = "$dir/actual";
    make_path($ref, $actual);
    write_file("$ref/result.dat", "0\n");
    write_file("$actual/result.dat", "5e-10\n");

    my ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 1, 'default tolerance rejects the near-zero difference');
    write_file("$ref/result.dat.tol", "  # Absolute tolerance\n\n1e-9\n");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 0, 'one tolerance overrides the absolute tolerance without an actual-side companion');

    write_file("$ref/result.dat", "1000\n");
    write_file("$actual/result.dat", "1000.009\n");
    write_file("$ref/result.dat.tol", "0\n");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 0, 'one tolerance retains the active relative tolerance');

    write_file("$actual/result.dat", "1015\n");
    write_file("$ref/result.dat.tol", "# Absolute\n0\n# Relative\n2e-2\n");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 0, 'two tolerances override absolute and relative tolerances');
    write_file("$ref/result.dat.tol", "0 1e-3\n");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 1, 'a tighter relative override rejects the same result');

    write_file("$ref/result.dat", "0\n");
    write_file("$actual/result.dat", "0\n");
    my @invalid = (
        ['empty', ''],
        ['comment-only', "# no tolerances\n"],
        ['too many values', "1 2 3\n"],
        ['negative value', "-1\n"],
        ['non-finite value', "nan\n"],
        ['overflowing value', "1e9999\n"],
        ['underflowing value', "1e-9999\n"],
        ['inline comment', "1e-9 # comments must occupy a full line\n"],
    );
    for my $case (@invalid) {
        write_file("$ref/result.dat.tol", $case->[1]);
        ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
        is($status, 2, "$case->[0] is a tolerance-file configuration error");
    }

    unlink("$ref/result.dat.tol") or die $!;
    write_file("$actual/result.dat", "5e-10\n");
    write_file("$actual/result.dat.tol", "1\n");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 1, 'an actual-side tolerance file is ignored');

    my $physical_dir = tempdir(CLEANUP => 1);
    my $physical_ref = "$physical_dir/ref";
    my $physical_actual = "$physical_dir/actual";
    make_path($physical_ref, $physical_actual);
    write_file("$physical_ref/td", "0\n");
    write_file("$physical_actual/td", "5e-10\n");
    ($status) = run_command($physical_dir, $^X, $compare, '--strict', '--actual', $physical_actual, $physical_ref);
    is($status, 0, 'strict physical comparison uses its normal absolute tolerance');
    write_file("$physical_ref/td.tol", "1e-12\n");
    ($status) = run_command($physical_dir, $^X, $compare, '--strict', '--actual', $physical_actual, $physical_ref);
    is($status, 1, 'a physical-output tolerance file can tighten the absolute tolerance');

    my $spectral_dir = tempdir(CLEANUP => 1);
    my $spectral_ref = "$spectral_dir/ref";
    my $spectral_actual = "$spectral_dir/actual";
    make_path($spectral_ref, $spectral_actual);
    my $spectral_name = 'spec_FDM_dens_x-x.dat';
    write_file("$spectral_ref/$spectral_name", "0 1\n");
    write_file("$spectral_actual/$spectral_name", "0 1.015\n");
    write_file("$spectral_ref/$spectral_name.tol", "0\n");
    ($status) = run_command($spectral_dir, $^X, $compare, '--strict', '--actual', $spectral_actual, $spectral_ref);
    is($status, 0, 'one spectral override retains the two-percent value tolerance');
    write_file("$spectral_ref/$spectral_name.tol", "0 1e-1\n");
    write_file("$spectral_actual/$spectral_name", "0 1.05\n");
    ($status) = run_command($spectral_dir, $^X, $compare, '--strict', '--actual', $spectral_actual, $spectral_ref);
    is($status, 0, 'a spectral override can loosen dependent-value comparison');
    write_file("$spectral_ref/$spectral_name.tol", "1e-3 1e-1\n");
    write_file("$spectral_actual/$spectral_name", "1e-8 1.05\n");
    ($status) = run_command($spectral_dir, $^X, $compare, '--strict', '--actual', $spectral_actual, $spectral_ref);
    is($status, 1, 'a spectral absolute override does not loosen the first-field grid comparison');
    write_file("$spectral_ref/$spectral_name", "1 1\n");
    write_file("$spectral_actual/$spectral_name", "1.5 1\n");
    write_file("$spectral_ref/$spectral_name.tol", "0 1\n");
    ($status) = run_command($spectral_dir, $^X, $compare, '--strict', '--actual', $spectral_actual, $spectral_ref);
    is($status, 1, 'a spectral relative override does not loosen the first-field grid comparison');

    my $association_dir = tempdir(CLEANUP => 1);
    my $association_ref = "$association_dir/ref";
    my $association_actual = "$association_dir/actual";
    make_path($association_ref, $association_actual);
    write_file("$association_ref/keep", "1\n");
    write_file("$association_actual/keep", "1\n");
    write_file("$association_ref/missing.tol", "1e-9\n");
    ($status) = run_command($association_dir, $^X, $compare, '--actual', $association_actual, $association_ref);
    is($status, 2, 'an orphan tolerance file is a configuration error');
    unlink("$association_ref/missing.tol") or die $!;

    write_file("$association_ref/.tol", "1e-9\n");
    ($status) = run_command($association_dir, $^X, $compare, '--actual', $association_actual, $association_ref);
    is($status, 2, 'an empty tolerance target is a configuration error');
    unlink("$association_ref/.tol") or die $!;
    write_file("$association_ref/keep.tol.tol", "1e-9\n");
    ($status) = run_command($association_dir, $^X, $compare, '--actual', $association_actual, $association_ref);
    is($status, 2, 'a recursive tolerance filename is a configuration error');
    unlink("$association_ref/keep.tol.tol") or die $!;
    make_path("$association_ref/keep.tol");
    ($status) = run_command($association_dir, $^X, $compare, '--actual', $association_actual, $association_ref);
    is($status, 2, 'a nonregular tolerance entry is a configuration error');
    rmdir("$association_ref/keep.tol") or die $!;

    write_file("$association_ref/payload.bin", "payload\n");
    write_file("$association_actual/payload.bin", "payload\n");
    write_file("$association_ref/payload.bin.tol", "1e-9\n");
    ($status) = run_command($association_dir, $^X, $compare, '--actual', $association_actual, $association_ref);
    is($status, 2, 'a tolerance file cannot target an ordinary binary reference');
    unlink("$association_ref/payload.bin", "$association_ref/payload.bin.tol", "$association_actual/payload.bin") == 3
        or die $!;

    write_file("$association_ref/result", "0\n");
    write_file("$association_actual/result", "5e-10\n");
    write_file("$association_ref/result.tol", "invalid\n");
    ($status) = run_command($association_dir, $^X, $compare, '--exclude', 'result',
                            '--actual', $association_actual, $association_ref);
    is($status, 0, 'excluding a reference also excludes its tolerance file');
    ($status) = run_command($association_dir, $^X, $compare, '--exclude', 'result.tol',
                            '--actual', $association_actual, $association_ref);
    is($status, 1, 'excluding only a malformed tolerance file uses the normal comparison policy');

    my $semantic_dir = tempdir(CLEANUP => 1);
    my $semantic_ref = "$semantic_dir/ref";
    my $semantic_actual = "$semantic_dir/actual";
    make_path($semantic_ref, $semantic_actual);
    write_file("$semantic_ref/.physical-outputs", "subspaces subspaces.dat\n");
    write_file("$semantic_ref/subspaces.dat.tol", "1e-9\n");
    write_file("$semantic_actual/subspaces.dat", "Iteration 0\nlen_dm=1\nI=0 kept=1 total=1\n");
    ($status) = run_command($semantic_dir, $^X, $compare, '--actual', $semantic_actual, $semantic_ref);
    is($status, 2, 'a tolerance file cannot target a semantically validated output');
};

subtest 'directory manifests and physical output cleanup' => sub {
    my $dir = tempdir(CLEANUP => 1);
    my $ref = "$dir/ref";
    my $actual = "$dir/actual";
    make_path($ref, $actual);
    write_file("$ref/result [1].dat", "1\n");
    write_file("$actual/result [1].dat", "1\n");
    write_file("$ref/.hidden", "2\n");
    write_file("$actual/.hidden", "2\n");
    my ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 0, 'special and hidden reference names are compared safely');
    unlink("$actual/.hidden") or die $!;
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $ref);
    is($status, 1, 'missing hidden results fail');

    my $empty_ref = "$dir/empty-ref";
    make_path($empty_ref);
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $empty_ref);
    is($status, 2, 'empty reference directories fail');
    my $nested_ref = "$dir/nested-ref";
    make_path("$nested_ref/nested");
    ($status) = run_command($dir, $^X, $compare, '--actual', $actual, $nested_ref);
    is($status, 2, 'nested reference entries fail explicitly');

    my $strict_ref = "$dir/strict-ref";
    my $strict_actual = "$dir/strict-actual";
    make_path($strict_ref, $strict_actual);
    write_file("$strict_ref/data", "1\n");
    write_file("$strict_actual/data", "1\n");
    write_file("$strict_actual/td", "# T value\n1 2\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'unreferenced physical text output fails strict comparison');
    unlink("$strict_actual/td") or die $!;

    write_file("$strict_ref/spec_FDM_dens_x-x.dat", "0 0.000001\n");
    write_file("$strict_actual/spec_FDM_dens_x-x.dat", "0 0.000001015\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'strict spectral comparison permits two-percent relative drift');
    write_file("$strict_actual/spec_FDM_dens_x-x.dat", "0 0.00000103\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'strict spectral comparison rejects drift above two percent');
    write_file("$strict_actual/spec_FDM_dens_x-x.dat", "0.000001 0.000001\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'strict spectral comparison keeps the frequency grid tight');
    unlink("$strict_ref/spec_FDM_dens_x-x.dat") or die $!;
    unlink("$strict_actual/spec_FDM_dens_x-x.dat") or die $!;

    write_file("$strict_actual/spec_FDM_dens_x-x.bin", "binary\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'unreferenced physical binary output fails strict comparison');
    write_file("$strict_ref/.physical-outputs", "binary-real spec_FDM_dens_x-x.bin\n");
    write_file("$strict_actual/spec_FDM_dens_x-x.bin", pack('d*', 1, 2, -1, 3));
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'a declared finite binary spectrum passes semantic validation');
    write_file("$strict_actual/spec_FDM_dens_x-x.bin", "partial");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a declared binary spectrum with a partial record fails');

    write_file("$strict_ref/.physical-outputs", "subspaces subspaces.dat\n");
    unlink("$strict_actual/spec_FDM_dens_x-x.bin") or die $!;
    write_file("$strict_actual/subspaces.dat", "Iteration 0\nlen_dm=1\nI=0 kept=1 total=1\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'a structurally valid subspace dump passes semantic validation');
    write_file("$strict_actual/subspaces.dat", "garbage\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a malformed subspace dump fails semantic validation');
    write_file("$strict_actual/subspaces.dat", "Iteration 0\nlen_dm=1\nI=0 kept=1 total=1\nIteration 0\nlen_dm=1\nI=0 kept=1 total=1\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a contradictory subspace dump fails semantic validation');
    unlink("$strict_actual/subspaces.dat") or die $!;

    write_file("$strict_ref/.physical-outputs", "states states.nrg\n");
    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0 1
Vectors:
vec(0)=[(1,0), (0,0)] norm-1=0
vec(1)=[(0,0), (1,0)] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'a complete finite state-vector dump passes semantic validation');
    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0 1
Vectors:
vec(0)=[1, 0] norm-1=0
vec(1)=[1, 0] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'non-orthogonal state vectors fail semantic validation');
    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[] norm-1=-1
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'an empty state vector fails semantic validation');
    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1e-9999] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'an underflowed state-vector coefficient fails semantic validation');

    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
STATES
    my $energy_dump = <<'ENERGIES';
===== Iteration number: 0
Subspace: 0 1
0
===== Iteration number: 1
Subspace: 0 1
0
ENERGIES
    write_file("$strict_ref/energies.nrg", $energy_dump);
    write_file("$strict_actual/energies.nrg", $energy_dump);
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a truncated state-vector dump cannot omit energy-dump iterations');

    $energy_dump = <<'ENERGIES';
===== Iteration number: 0
Subspace: 0 1
0
Subspace: 1 1
1
ENERGIES
    write_file("$strict_ref/energies.nrg", $energy_dump);
    write_file("$strict_actual/energies.nrg", $energy_dump);
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a state-vector dump cannot omit an energy-dump subspace');

    $energy_dump = <<'ENERGIES';
===== Iteration number: 0
Subspace: 0 1
0 1
ENERGIES
    write_file("$strict_ref/energies.nrg", $energy_dump);
    write_file("$strict_actual/energies.nrg", $energy_dump);
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'a state-vector dump cannot omit an energy-dump eigenpair');

    $energy_dump = <<'ENERGIES';
===== Iteration number: 0
Subspace: 0 1
1
ENERGIES
    write_file("$strict_ref/energies.nrg", $energy_dump);
    write_file("$strict_actual/energies.nrg", $energy_dump);
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'state-vector energies must match the energy dump');

    $energy_dump = <<'ENERGIES';
===== Iteration number: 0
Subspace: 0 1
0
Subspace: 1 1
0
===== Iteration number: 0
Subspace: 2 1
0
ENERGIES
    write_file("$strict_ref/energies.nrg", $energy_dump);
    write_file("$strict_actual/energies.nrg", $energy_dump);
    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
===== Iteration number: 0
Subspace: 1 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
Subspace: 2 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'state-vector subspaces cannot move across duplicate iteration-zero sections');
    unlink("$strict_ref/energies.nrg") or die $!;
    unlink("$strict_actual/energies.nrg") or die $!;

    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
Subspace: +0 01
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'equivalent integer spellings cannot duplicate a state-vector subspace');

    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 5
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
===== Iteration number: 5
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
STATES
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 1, 'only the initial zero iteration may be duplicated');

    write_file("$strict_actual/states.nrg", <<'STATES');
===== Iteration number: 0
Subspace: 0 1
Energies (rel): 0
Vectors:
vec(0)=[1] norm-1=0
STATES
    write_file("$strict_actual/energies.nrg", "stale\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--exclude', 'energies.nrg',
                            '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'an excluded energy dump does not constrain state-vector validation');
    unlink("$strict_actual/energies.nrg") or die $!;
    unlink("$strict_actual/states.nrg") or die $!;

    write_file("$strict_ref/.physical-outputs", "binary-real spec_x/../../outside.dat\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 2, 'manifest output names must be basenames');

    write_file("$strict_ref/.physical-outputs", "hdf5 raw.h5\n");
    write_file("$strict_actual/raw.h5", "not really hdf5\n");
    {
        local $ENV{PATH} = "$dir/no-programs";
        make_path($ENV{PATH});
        ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    }
    is($status, 2, 'missing semantic-validator tooling is an input error');
    unlink("$strict_actual/raw.h5") or die $!;
    unlink("$strict_ref/.physical-outputs") or die $!;
    write_file("$strict_actual/param", "control\n");
    write_file("$strict_actual/log", "diagnostic\n");
    write_file("$strict_actual/DONE", "\n");
    write_file("$strict_actual/input.dat", "input\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'known controls and non-output inputs are outside the strict manifest');

    write_file("$strict_ref/td", "# T value\n1 2\n");
    write_file("$strict_actual/td", "# T value\n1 2\n");
    ($status) = run_command($dir, $^X, $compare, '--strict', '--exclude', 'data', '--actual', $strict_actual, $strict_ref);
    is($status, 0, 'an explicitly excluded generated input can use another comparator policy');

    write_file("$strict_actual/custom", "1\n");
    write_file("$strict_actual/raw.h5", "not really hdf5\n");
    ($status) = run_command($strict_actual, $^X, $clean);
    is($status, 0, 'physical output cleanup succeeds');
    ok(!-e "$strict_actual/td" && !-e "$strict_actual/custom" && !-e "$strict_actual/raw.h5" && !-e "$strict_actual/DONE",
       'physical outputs and completion marker are removed');
    is(read_file("$strict_actual/data"), "1\n", 'generated input is preserved by cleanup');
    is(read_file("$strict_actual/param"), "control\n", 'parameter file is preserved by cleanup');

    make_path("$strict_actual/td");
    ($status) = run_command($strict_actual, $^X, $clean);
    isnt($status, 0, 'cleanup refuses an output-like directory');
    rmdir("$strict_actual/td") or die $!;
};

done_testing();

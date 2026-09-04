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

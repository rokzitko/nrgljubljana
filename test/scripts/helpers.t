#!/usr/bin/env perl

use strict;
use warnings;

use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir tempfile);
use POSIX ();
use Test::More;

my ($scripts, $tools) = @ARGV;
defined($scripts) && defined($tools) or die "Usage: helpers.t SCRIPTS TOOLS\n";
$scripts = File::Spec->rel2abs($scripts);
$tools = File::Spec->rel2abs($tools);

sub write_file {
    my ($path, $content) = @_;
    open(my $fh, '>', $path) or die "Can't write $path: $!\n";
    print {$fh} $content or die "Can't write $path: $!\n";
    close($fh) or die "Can't finish $path: $!\n";
}

sub read_file {
    my ($path) = @_;
    open(my $fh, '<', $path) or die "Can't read $path: $!\n";
    my $content = do { local $/; <$fh> };
    close($fh) or die "Can't finish $path: $!\n";
    return $content;
}

sub run_command {
    my ($directory, @command) = @_;
    my ($output, $output_path) = tempfile();
    my ($error, $error_path) = tempfile();
    my $pid = fork();
    defined($pid) or die "Can't fork: $!\n";
    if ($pid == 0) {
        chdir($directory) or POSIX::_exit(126);
        open(STDIN, '<', File::Spec->devnull()) or POSIX::_exit(126);
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
    my $signal = $child_status & 127;
    my $status = $signal ? 128 + $signal : $child_status >> 8;
    return ($status, $stdout, $stderr);
}

sub script {
    return File::Spec->catfile($scripts, shift);
}

subtest 'parameter and iteration state' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/loop", "\$Nz=4;\n");
    my ($status, $out) = run_command($dir, $^X, script('getNz'), 'loop');
    is($status, 0, 'getNz accepts an explicit filename');
    is($out, "4\n", 'getNz extracts one value');

    write_file("$dir/param", "a.b=7\naxb=8\n");
    ($status, $out) = run_command($dir, $^X, script('getparam'), 'a.b');
    is($status, 0, 'literal parameter name succeeds');
    is($out, "7\n", 'regex characters are literal');

    ($status, $out) = run_command($dir, $^X, script('getiter'));
    is_deeply([$status, $out], [0, '0'], 'getiter initializes zero without a newline');
    ($status) = run_command($dir, $^X, script('newiter'));
    is($status, 0, 'newiter increments under lock');
    ($status, $out) = run_command($dir, $^X, script('getiter'));
    is_deeply([$status, $out], [0, '1'], 'increment is persisted');
    chmod(0444, "$dir/ITER") or die "Can't make ITER read-only: $!\n";
    ($status, $out) = run_command($dir, $^X, script('getiter'));
    is_deeply([$status, $out], [0, '1'], 'getiter accepts read-only state');
};

subtest 'destructive helpers stay contained' => sub {
    my $dir = tempdir(CLEANUP => 1);
    for my $name (1, 3) {
        make_path("$dir/$name");
        write_file("$dir/$name/param", "p\n");
        write_file("$dir/$name/data", "d\n");
    }
    make_path("$dir/2");
    write_file("$dir/2/param", "keep\n");
    my ($status) = run_command($dir, $^X, script('delete_NRG_output'));
    is($status, 0, 'cleanup succeeds across numbering gaps');
    ok(!-e "$dir/1" && !-e "$dir/3", 'marked directories removed');
    ok(-d "$dir/2", 'unmarked directory retained');

    make_path("$dir/job/dmft", "$dir/job/res");
    write_file("$dir/job/param.loop", "loop\n");
    write_file("$dir/job/dmft/value", "1\n");
    write_file("$dir/job/res/value", "2\n");
    ($status) = run_command($dir, $^X, script('compress_leafs_DMFT'));
    is($status, 0, 'DMFT archive succeeds');
    ok(!-e "$dir/job" && -s "$dir/job.tgz", 'source replaced by verified archive');
};

subtest 'normalization and table pairing' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/table", "# data\n0 2\n1 2\n");
    local $ENV{PATH} = "$tools:$scripts:$ENV{PATH}";
    my ($status, $out) = run_command($dir, $^X, script('normalize'), 'table');
    is($status, 0, 'normalize succeeds');
    is($out, "# data\n0 1\n1 1\n", 'normalize preserves metadata and scales y');

    write_file("$dir/a", "# first\n0 4\n\n1 6\n");
    write_file("$dir/b", "0 2\n# second\n1 3\n");
    ($status, $out) = run_command($dir, $^X, script('divy'), 'a', 'b');
    is($status, 0, 'pairwise comments are independent');
    is($out, "# first\n0 2\n\n1 2\n", 'first-file layout is preserved');

    write_file("$dir/b", "0 2\n2 3\n");
    ($status, $out) = run_command($dir, $^X, script('divy'), 'a', 'b');
    isnt($status, 0, 'coordinate mismatch fails');
    is($out, '', 'pairwise failure emits no partial table');
};

subtest 'averages and exact moments' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/a", "0 1\n1 3\n");
    write_file("$dir/b", "0 3\n1 5\n");
    my ($status, $out) = run_command($dir, $^X, script('averagefiles'), 'a', 'b');
    is($status, 0, 'averagefiles succeeds');
    is($out, "0 2\n1 4\n", 'only dependent values are averaged');

    ($status, $out) = run_command($dir, $^X, script('averagefiles_mc'), 'a', 'b');
    is($status, 0, 'MC average succeeds with two samples');
    is($out, "0 2 1\n1 4 1\n", 'sample SEM uses Welford M2');

    write_file("$dir/a", "0 1e308\n");
    write_file("$dir/b", "0 -1e308\n");
    ($status, $out) = run_command($dir, $^X, script('averagefiles'), 'a', 'b');
    is_deeply([$status, $out], [0, "0 0\n"], 'scaled averaging avoids overflow');

    write_file("$dir/linear", "0 0\n1 1\n");
    ($status, $out) = run_command($dir, $^X, script('moments'), 'linear');
    is($status, 0, 'moments succeeds');
    my @found = split /\n/, $out;
    my @expected = (1 / 3, 1 / 4, 1 / 5, 1 / 6);
    ok(abs($found[$_] - $expected[$_]) < 1e-15, "moment " . ($_ + 1)) for 0 .. 3;
    write_file("$dir/linear", "0 1e308\n1 -1e308\n");
    ($status, $out) = run_command($dir, $^X, script('moments'), 'linear');
    is($status, 0, 'finite extreme moments do not overflow endpoint differences');
    unlike($out, qr/(?:inf|nan)/i, 'extreme moment output remains finite');
};

subtest 'collectors publish only valid complete results' => sub {
    my $dir = tempdir(CLEANUP => 1);
    for my $run (1, 2) {
        make_path("$dir/$run");
        write_file("$dir/$run/DONE", '');
        write_file("$dir/$run/custom", "# x y\n0 $run\n1 " . ($run + 1) . "\n");
    }
    my ($status) = run_command($dir, $^X, script('gatherlastlines'));
    is($status, 0, 'last-line collector succeeds');
    is(read_file("$dir/custom"), "# x y\n1 2\n1 3\n", 'one header and last logical rows emitted');
    for my $run (1, 2) {
        write_file("$dir/$run/param", "P=" . (10 * $run) . "\n");
    }
    ($status) = run_command($dir, $^X, script('gather'), 'P');
    is($status, 0, 'parameter collector succeeds');
    like(read_file("$dir/custom.dat"), qr/\A# P y\n/, 'gather relabels the replaced coordinate');

    write_file("$dir/td.dat", "sentinel\n");
    write_file("$dir/1/td", "# x y\n1 1\n");
    write_file("$dir/2/td", "# x y\nmalformed\n");
    ($status) = run_command($dir, $^X, script('gathertd'));
    isnt($status, 0, 'malformed collector input fails');
    is(read_file("$dir/td.dat"), "sentinel\n", 'failed collector preserves old result');

    make_path("$dir/A=2", "$dir/A=1");
    write_file("$dir/A=2/result", "# x y\n0 20\n1 21\n");
    write_file("$dir/A=1/result", "# x y\n0 10\n1 11\n");
    ($status) = run_command($dir, $^X, script('gatherfiles'), 'result');
    is($status, 0, 'gatherfiles succeeds');
    is(read_file("$dir/result"), "# A x y\n1 0 10\n1 1 11\n2 0 20\n2 1 21\n",
       'every source row is prefixed and sorted');
};

subtest 'symmetry, mesh, and reports' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/symmetric", "# x y\n-1 2\n0 7\n1 4\n");
    my ($status, $out) = run_command($dir, $^X, script('symmetrize'), 'symmetric');
    is($status, 0, 'odd symmetric grid succeeds');
    is(read_file("$dir/symmetric"), "# x y\n-1 3\n0 7\n1 3\n", 'center retained and pairs averaged');

    write_file("$dir/zero", "0 0\n1 0\n2 0\n");
    ($status, $out) = run_command($dir, $^X, script('optimize_mesh'), 'zero', '0.1');
    is($status, 0, 'zero-valued mesh does not divide by zero');
    is($out, "0 0\n1 0\n2 0\n", 'zero-valued mesh is unchanged');
    write_file("$dir/large", "0 1e308\n1 1e308\n2 1e308\n");
    ($status) = run_command($dir, $^X, script('optimize_mesh'), 'large', '0.1');
    is($status, 0, 'large constant mesh avoids scale overflow');

    write_file("$dir/m1", "1 0\n1.000000005 0\n");
    write_file("$dir/m2", "0 0\n");
    ($status, $out) = run_command($dir, $^X, script('mesh_union'), 'm1', 'm2');
    is($status, 0, 'mesh union succeeds');
    is($out, "-1 0\n0 0\n1 0\n", 'deduplicated mesh is exactly symmetric');

    write_file("$dir/custom.dat", "# 1 2 3\n# T A B\n0 10 20\n1 11 21\n");
    ($status) = run_command($dir, $^X, script('report'), 'custom.dat', "$dir/out");
    is($status, 0, 'report succeeds');
    ok(!-e "$dir/out-T.dat", 'coordinate output is skipped');
    is(read_file("$dir/out-A.dat"), "# T A\n0 10\n1 11\n", 'value output is complete');
    is(read_file("$dir/out-B.dat"), "# T B\n0 20\n1 21\n", 'last value output is complete');
    ($status) = run_command($dir, $^X, script('report'), 'custom.dat', "$dir/out", 'T');
    isnt($status, 0, 'explicit coordinate selection fails');
};

done_testing();

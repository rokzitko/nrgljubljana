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
    my ($input, $input_path);
    if (@command && ref($command[0]) eq 'SCALAR') {
        my $content = ${shift @command};
        ($input, $input_path) = tempfile();
        print {$input} $content or die "Can't write $input_path: $!\n";
        seek($input, 0, 0) or die "Can't rewind $input_path: $!\n";
    }
    my ($output, $output_path) = tempfile();
    my ($error, $error_path) = tempfile();
    my $pid = fork();
    defined($pid) or die "Can't fork: $!\n";
    if ($pid == 0) {
        chdir($directory) or POSIX::_exit(126);
        open(STDIN, '<', $input_path // File::Spec->devnull()) or POSIX::_exit(126);
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

subtest 'setparam literal and atomic updates' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/params", "a.b=1\naxb=2\n");
    chmod(0640, "$dir/params") or die "Can't set params mode: $!\n";
    my ($status, $out, $err) = run_command(
        $dir, $^X, script('setparam'), 'a.b', '9', 'params');
    is_deeply([$status, $out, $err],
              [0, "setparam: keyword=a.b value=9 fn=params\n"
                  . "setparam done. nrsubst=1\n", ''],
              'literal replacement has exact successful diagnostics');
    is(read_file("$dir/params"), "a.b=9\naxb=2\n", 'axb is not replaced');
    my @stat = stat "$dir/params";
    is($stat[2] & 07777, 0640, 'replacement preserves file mode');

    write_file("$dir/param file", 'keep=1');
    ($status, $out, $err) = run_command(
        $dir, $^X, script('setparam'), 'new.key', '3', 'param file');
    is_deeply([$status, $out, $err],
              [0, "setparam: keyword=new.key value=3 fn=param file\n"
                  . "setparam done. nrsubst=0\n", ''],
              'unterminated file with a space is updated successfully');
    is(read_file("$dir/param file"), "keep=1\nnew.key=3\n",
       'append supplies the missing newline');

    my $duplicate = "a.b=1\na.b=2\n";
    write_file("$dir/duplicate", $duplicate);
    ($status, $out) = run_command(
        $dir, $^X, script('setparam'), 'a.b', '9', 'duplicate');
    isnt($status, 0, 'duplicate key is rejected');
    is($out, '', 'duplicate rejection emits no success diagnostics');
    is(read_file("$dir/duplicate"), $duplicate, 'duplicate rejection preserves content');

    ($status, $out, $err) = run_command(
        $dir, $^X, script('setparam'), "a.b\naxb", '9', 'params');
    isnt($status, 0, 'newline in a literal key is rejected');
    is_deeply([$out, $err], ['', "setparam: keyword contains CR or LF\n"],
              'newline failure occurs before output');
};

subtest 'custparam logical tail and literal parameter' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/custom.dat",
               "# stale\n0 1\n1 2\n\n# final comment\n\n");
    write_file("$dir/param.loop", "a.b=7\naxb=8\n");
    my ($status, $out, $err) = run_command(
        $dir, $^X, script('custparam'), 'a.b');
    is_deeply([$status, $out, $err],
              [0, "# final comment\n1 2   7\n", ''],
              'final logical row, comment, and literal value are selected');

    write_file("$dir/param.loop", "a.b=7\na.b=8\n");
    ($status, $out) = run_command($dir, $^X, script('custparam'), 'a.b');
    isnt($status, 0, 'duplicate parameter fails');
    is($out, '', 'duplicate parameter failure emits no stdout');
};

subtest 'literal numeric transforms and domains' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/table", "# values\n2 3 keep\n\n-1 4 tail\n");
    my ($status, $out, $err) = run_command(
        $dir, $^X, script('scalex'), '0.5', 'table');
    is_deeply([$status, $out, $err],
              [0, "# values\n1 3 keep\n\n-0.5 4 tail\n", ''],
              'scalex preserves metadata and trailing tokens');
    ($status, $out) = run_command(
        $dir, $^X, script('scalex'), '1/2', 'table');
    isnt($status, 0, 'numeric expression is rejected');
    is($out, '', 'invalid factor emits no table');

    write_file("$dir/unterminated", '# first source');
    write_file("$dir/following", "2 3\n");
    ($status, $out) = run_command(
        $dir, $^X, script('scalex'), '1', 'unterminated', 'following');
    is_deeply([$status, $out], [0, "# first source\n2 3\n"],
              'multiple sources cannot merge across a missing newline');

    write_file("$dir/band", "# band\n-2 a\n-1 b\n0 c\n1 d\n2 e\n");
    for my $case (
        ['cutoutlow',  "# band\n-2 a\n-1 b\n1 d\n2 e\n"],
        ['cutouthigh', "# band\n-1 b\n0 c\n1 d\n"],
    ) {
        ($status, $out) = run_command(
            $dir, $^X, script($case->[0]), '1', 'band');
        is_deeply([$status, $out], [0, $case->[1]],
                  "$case->[0] includes the cutoff boundary");
        ($status, $out) = run_command(
            $dir, $^X, script($case->[0]), '-1', 'band');
        isnt($status, 0, "$case->[0] rejects a negative cutoff");
        is($out, '', "$case->[0] rejection emits no table");
    }

    ($status, $out) = run_command($dir, $^X, script('tk'), '2', '0.5');
    is($status, 0, 'tk accepts an ordinary positive domain');
    like($out, qr/\Au=2\ngamma=0\.5\nrho_0 J=[^\n]+\nDeff=[^\n]+\nT_K=[^\n]+\n\z/,
         'tk reports all output labels');
    ($status, $out) = run_command(
        $dir, $^X, script('tk'), '5e-324', '1e-308');
    is($status, 0, 'tk avoids an underflowing effective-bandwidth intermediate');
    my ($tiny_tk) = $out =~ /T_K=([^\n]+)/;
    cmp_ok($tiny_tk // 0, '>', 0, 'tk retains a representable tiny temperature');
    my $input_u = 0 + '5e-324';
    my $input_gamma = 0 + '1e-308';
    my $expected_rho = (8 / 3.14159) * ($input_gamma / $input_u);
    my $expected_tk = exp(
        log(0.182) + log($input_u) + 0.5 * log($expected_rho) - 1 / $expected_rho);
    cmp_ok(abs($tiny_tk - $expected_tk) / $expected_tk, '<', 1e-12,
           'tk consistently uses the parsed input values');
    ($status, $out) = run_command(
        $dir, $^X, script('tk'), '2', '5e-324');
    is($status, 0, 'tk avoids an underflowing exchange-ratio intermediate');
    my ($tiny_rho) = $out =~ /rho_0 J=([^\n]+)/;
    cmp_ok($tiny_rho // 0, '>', 0, 'tk retains a representable tiny exchange ratio');
    for my $case (['U', '0', '1'], ['Gamma', '1', '0']) {
        ($status, $out) = run_command(
            $dir, $^X, script('tk'), $case->[1], $case->[2]);
        isnt($status, 0, "tk rejects nonpositive $case->[0]");
        is($out, '', "invalid $case->[0] emits no output");
    }
};

subtest 'weighted table mixing' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/a", "# mixed\n0 2\n\n# middle\n1 4\n");
    write_file("$dir/b", "0 10\n# ignored layout\n1 20\n");
    my ($status, $out) = run_command(
        $dir, $^X, script('mixfiles'), '2', 'a', '-0.5', 'b');
    is_deeply([$status, $out],
              [0, "# mixed\n0 -1\n\n# middle\n1 -2\n"],
              'arbitrary weights preserve first-file layout');

    write_file("$dir/b", "0 10\n2 20\n");
    ($status, $out) = run_command(
        $dir, $^X, script('mixfiles'), '1', 'a', '1', 'b');
    isnt($status, 0, 'grid mismatch fails');
    is($out, '', 'grid mismatch emits no partial table');

    write_file("$dir/huge-a", "0 1e308\n");
    write_file("$dir/huge-b", "0 1e308\n");
    ($status, $out) = run_command(
        $dir, $^X, script('mixfiles'), '1e308', 'huge-a', '-1e308', 'huge-b');
    is_deeply([$status, $out], [0, "0 0\n"],
              'extreme weighted cancellation remains finite');

    write_file("$dir/huge-b", "0 1e-308\n");
    ($status, $out) = run_command(
        $dir, $^X, script('mixfiles'), '1e-308', 'huge-a', '1e308', 'huge-b');
    is_deeply([$status, $out], [0, "0 2\n"],
              'mixed-scale products are accumulated before conversion');
};

subtest 'stable mean and sign filters' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/values", "# values\n1\n\n# more\n3\n");
    my ($status, $out) = run_command($dir, $^X, script('mean'), 'values');
    is_deeply([$status, $out], [0, "2\n"], 'mean ignores comments and blanks');
    write_file("$dir/values", "1e308\n-1e308\n");
    ($status, $out) = run_command($dir, $^X, script('mean'), 'values');
    is_deeply([$status, $out], [0, "0\n"], 'extreme mean is stable');
    write_file("$dir/values", "1e308\n3\n-1e308\n");
    ($status, $out) = run_command($dir, $^X, script('mean'), 'values');
    is_deeply([$status, $out], [0, "1\n"],
              'extreme cancellation retains a representable residual');

    write_file("$dir/ratio", "# ratio\n2 4 keep\n0 1 stop\n");
    ($status, $out) = run_command($dir, $^X, script('divybyx'), 'ratio');
    isnt($status, 0, 'divybyx rejects a zero coordinate');
    is($out, '', 'division failure emits no partial table');

    write_file("$dir/signs", "# signs\n-1 neg\n\n0 zero\n# middle\n1 pos\n");
    for my $case (
        ['positivepart', "# signs\n\n# middle\n1 pos\n"],
        ['negativepart', "# signs\n-1 neg\n\n# middle\n"],
    ) {
        ($status, $out) = run_command($dir, $^X, script($case->[0]), 'signs');
        is_deeply([$status, $out], [0, $case->[1]],
                  "$case->[0] preserves metadata and excludes zero");
    }
};

subtest 'positional column extraction' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/table data", "# labels\n1 a X\n\n2 b Y\n");
    my ($status, $out) = run_command(
        $dir, $^X, script('columns'), 'table data', '3', '1', '3');
    is_deeply([$status, $out], [0, "X 1 X\nY 2 Y\n"],
              'columns supports reordered repeated extraction');

    write_file("$dir/ragged", "1 a\n2 b extra\n");
    for my $case (
        ['invalid index', 'table data', '0'],
        ['out-of-range index', 'table data', '4'],
        ['ragged input', 'ragged', '1'],
    ) {
        ($status, $out) = run_command(
            $dir, $^X, script('columns'), $case->[1], $case->[2]);
        isnt($status, 0, "$case->[0] fails");
        is($out, '', "$case->[0] emits no output");
    }

    ($status, $out) = run_command(
        $dir, $^X, script('column'), 'table data', '2');
    is_deeply([$status, $out], [0, "a\nb\n"], 'column delegates to columns');
    my $stdin = "# stdin\n5 five\n6 six\n";
    ($status, $out) = run_command($dir, \$stdin, $^X, script('column1'));
    is_deeply([$status, $out], [0, "5\n6\n"], 'column1 reads standard input');
    ($status, $out) = run_command(
        $dir, $^X, script('column1'), 'table data');
    is_deeply([$status, $out], [0, "1\n2\n"],
              'column1 accepts a filename containing spaces');
    ($status, $out) = run_command(
        $dir, $^X, script('column2'), 'table data');
    is_deeply([$status, $out], [0, "a\nb\n"], 'column2 extracts column two');

    write_file("$dir/one-column", "9\n");
    ($status, $out) = run_command(
        $dir, $^X, script('column2'), 'table data', 'one-column');
    isnt($status, 0, 'column2 fails if a later input is invalid');
    is($out, '', 'column2 buffers earlier inputs until all validate');
    ($status, $out) = run_command(
        $dir, \$stdin, $^X, script('column1'), '-', '-');
    isnt($status, 0, 'column1 rejects repeated standard input');
    is($out, '', 'repeated standard input emits no partial output');
};

subtest 'named column extraction' => sub {
    my $dir = tempdir(CLEANUP => 1);
    write_file("$dir/named",
               "# stale stale stale\n# x alpha beta\n0 10 20\n"
               . "# beta beta beta\n1 11 21\n");
    my ($status, $out) = run_command(
        $dir, $^X, script('extractcolumns'), 'named', 'beta', 'x', 'beta');
    is_deeply([$status, $out], [0, "20 0 20\n21 1 21\n"],
              'final leading header controls ordered repeated extraction');

    ($status, $out) = run_command(
        $dir, $^X, script('extractcolumns'), 'named', 'alpha');
    my ($delegate_status, $delegate_out) = run_command(
        $dir, $^X, script('extractcolumn'), 'named', 'alpha');
    is_deeply([$status, $out], [0, "10\n11\n"],
              'comments after data are ignored');
    is_deeply([$delegate_status, $delegate_out], [$status, $out],
              'extractcolumn is equivalent to extractcolumns');

    write_file("$dir/no-header", "0 1\n");
    write_file("$dir/duplicate-header", "# x x\n0 1\n");
    for my $case (
        ['missing header', 'no-header'],
        ['duplicate header', 'duplicate-header'],
    ) {
        ($status, $out) = run_command(
            $dir, $^X, script('extractcolumns'), $case->[1], 'x');
        isnt($status, 0, "$case->[0] is rejected");
        is($out, '', "$case->[0] rejection emits no output");
    }
};

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

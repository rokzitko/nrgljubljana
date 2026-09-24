#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Spec;
use JSON::PP qw(decode_json);

my ($root, $cmake, $ctest) = @ARGV;
die "Usage: chain-registration.t TEST_SOURCE CMAKE CTEST\n" unless defined $ctest;
$root = File::Spec->rel2abs($root);
my $scratch = tempdir(CLEANUP => 1);
my $source = "$scratch/source";
mkdir $source or die $!;
mkdir "$source/fixture" or die $!;
symlink $root, "$source/test" or die $!;
sub write_text {
    open my $out, '>', $_[0] or die $!;
    print {$out} $_[1];
    close $out or die $!;
}
sub read_text {
    open my $in, '<', $_[0] or die $!;
    local $/;
    return <$in>;
}
sub run {
    my ($log, @command) = @_;
    my $pid = fork();
    die $! unless defined $pid;
    if (!$pid) {
        open STDOUT, '>', "$scratch/$log" or die $!;
        open STDERR, '>&', \*STDOUT or die $!;
        exec {$command[0]} @command;
        die $!;
    }
    waitpid($pid, 0);
    return $?;
}
write_text("$source/fixture/param", "[param]\ntri=old\ntridiag_method=lanczos\n");
write_text("$source/fixture/input.txt", "immutable fixture\n");
write_text("$source/CMakeLists.txt", <<'CMAKE');
cmake_minimum_required(VERSION 3.25)
project(chain_registration NONE)
enable_testing()
include(${PROJECT_SOURCE_DIR}/test/ChainBackendTests.cmake)
foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
  add_chain_test(example ${backend} MODE initializer SOURCE ${PROJECT_SOURCE_DIR}/fixture
    INPUTS param input.txt COMMAND ${CMAKE_COMMAND} -E true)
endforeach()
CMAKE

for my $legacy (qw(ON OFF)) {
    for my $rkpw (qw(ON OFF)) {
        my $build = "$scratch/$legacy-$rkpw";
        is(run('configure.log', $cmake, '-S', $source, '-B', $build,
               "-DTEST_CHAIN_LEGACY=$legacy", "-DTEST_CHAIN_RKPW=$rkpw"), 0, "configure $legacy/$rkpw")
            or diag read_text("$scratch/configure.log");
        is(run('manifest.json', $ctest, '--test-dir', $build, '--show-only=json-v1'), 0, 'read registration manifest');
        my $manifest = decode_json(read_text("$scratch/manifest.json"));
        my @expected = (($legacy eq 'ON' ? ('example_legacy') : ()), ($rkpw eq 'ON' ? ('example_rkpw') : ()));
        is_deeply([sort map { $_->{name} } @{$manifest->{tests}}], \@expected, 'only enabled independent registrations exist');
        for my $test (@{$manifest->{tests}}) {
            my $name = $test->{name};
            my ($backend) = $name =~ /_(legacy|rkpw)$/;
            my %properties = map { $_->{name} => $_->{value} } @{$test->{properties}};
            is($properties{WORKING_DIRECTORY}, "$build/$name", 'isolated suffixed working directory');
            ok(grep($_ eq "chain-$backend", @{$properties{LABELS}}), 'backend selection label');
            ok(!exists($properties{DEPENDS}) && !exists($properties{FIXTURES_REQUIRED}), 'no dependency on another backend test');
            my $tri = $backend eq 'legacy' ? 'old' : 'rkpw';
            like(read_text("$build/$name/param"), qr/^tri=$tri$/m, 'explicit staged initializer selection');
            is(run('test.log', $ctest, '--test-dir', $build, '-R', "^$name\$", '--output-on-failure', '--no-tests=error'),
               0, 'test runs on its own') or diag read_text("$scratch/test.log");
        }
    }
}
is(read_text("$source/fixture/param"), "[param]\ntri=old\ntridiag_method=lanczos\n", 'shared fixture remains unchanged');
done_testing();

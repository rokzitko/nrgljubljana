#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Test::More;

my $source = abs_path(shift @ARGV // die "Usage: $0 SOURCE_DIR\n");
my $work = tempdir(CLEANUP => 1);
make_path("$work/source/nrginit", "$work/fixture/ref", "$work/run");
symlink "$source/test", "$work/source/test" or die "symlink test: $!";

sub write_file {
    my ($path, $text) = @_;
    open my $out, '>', $path or die "$path: $!";
    print {$out} $text;
    close $out or die "$path: $!";
}

# A private fake source tree keeps the installed/local initializer untouched.
write_file("$work/source/nrginit/nrginit", <<'SH');
#!/bin/sh
set -eu
touch producer-ran
case "$1" in
  no-output) ;;
  primary-only) cp "$MOCK_REF/$MOCK_PRIMARY" "$MOCK_PRIMARY" ;;
  complete) cp "$MOCK_REF/$MOCK_PRIMARY" "$MOCK_PRIMARY"; cp "$MOCK_REF/ham_1" ham_1 ;;
  *) exit 2 ;;
esac
SH
chmod 0755, "$work/source/nrginit/nrginit" or die "chmod: $!";
chdir "$work/run" or die "chdir: $!";

for my $primary (qw(data data.in)) {
    unlink "$work/fixture/ref/data", "$work/fixture/ref/data.in";
    write_file("$work/fixture/ref/$primary", "1.0\n");
    write_file("$work/fixture/ref/ham_1", "2.0\n");
    local $ENV{MOCK_REF} = "$work/fixture/ref";
    local $ENV{MOCK_PRIMARY} = $primary;
    for my $mode (qw(no-output primary-only complete)) {
        write_file('data', "1.0\n");
        write_file('data.in', "1.0\n");
        write_file('ham_1', "2.0\n");
        unlink 'producer-ran';
        my $pid = fork();
        die "fork: $!" unless defined $pid;
        if (!$pid) {
            open STDOUT, '>', 'runner.log' or die $!;
            open STDERR, '>&', \*STDOUT or die $!;
            exec 'bash', "$source/test/nrginit/runtest", "$work/source", "$work/fixture", $mode;
            die "exec: $!";
        }
        waitpid($pid, 0);
        my $status = $?;
        ok(-e 'producer-ran', "$primary $mode producer ran successfully");
        if ($mode eq 'complete') {
            is($status, 0, "$primary accepts fresh complete outputs");
        } else {
            isnt($status, 0, "$primary rejects $mode despite seeded stale outputs");
            ok(!-e 'ham_1', "$primary stale auxiliary output removed");
            ok(!-e $primary, "$primary stale primary removed") if $mode eq 'no-output';
        }
    }
}
chdir $source or die "restore cwd: $!";
done_testing();

#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use JSON::PP qw(decode_json);
use File::Spec;

my ($executable, $backend, $extended) = @ARGV;
die "Usage: qualification-report.t DRIVER legacy|rkpw [--extended]\n"
    unless defined($backend) && $backend =~ /\A(?:legacy|rkpw)\z/
        && (!defined($extended) || $extended eq '--extended') && @ARGV <= 3;
$executable = File::Spec->rel2abs($executable);
my $root = tempdir(CLEANUP => 1);
my $next = 0;

sub run {
    my (@args) = @_;
    my $directory = "$root/" . ++$next;
    mkdir $directory or die $!;
    my $pid = fork(); die $! unless defined $pid;
    if (!$pid) {
        chdir $directory or die $!;
        open STDOUT, '>', 'log' or die $!;
        open STDERR, '>&', \*STDOUT or die $!;
        exec {$executable} $executable, @args;
        die $!;
    }
    waitpid($pid, 0);
    return ($?, $directory);
}
sub report {
    my ($directory, $tier) = @_;
    $tier //= 'compact';
    open my $input, '<', "$directory/scalar_chain_qualification_${backend}_${tier}.tsv.status.json" or die $!;
    local $/;
    return decode_json(<$input>);
}

for my $args (['--gtest_list_tests'], ['--gtest_repeat=0'], ['--gtest_filter=NoSuchQualification.*']) {
    my ($status, $directory) = run('--backend', $backend, @$args);
    is($status, 0, 'GoogleTest inspection/no-op succeeds');
    my $report = report($directory);
    is($report->{status}, 'not-run', 'no-op is not a successful qualification');
    is($report->{tests_run}, 0, 'counts actual executions');
    is($report->{advisories}, 0, 'no-op cannot report an advisory');
}
{
    local $ENV{GTEST_REPEAT} = '0';
    my ($status, $directory) = run('--backend', $backend);
    is($status, 0, 'zero repeat from environment');
    is(report($directory)->{status}, 'not-run', 'inherited no-op is not qualification');
}
for my $args ([], ['--backend', 'unknown'], ['--backend', $backend, '--extende']) {
    my ($status, $directory) = run(@$args);
    is($status >> 8, 2, 'invalid driver invocation fails before qualification');
    ok(!-e "$directory/scalar_chain_qualification_${backend}_compact.tsv.status.json", 'no misleading completion report');
}
{
    my ($status, $directory) = run('--backend', $backend, '--gtest_filter=ScalarChainQualification.FinitePhysicalGreenFunction',
                                 '--gtest_repeat=2');
    is($status, 0, 'selected backend can complete repeated qualification');
    my $report = report($directory);
    is($report->{status}, 'passed', 'completed qualification is marked passed');
    is($report->{tests_run}, 2, 'repeated executions counted rather than selected tests');
    is($report->{failed_tests}, 0, 'no failed executions');
    is($report->{skipped_tests}, 0, 'no skipped executions');
    is($report->{advisories}, 0, 'ordinary passes have no advisories');
}
if ($extended) {
    my ($status, $directory) = run('--backend', $backend, '--extended', '--gtest_filter=ScalarChainQualification.SameStarMatrix');
    is($status, 0, 'extended matrix succeeds with only the known accuracy miss advisory');
    my $report = report($directory, 'extended');
    open my $input, '<', "$directory/scalar_chain_qualification_${backend}_extended.tsv" or die $!;
    chomp(my $header = <$input>);
    my @columns = split /\t/, $header;
    my ($advisories, $known, $remaining) = (0, 0, 0);
    while (my $line = <$input>) {
        chomp $line;
        my %row;
        @row{@columns} = split /\t/, $line;
        isnt($row{status}, 'fail', 'no blocking row failed');
        if ($row{metric} eq 'known_gap_hop') {
            ++$known;
            is($backend, 'rkpw', 'advisory policy never covers legacy');
            is($row{case}, 'gap_large_lambda', 'advisory covers only the known case');
            is($row{gated}, 0, 'known hopping accuracy is non-blocking');
            is($row{hop_index}, 19, 'only hopping 19 is advisory');
            is($row{budget}, 2e-12, 'nominal accuracy target is unchanged');
            is($row{status}, $row{max_relative_hop} > $row{budget} ? 'advisory' : 'pass', 'measured exceedance is reported');
        } elsif ($row{case} eq 'gap_large_lambda' && $row{metric} eq 'reconstruction') {
            ++$remaining;
            is($row{gated}, 1, 'other coefficients remain blocking');
            is($row{status}, 'pass', 'remaining coefficient checks passed');
            isnt($row{hop_index}, 19, 'advisory hop is split out of gated maximum') if $backend eq 'rkpw';
        }
        if ($row{status} eq 'advisory') {
            ++$advisories;
            is($row{metric}, 'known_gap_hop', 'no unrelated metric can be advisory');
        }
    }
    close $input;
    is($known, $backend eq 'rkpw' ? 3 : 0, 'one known-hop row per RKPW frontend');
    is($remaining, 3, 'all three frontends retain blocking checks');
    is($report->{advisories}, $advisories, 'completion counts advisory rows');
    is($report->{status}, $advisories ? 'passed-with-advisories' : 'passed', 'completion distinguishes advisory success');
    is($report->{tests_run}, 1, 'extended matrix ran');
    is($report->{failed_tests}, 0, 'no failed executions');
    is($report->{skipped_tests}, 0, 'no skipped executions');
    open my $log, '<', "$directory/log" or die $!;
    my $warnings = grep { /^ADVISORY:/ } <$log>;
    close $log;
    is($warnings, $advisories, 'stdout exposes every advisory');
}
done_testing();

#!/usr/bin/env perl
# Compare text files containing numeric data with round-off differences.
# Part of "NRG Ljubljana"
# (C) Rok Zitko, rok.zitko@ijs.si, 2006-2026

use strict;
use warnings;

use Getopt::Long qw(GetOptions);

my $ATOL = 1e-12;
my $RTOL = 1e-5;
my $MAX_DIAGNOSTICS = 10;
my $ignore_signs = 0;
my $check_comments = 0;
my ($atol_arg, $rtol_arg, $first_atol_arg, $first_rtol_arg);

sub is_finite {
    my ($value) = @_;
    return $value == $value && $value - $value == 0;
}

sub fail_usage {
    my ($message) = @_;
    print STDERR "$message\n" if defined $message;
    print STDERR "Usage: mycomp.pl [--ignore-signs] [--check-comments] [--atol VALUE] [--rtol VALUE] [--first-atol VALUE] [--first-rtol VALUE] FILE1 FILE2\n";
    exit 2;
}

sub fail_input {
    my ($message) = @_;
    print STDERR "$message\n";
    exit 2;
}

GetOptions(
    'ignore-signs'  => \$ignore_signs,
    'check-comments' => \$check_comments,
    'atol=s'        => \$atol_arg,
    'rtol=s'        => \$rtol_arg,
    'first-atol=s'  => \$first_atol_arg,
    'first-rtol=s'  => \$first_rtol_arg,
) or fail_usage();

@ARGV == 2 or fail_usage('Two file arguments are required.');
my ($fn1, $fn2) = @ARGV;

my $FLOAT = qr/[+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?/;

sub finite_number {
    my ($text) = @_;
    return (0, undef, 'not a number') unless $text =~ /\A$FLOAT\z/;

    my $value = 0 + $text;
    return (0, undef, 'not finite') unless is_finite($value);

    # Reject values that silently underflow during Perl numeric conversion.
    my $mantissa = $text;
    $mantissa =~ s/[Ee].*\z//;
    $mantissa =~ s/[^0-9]//g;
    if ($value == 0 && $mantissa =~ /[1-9]/) {
        return (0, undef, 'not representable');
    }
    return (1, $value, undef);
}

sub parse_tolerance {
    my ($name, $text, $default) = @_;
    return $default unless defined $text;
    my ($ok, $value) = finite_number($text);
    fail_usage("$name must be a finite, nonnegative number") unless $ok && $value >= 0;
    return $value;
}

$ATOL = parse_tolerance('--atol', $atol_arg, $ATOL);
$RTOL = parse_tolerance('--rtol', $rtol_arg, $RTOL);
my $FIRST_ATOL = parse_tolerance('--first-atol', $first_atol_arg, $ATOL);
my $FIRST_RTOL = parse_tolerance('--first-rtol', $first_rtol_arg, $RTOL);

for my $filename ($fn1, $fn2) {
    -f $filename or fail_input("$filename is not a regular file.");
    -r $filename or fail_input("$filename is not readable.");
    -s $filename or fail_input("$filename is empty.");
}

open(my $F1, '<', $fn1) or fail_input("Can't open $fn1 for reading: $!");
open(my $F2, '<', $fn2) or fail_input("Can't open $fn2 for reading: $!");
print "Comparing $fn1 and $fn2.\n";

sub trim {
    my ($text) = @_;
    $text =~ s/^\s+//;
    $text =~ s/\s+$//;
    return $text;
}

sub next_record {
    my ($fh, $filename, $state) = @_;
    while (1) {
        $! = 0;
        my $line = <$fh>;
        if (!defined $line) {
            fail_input("Error reading $filename: $!") if $!;
            return undef;
        }
        $state->{physical_line}++;
        $line =~ s/\r?\n\z//;
        $line = trim($line);
        next if $line eq '';
        if ($line =~ /^#/) {
            next unless $check_comments;
            return [$line, $state->{physical_line}, 1];
        }
        $state->{data_records}++;
        return [$line, $state->{physical_line}, 0];
    }
}

sub parse_field {
    my ($token, $filename, $line_number, $field_number) = @_;
    my ($label, $payload) = ('', $token);
    if ($token =~ /\A([^=\s]+)=(.*)\z/) {
        ($label, $payload) = ($1, $2);
    }

    my ($ok, $value, $reason) = finite_number($payload);
    if ($ok) {
        return {label => $label, kind => 'number', value => $value, raw => $token};
    }
    if ($payload =~ /\A$FLOAT\z/) {
        fail_input("$filename:$line_number: field $field_number is $reason: '$token'.");
    }

    if ($payload =~ /\A\(($FLOAT),($FLOAT)\)\z/) {
        my ($real_text, $imag_text) = ($1, $2);
        my ($real_ok, $real, $real_reason) = finite_number($real_text);
        my ($imag_ok, $imag, $imag_reason) = finite_number($imag_text);
        fail_input("$filename:$line_number: field $field_number is invalid: $real_reason in '$token'.") unless $real_ok;
        fail_input("$filename:$line_number: field $field_number is invalid: $imag_reason in '$token'.") unless $imag_ok;
        return {label => $label, kind => 'complex', value => [$real, $imag], raw => $token};
    }

    if ($payload =~ /\A[+-]?(?:nan(?:\(.*\))?|inf(?:inity)?)\z/i
        || $payload =~ /\A[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))[Ee][+-]?\z/
        || $payload =~ /\A\([^\s]*,[^\s]*\).*\z/) {
        fail_input("$filename:$line_number: field $field_number is not a finite representable number: '$token'.");
    }

    return {label => $label, kind => 'text', value => $payload, raw => $token};
}

sub parse_record {
    my ($record, $filename) = @_;
    my ($line, $line_number) = @$record;
    my @tokens = split(/\s+/, $line);
    my @fields;
    for my $index (0 .. $#tokens) {
        push @fields, parse_field($tokens[$index], $filename, $line_number, $index + 1);
    }
    return \@fields;
}

sub close_number {
    my ($left, $right, $atol, $rtol) = @_;
    my $difference = abs($left - $right);
    return 0 unless is_finite($difference);
    my $scale = abs($left) > abs($right) ? abs($left) : abs($right);
    return $difference <= $atol + $rtol * $scale;
}

sub fields_match {
    my ($left, $right, $atol, $rtol) = @_;
    return 0 unless $left->{label} eq $right->{label};
    return 0 unless $left->{kind} eq $right->{kind};

    if ($left->{kind} eq 'number') {
        return close_number(abs($left->{value}), abs($right->{value}), $atol, $rtol) if $ignore_signs;
        return close_number($left->{value}, $right->{value}, $atol, $rtol);
    }
    if ($left->{kind} eq 'complex') {
        my ($lr, $li) = @{$left->{value}};
        my ($rr, $ri) = @{$right->{value}};
        my $same = close_number($lr, $rr, $atol, $rtol) && close_number($li, $ri, $atol, $rtol);
        return $same unless $ignore_signs;
        my $opposite = close_number($lr, -$rr, $atol, $rtol) && close_number($li, -$ri, $atol, $rtol);
        return $same || $opposite;
    }
    return $left->{value} eq $right->{value};
}

my %state1 = (physical_line => 0, data_records => 0);
my %state2 = (physical_line => 0, data_records => 0);
my $records = 0;
my $differences = 0;
my $shown = 0;

while (1) {
    my $record1 = next_record($F1, $fn1, \%state1);
    my $record2 = next_record($F2, $fn2, \%state2);
    last if !defined($record1) && !defined($record2);
    $records++;

    if (!defined($record1) || !defined($record2)) {
        print "Comparable record counts differ.\n" if $shown++ < $MAX_DIAGNOSTICS;
        $differences++;
        last;
    }

    my $fields1 = parse_record($record1, $fn1);
    my $fields2 = parse_record($record2, $fn2);
    if (@$fields1 != @$fields2) {
        print "$fn1:$record1->[1] and $fn2:$record2->[1] have different field counts ("
              . scalar(@$fields1) . " vs " . scalar(@$fields2) . ").\n"
            if $shown++ < $MAX_DIAGNOSTICS;
        $differences++;
        next;
    }

    for my $index (0 .. $#$fields1) {
        my ($atol, $rtol) = $index == 0 ? ($FIRST_ATOL, $FIRST_RTOL) : ($ATOL, $RTOL);
        next if fields_match($fields1->[$index], $fields2->[$index], $atol, $rtol);
        print "$fn1:$record1->[1]:" . ($index + 1) . " '$fields1->[$index]{raw}' != "
              . "$fn2:$record2->[1]:" . ($index + 1) . " '$fields2->[$index]{raw}'\n"
            if $shown++ < $MAX_DIAGNOSTICS;
        $differences++;
    }
}

close($F1) or fail_input("Can't finish reading $fn1: $!");
close($F2) or fail_input("Can't finish reading $fn2: $!");
fail_input("$fn1 contains no data records.") unless $state1{data_records};
fail_input("$fn2 contains no data records.") unless $state2{data_records};

if ($differences) {
    print "$records comparable records; $differences differences.\n";
    print "Only the first $MAX_DIAGNOSTICS differences are shown.\n" if $differences > $MAX_DIAGNOSTICS;
    exit 1;
}

exit 0;

#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(isfinite);

@ARGV == 4 or die "Usage: check_scalar.pl expected-file actual-file abs-tolerance rel-tolerance\n";
my ($expected_file, $actual_file, $abs_tolerance, $rel_tolerance) = @ARGV;
my $number = qr/[+-]?(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee][+-]?\d+)?/;

sub read_scalar {
  my ($path) = @_;
  open(my $input, '<', $path) or die "Can't open $path: $!\n";
  local $/;
  my $text = <$input>;
  close($input) or die "Can't close $path: $!\n";
  defined($text) && $text =~ /\A($number)\n\z/
    or die "$path does not contain exactly one finite numeric scalar followed by a newline\n";
  my $value = 0.0 + $1;
  isfinite($value) or die "$path does not contain a finite numeric scalar\n";
  return $value;
}

$abs_tolerance =~ /\A$number\z/ && $abs_tolerance >= 0.0
  or die "Invalid absolute tolerance: $abs_tolerance\n";
$rel_tolerance =~ /\A$number\z/ && $rel_tolerance >= 0.0
  or die "Invalid relative tolerance: $rel_tolerance\n";
isfinite(0.0 + $abs_tolerance) && isfinite(0.0 + $rel_tolerance)
  or die "Tolerances must be finite\n";

my $expected = read_scalar($expected_file);
my $actual = read_scalar($actual_file);
my $scale = abs($expected) > abs($actual) ? abs($expected) : abs($actual);
my $limit = $abs_tolerance + $rel_tolerance * $scale;
my $difference = abs($actual - $expected);
if ($difference > $limit) {
  die "$actual_file: $actual differs from $expected by $difference (limit $limit)\n";
}

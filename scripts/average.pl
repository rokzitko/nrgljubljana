#!/usr/bin/env perl
# Spectral function averaging
# Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Oct 2008
# $Id: average.pl,v 1.1 2008/10/08 08:11:40 rzitko Exp rzitko $

use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(dirname);
use lib $Bin, dirname(__FILE__);
use TableIO qw(atomic_write data_rows parse_number read_records scaled_mean);

$|=1;

sub avg {
    my ($Nz, $suffix, $output) = @_;
    defined($Nz) or die "avg: missing Nz\n";
    my $file_count = parse_number("$Nz", 'avg: Nz');
    $file_count > 0 && $file_count == int($file_count)
        or die "avg: Nz must be a positive integer\n";
    defined($suffix) && defined($output)
        or die "avg: suffix and output are required\n";

    my @filenames = map { "$_$suffix" } 1 .. $file_count;
    my @records = map { read_records($_, columns => 2) } @filenames;
    my @tables = map { data_rows($_) } @records;
    @{$tables[0]} or die "$filenames[0]: expected at least one data row\n";

    my $row_count = @{$tables[0]};
    for my $file (1 .. $#tables) {
        my $found = @{$tables[$file]};
        $found == $row_count
            or die "$filenames[$file]: expected $row_count data rows "
                . "to match $filenames[0], found $found\n";
        for my $row (0 .. $row_count - 1) {
            my ($reference, $current) = ($tables[0]->[$row], $tables[$file]->[$row]);
            my ($x, $other_x) = ($reference->{values}[0], $current->{values}[0]);
            $x == $other_x or die "$current->{path}:$current->{line}: column 1 ($other_x) "
                . "does not match $reference->{path}:$reference->{line} ($x)\n";
        }
    }
    my @means = map {
        my $row = $_;
        scaled_mean([map { $_->[$row]{values}[1] } @tables], "avg: row " . ($row + 1));
    } 0 .. $row_count - 1;

    my @results;
    for my $row (0 .. $row_count - 1) {
        push @results, [$tables[0]->[$row]{tokens}[0], $means[$row]];
    }

    my ($content, $row) = ('', 0);
    for my $record (@{$records[0]}) {
        $content .= $record->{type} eq 'data'
            ? join(' ', @{$results[$row++]}) . "\n"
            : $record->{raw};
    }

    print "$suffix -> $output (Nz=$Nz)\n";
    print "$_ " for 1 .. $file_count;
    atomic_write($output, $content);
    print "\n";
}

1;

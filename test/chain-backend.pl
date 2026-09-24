#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempfile);
use Cwd qw(abs_path);

Getopt::Long::Configure('require_order', 'no_auto_abbrev');
my ($backend, $mode, $input, $output, $source, $work, @params);
GetOptions('backend=s' => \$backend, 'mode=s' => \$mode,
           'input=s' => \$input, 'output=s' => \$output,
           'source=s' => \$source, 'work=s' => \$work, 'param=s@' => \@params)
    or die "Invalid chain-backend arguments\n";
die "Expected --backend legacy or rkpw\n" unless defined($backend) && $backend =~ /\A(?:legacy|rkpw)\z/;
die "Expected --mode initializer, tool or runtime\n" unless defined($mode) && $mode =~ /\A(?:initializer|tool|runtime)\z/;
my $method = $backend eq 'legacy' ? 'lanczos' : 'rkpw';

sub stage {
    my ($from, $to) = @_;
    open my $in, '<', $from or die "$from: $!\n";
    my @lines = <$in>;
    close $in or die "$from: $!\n";
    my ($start, $end, %found);
    my $in_param = 0;
    for my $i (0 .. $#lines) {
        if ($lines[$i] =~ /^\s*\[([^\]]+)\]\s*(?:#.*)?$/) {
            $end = $i if $in_param;
            $in_param = $1 eq 'param';
            if ($in_param) {
                die "$from: repeated [param] block\n" if defined $start;
                $start = $i;
                # The tool parser finds the literal header, unlike Mathematica.
                $lines[$i] = "[param]\n";
            }
        } elsif ($in_param && $lines[$i] =~ /^\s*(tri|tridiag_method)\s*=\s*([^#\r\n]*)(?:#.*)?$/) {
            my ($key, $value) = ($1, $2);
            $value =~ s/\s+$//;
            die "$from: duplicate $key in [param]\n" if exists $found{$key};
            $found{$key} = {index => $i, value => $value};
        }
    }
    die "$from: missing [param] block\n" unless defined $start;
    $end = scalar @lines unless defined $end;
    my %settings = (tridiag_method => $method);
    if ($mode eq 'initializer') {
        my $tri = exists($found{tri}) ? $found{tri}{value} : 'old';
        if ($tri eq 'old' || $tri eq 'rkpw') {
            $settings{tri} = $backend eq 'legacy' ? 'old' : 'rkpw';
        } elsif ($tri eq 'cpp' || $tri eq 'none') {
            $settings{tri} = $tri;
        } else {
            die "$from: initializer mode cannot select scalar backend for tri=$tri\n";
        }
    }
    my @insert;
    for my $key (sort keys %settings) {
        my $line = "$key=$settings{$key}\n";
        if (exists $found{$key}) {
            $lines[$found{$key}{index}] = $line;
        } else {
            push @insert, $line;
        }
    }
    $lines[$end - 1] .= "\n" if $end > 0 && $lines[$end - 1] !~ /\n\z/;
    splice @lines, $end, 0, @insert;
    # Replace the output itself, never follow an old parameter-file symlink.
    my ($out, $temporary) = tempfile('.chain-param-XXXXXX', DIR => dirname($to), UNLINK => 1);
    print {$out} @lines or die "$temporary: $!\n";
    close $out or die "$temporary: $!\n";
    rename $temporary, $to or die "Cannot publish $to: $!\n";
}

if (defined($input) || defined($output)) {
    die "Use --input FILE --output FILE without a command\n"
        unless defined($input) && defined($output) && !defined($source) && !defined($work) && !@params && !@ARGV;
    stage($input, $output);
    exit 0;
}

die "Use --source DIR --work DIR [--param FILE ...] -- COMMAND ...\n"
    unless defined($source) && defined($work) && @ARGV;
$source = abs_path($source) // die "Missing fixture directory\n";
$work = abs_path($work) // die "Missing test working directory\n";
die "Refusing to run in source fixture directory\n" if $source eq $work;
@params = ('param') unless @params;
for my $param (@params) {
    die "Expected fixture-relative parameter filename\n"
        if File::Spec->file_name_is_absolute($param) || $param =~ m{(?:\A|/)\.\.(?:/|\z)};
    stage("$source/$param", "$work/$param");
}
chdir $work or die "$work: $!\n";
$ENV{CHAIN_TEST_BACKEND} = $backend;
print "Scalar-chain regression backend: $backend ($method), mode: $mode\n";
$| = 1;
exec {$ARGV[0]} @ARGV;
die "Cannot execute $ARGV[0]: $!\n";

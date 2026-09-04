package PhysicalOutput;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(is_physical_output is_spectral_output is_text_physical_output);

my %EXACT_OUTPUT = map { $_ => 1 } qw(
    absolute_energies.dat
    annotated.dat
    custom
    customfdm
    energies.nrg
    raw-dm.h5
    raw.h5
    report.nrg
    states.nrg
    subspaces.dat
    td
    tdfdm
);

sub is_spectral_output {
    my ($name) = @_;
    return $name =~ m{\A(?:chit|corr|gt|i1t|i2t|orbspin|spec|specq|spin)_[^/\\]+\.(?:bin|dat)\z};
}

sub is_physical_output {
    my ($name) = @_;
    return $EXACT_OUTPUT{$name} || is_spectral_output($name);
}

sub is_text_physical_output {
    my ($name) = @_;
    return is_physical_output($name) && $name !~ /\.(?:bin|h5)\z/;
}

1;

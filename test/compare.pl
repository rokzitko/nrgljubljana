#!/usr/bin/env perl
# Compare reference files with results in an explicit or current directory.

use strict;
use warnings;

use FindBin qw($RealBin);
use File::Spec;
use Getopt::Long qw(GetOptions);
use lib $RealBin;
use PhysicalOutput qw(is_physical_output is_spectral_output is_text_physical_output);

my $VALIDATION_MANIFEST = '.physical-outputs';
my $actual_dir = '.';
my $ignore_signs = 0;
my $strict = 0;
my @exclude;

sub is_finite {
    my ($value) = @_;
    return $value == $value && $value - $value == 0;
}

sub usage_error {
    my ($message) = @_;
    print STDERR "$message\n" if defined $message;
    print STDERR "Usage: compare.pl [--actual DIR] [--exclude FILE] [--ignore-signs] [--strict] REF_DIR\n";
    exit 2;
}

sub input_error {
    my ($message) = @_;
    print STDERR "$message\n";
    exit 2;
}

GetOptions(
    'actual=s'      => \$actual_dir,
    'exclude=s@'    => \@exclude,
    'ignore-signs'  => \$ignore_signs,
    'strict'        => \$strict,
) or usage_error();

@ARGV == 1 or usage_error('One reference directory is required.');
my $reference_dir = File::Spec->rel2abs($ARGV[0]);
$actual_dir = File::Spec->rel2abs($actual_dir);
-d $reference_dir or input_error("Reference directory $reference_dir does not exist.");
-d $actual_dir or input_error("Result directory $actual_dir does not exist.");

my %excluded = map { $_ => 1 } @exclude;

sub directory_entries {
    my ($directory) = @_;
    opendir(my $handle, $directory) or input_error("Can't read directory $directory: $!");
    my @entries = sort grep { $_ ne '.' && $_ ne '..' } readdir($handle);
    closedir($handle) or input_error("Can't finish reading directory $directory: $!");
    return @entries;
}

my @reference_entries = directory_entries($reference_dir);
my @reference_files;
for my $name (@reference_entries) {
    next if $name eq $VALIDATION_MANIFEST;
    next if $excluded{$name};
    my $path = File::Spec->catfile($reference_dir, $name);
    -f $path or input_error("Unsupported reference entry $path; only regular files are allowed.");
    push @reference_files, $name;
}
my %reference_file = map { $_ => 1 } @reference_files;
my %validated_file;
my $manifest_path = File::Spec->catfile($reference_dir, $VALIDATION_MANIFEST);
if (-e $manifest_path) {
    -f $manifest_path or input_error("Validation manifest $manifest_path is not a regular file.");
    open(my $manifest, '<', $manifest_path) or input_error("Can't read $manifest_path: $!");
    my $line_number = 0;
    while (my $line = <$manifest>) {
        $line_number++;
        $line =~ s/\r?\n\z//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;
        my ($validator, $name, @extra) = split(/\s+/, $line);
        input_error("$manifest_path:$line_number: expected VALIDATOR FILE.")
            unless defined($validator) && defined($name) && !@extra;
        input_error("$manifest_path:$line_number: unsupported validator '$validator'.")
            unless $validator =~ /\A(?:binary-complex|binary-real|hdf5|subspaces)\z/;
        input_error("$manifest_path:$line_number: '$name' is not a physical output name.")
            unless is_physical_output($name);
        input_error("$manifest_path:$line_number: validator '$validator' is incompatible with '$name'.")
            if (($validator =~ /\Abinary-/ && ($name !~ /\.bin\z/ || !is_spectral_output($name)))
                || ($validator eq 'hdf5' && $name !~ /\.h5\z/)
                || ($validator eq 'subspaces' && $name ne 'subspaces.dat'));
        next if $excluded{$name};
        input_error("$manifest_path:$line_number: duplicate output '$name'.")
            if $reference_file{$name} || $validated_file{$name};
        $validated_file{$name} = $validator;
    }
    close($manifest) or input_error("Can't finish reading $manifest_path: $!");
}

(@reference_files || keys(%validated_file))
    or input_error("Reference directory $reference_dir contains no files or validations.");

my $differences = 0;
my $double_size = length(pack('d', 0));

sub validate_binary {
    my ($path, $columns) = @_;
    -f $path or return "Result $path does not exist or is not a regular file.";
    -s $path or return "Binary result $path is empty.";
    open(my $input, '<:raw', $path) or input_error("Can't read binary result $path: $!");
    my $contents = do { local $/; <$input> };
    defined($contents) or input_error("Can't read binary result $path: $!");
    close($input) or input_error("Can't finish reading binary result $path: $!");
    my $record_size = $columns * $double_size;
    return "Binary result $path has a partial record."
        if length($contents) % $record_size;
    my @values = unpack('d*', $contents);
    return "Binary result $path contains no records." unless @values;
    for my $index (0 .. $#values) {
        return "Binary result $path has a non-finite value at record "
            . int($index / $columns) . ', field ' . ($index % $columns + 1) . '.'
            unless is_finite($values[$index]);
    }
    return undef;
}

sub validate_hdf5 {
    my ($path) = @_;
    -f $path or return "Result $path does not exist or is not a regular file.";
    -s $path or return "HDF5 result $path is empty.";
    my ($h5dump) = grep { -f $_ && -x _ }
        map { File::Spec->catfile($_, 'h5dump') } File::Spec->path();
    defined($h5dump) or input_error('h5dump is required to validate HDF5 output.');
    open(my $dump, '-|', $h5dump, '-y', '-w', '0', $path)
        or input_error("Can't execute $h5dump for $path: $!");
    my $datasets = 0;
    my $nonfinite = 0;
    while (my $line = <$dump>) {
        $datasets++ if $line =~ /\bDATASET\b/;
        $nonfinite = 1
            if $line =~ /(?:^|[^A-Za-z])(?:[+-]?(?:nan|inf(?:inity)?))(?:[^A-Za-z]|$)/i;
    }
    my $closed = close($dump);
    if (!$closed) {
        input_error("$h5dump terminated while reading $path.") if $? == -1 || ($? & 127);
        return "h5dump could not read $path.";
    }
    return "HDF5 result $path contains no datasets." unless $datasets;
    return "HDF5 result $path contains a non-finite value." if $nonfinite;
    return undef;
}

sub validate_subspaces {
    my ($path) = @_;
    -f $path or return "Result $path does not exist or is not a regular file.";
    -s $path or return "Subspace result $path is empty.";
    open(my $input, '<', $path) or input_error("Can't read subspace result $path: $!");

    my ($active, $expected_rows, $rows, $iterations, $last_iteration, $invariant_arity)
        = (0, undef, 0, 0, undef, undef);
    my %seen_invariant;
    my $line_number = 0;
    while (my $line = <$input>) {
        $line_number++;
        $line =~ s/\r?\n\z//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '';

        if ($line =~ /\AIteration (\d+)\z/) {
            my $iteration = 0 + $1;
            if ($active && (!defined($expected_rows) || $rows != $expected_rows)) {
                close($input);
                return "Subspace result $path has an incomplete iteration before line $line_number.";
            }
            if (defined($last_iteration) && $iteration <= $last_iteration) {
                close($input);
                return "Subspace result $path has a non-increasing iteration at line $line_number.";
            }
            $last_iteration = $iteration;
            ($active, $expected_rows, $rows) = (1, undef, 0);
            %seen_invariant = ();
            $iterations++;
        } elsif ($line =~ /\Alen_dm=(\d+)\z/) {
            if (!$active || defined($expected_rows) || $1 == 0) {
                close($input);
                return "Subspace result $path has an invalid len_dm at line $line_number.";
            }
            $expected_rows = 0 + $1;
        } elsif ($line =~ /\AI=([+-]?\d+(?:\s+[+-]?\d+)*)\s+kept=(\d+)\s+total=(\d+)\z/) {
            my ($invariant, $kept, $total) = ($1, 0 + $2, 0 + $3);
            my $arity = scalar split(/\s+/, $invariant);
            $invariant_arity = $arity unless defined $invariant_arity;
            if (!$active || !defined($expected_rows) || $kept > $total
                || $seen_invariant{$invariant}++ || $arity != $invariant_arity
                || ++$rows > $expected_rows) {
                close($input);
                return "Subspace result $path has an invalid subspace at line $line_number.";
            }
        } else {
            close($input);
            return "Subspace result $path has a malformed record at line $line_number.";
        }
    }
    close($input) or input_error("Can't finish reading subspace result $path: $!");
    return "Subspace result $path contains no iterations." unless $iterations;
    return "Subspace result $path has an incomplete final iteration."
        unless defined($expected_rows) && $rows == $expected_rows;
    return undef;
}

for my $name (sort keys %validated_file) {
    my $path = File::Spec->catfile($actual_dir, $name);
    my $validator = $validated_file{$name};
    my $message;
    if ($validator eq 'binary-real') {
        $message = validate_binary($path, 2);
    } elsif ($validator eq 'binary-complex') {
        $message = validate_binary($path, 3);
    } elsif ($validator eq 'hdf5') {
        $message = validate_hdf5($path);
    } else {
        $message = validate_subspaces($path);
    }
    if (defined $message) {
        print STDERR "$message\n";
        $differences++;
    } else {
        print "Validated $name ($validator)\n";
    }
}

for my $name (@reference_files) {
    my $expected = File::Spec->catfile($reference_dir, $name);
    my $actual = File::Spec->catfile($actual_dir, $name);
    if (!-f $actual) {
        print STDERR "Result $actual does not exist or is not a regular file.\n";
        $differences++;
        next;
    }
    if ($name =~ /\.(?:bin|h5)\z/) {
        input_error("Reference $expected requires a semantic binary-file validator.");
    }

    print "Comparing $name\n";
    my @command = ($^X, "$RealBin/mycomp.pl");
    push @command, '--ignore-signs' if $ignore_signs;
    push @command, '--check-comments' if $strict && is_text_physical_output($name);
    if ($strict && is_spectral_output($name)) {
        push @command, '--atol', '1e-12', '--rtol', '2e-2',
            '--first-atol', '1e-12', '--first-rtol', '1e-5';
    } elsif ($strict && is_physical_output($name)) {
        push @command, '--atol', '1e-9', '--rtol', '1e-5';
    }
    push @command, $expected, $actual;
    my $status = system(@command);
    input_error("Can't execute mycomp.pl: $!") if $status == -1;
    input_error("mycomp.pl terminated by signal " . ($status & 127)) if $status & 127;
    my $exit = $status >> 8;
    $differences++ if $exit == 1;
    input_error("mycomp.pl could not compare $name (exit $exit).") if $exit > 1;
}

if ($strict) {
    for my $name (directory_entries($actual_dir)) {
        next if $excluded{$name} || $reference_file{$name} || $validated_file{$name};
        next unless is_physical_output($name);
        my $path = File::Spec->catfile($actual_dir, $name);
        print STDERR "Physical result $path has no reference.\n";
        $differences++;
    }
}

if ($differences) {
    print STDERR "Comparison failed with $differences unmatched result(s).\n";
    exit 1;
}

print "OK!\n";
exit 0;

#!/usr/bin/env perl
# The spin expectation values of the rotated run against those of the diagonal one.
#
#   check_rotation.pl PHI [TOLERANCE]      (radians; default tolerance 1e-8, relative to |<S_z>| of the diagonal run)
#
# The rotated problem is the diagonal one seen in a spin basis turned by phi about y, bath and field together, so its
# spin expectation value is the rotation of the other:
#
#   <S_z>_rot = cos(phi) <S_z>_diag,     <S_x>_rot = sin(phi) <S_z>_diag.
#
# This is what the off-diagonal coefficient sets have to get right. Writing the chain as its own transpose swaps the
# UPDO and DOUP sets, which is the same as turning by -phi: <S_z> survives it and <S_x> changes sign. The test
# therefore reports both signs, so a failure says which convention the data follow rather than only that they differ.

use strict;
use warnings;

my $phi = shift // die "usage: check_rotation.pl PHI [TOLERANCE]\n";
my $tolerance = shift // 1e-8;

# The custom file: two comment lines, the second naming the columns, then one row per temperature.
sub read_custom {
    my $path = shift;
    open(my $in, "<", $path) or die "check_rotation: cannot read $path: $!\n";
    my (@names, @rows);
    while (defined(my $line = <$in>)) {
        chomp $line;
        if ($line =~ /^#/) {
            my @fields = split ' ', $line;
            shift @fields;                       # the '#'
            @names = @fields if @fields && $fields[0] eq "T";
            next;
        }
        next unless $line =~ /\S/;
        push @rows, [ split ' ', $line ];
    }
    close $in;
    die "check_rotation: $path has no column names.\n" unless @names;
    die "check_rotation: $path has no data rows.\n" unless @rows;
    my %column;
    $column{ $names[$_] } = $_ for 0 .. $#names;
    return (\%column, \@rows);
}

my ($diagonal_columns, $diagonal_rows) = read_custom("diagonal/custom");
my ($rotated_columns,  $rotated_rows)  = read_custom("rotated/custom");

for my $name (qw(T SZd SXd)) {
    die "check_rotation: column $name is missing.\n"
      unless defined $diagonal_columns->{$name} && defined $rotated_columns->{$name};
}
die "check_rotation: the two runs have " . scalar(@$diagonal_rows) . " and " . scalar(@$rotated_rows) . " rows.\n"
  if @$diagonal_rows != @$rotated_rows;

my $cosine = cos($phi);
my $sine   = sin($phi);
my ($worst, $worst_flipped, $where, $scale) = (0, 0, "", 0);

for my $row (0 .. $#$diagonal_rows) {
    my $sz_diagonal = $diagonal_rows->[$row][ $diagonal_columns->{SZd} ];
    my $sz_rotated  = $rotated_rows->[$row][ $rotated_columns->{SZd} ];
    my $sx_rotated  = $rotated_rows->[$row][ $rotated_columns->{SXd} ];
    my $temperature = $diagonal_rows->[$row][ $diagonal_columns->{T} ];

    # The scale of the comparison is the quantity itself, not the rounding of a vanishing component.
    my $size = abs($sz_diagonal);
    $scale = $size if $size > $scale;
    next unless $size > 0;

    my $deviation = (abs($sz_rotated - $cosine * $sz_diagonal) + abs($sx_rotated - $sine * $sz_diagonal)) / $size;
    my $flipped   = (abs($sz_rotated - $cosine * $sz_diagonal) + abs($sx_rotated + $sine * $sz_diagonal)) / $size;
    $worst_flipped = $flipped if $flipped > $worst_flipped;
    if ($deviation > $worst) {
        $worst = $deviation;
        $where = sprintf "T=%.6g: <S_z> %.10g -> (%.10g, %.10g), expected (%.10g, %.10g)",
          $temperature, $sz_diagonal, $sz_rotated, $sx_rotated, $cosine * $sz_diagonal, $sine * $sz_diagonal;
    }
}

printf "rotation by phi=%.10g: largest deviation %.2e (with <S_x> flipped it would be %.2e), largest |<S_z>| = %.6g\n",
  $phi, $worst, $worst_flipped, $scale;

if ($worst <= $tolerance) {
    print "ok       the rotated run is the diagonal one turned by phi\n";
    exit 0;
}
print "ABOVE    $where\n";
print "         <S_x> matches the opposite sign, so the two off-diagonal sets are swapped: xi3 and xi4, zeta3 and zeta4.\n"
  if $worst_flipped <= $tolerance;
exit 1;

#!/usr/bin/env perl
# Process a 'data.in' template file: fill in the numerical values,
# diagonalize the Hamiltonian matrices, perform the transformations.
# Rok Zitko, rok.zitko@ijs.si, June 2009, May 2010

# CHANGE LOG
# 24.5.2010 - improved error reporting
#           - code cleanup, strictness

use strict;
use warnings;
use File::Copy;
use Math::Trig;

# Paths to tools used
my $matrix = "matrix";
my $diag = "diag -B";
my $unitary = "unitary -B";

my $verbose = 1;

my $fnin = "data.in";
open(F, "<$fnin") or die "Can't open $fnin.";

my $fnout = "data";
open(OUT, ">$fnout") or die "Can't open $fnout.";

my $mode = 0; # Current parser mode
my $factor = -1; # Rescale factor for energies
my $smallest = 1e9; # Lowest rescaled energy (Egs=smallest/factor)

my %dim; # Dimensions of invariant subspaces

my $subspaces = -1; # Number of different subspaces
my $Egs = -999; 

while (<F>) {
    print;
    
    if (/^#/) {
	if (/^# SCALE\s+(\S+)$/) {
	    # Extract the rescale factor
	    $factor = 1.0/$1;
	}
	# Copy comments verbatim
	print OUT;
	next;
    }
    
    if ($mode == 0) {
	if (!/(\d+)\s+(\d+)\s+(\d+)/) {
	    die "Parse error 1: $_";
	}
	my $channels = $1;
	my $impurities = $2;
	$subspaces = $3;
	print "ch=$channels im=$impurities sub=$subspaces\n";
	print OUT "$channels $impurities $subspaces\n";
	$mode = 1;
	next;
    }

    if ($mode == 1) {
	my @buffer;
	my $i;
	for ($i = 1; $i <= $subspaces; $i++) {
	    if ($i != 1) {
		$_ = <F>;
	    }
	    if (!/(\S+)\s+(\S+)/) {
		die "Parse error 2: $_";
	    }

	    my $qnspaces = "$1 $2";
	    my $qn = "$1.$2";
	    my $size;
	    chomp($size = <F>);
	    $dim{$qn} = $size; # store size!
	    print "[$i] $qn dim=$size\n";
	    
	    my $line1;
	    chomp($line1 = <F>);
	    $line1 =~ /DIAG\s+(\S+)$/;
	    my $fn = $1;

	    system("$matrix param $fn >ham");
	    if ($? != 0) {
		die "matrix call error";
	    }
	    copy("ham", "ham.$qn");

	    ($factor > 0) or die "no factor?";
	    system("$diag -q -s $factor -o val -O vec ham");

	    my $eig = `cat val`;
	    push(@buffer,  [ "$qnspaces\n$size\n", $eig ]);
	    if ($eig =~ /^([-0-9.]+)(\s|\n)/) {
		my $min = $1;
		if ($smallest > $min) { $smallest = $min; }
	    }
	    copy("vec", "vec.$qn");
	}
	$Egs=$smallest/$factor;
	print "E_gs=$Egs\n";

	@buffer = map { $_->[0].subtract($_->[1]) } @buffer;
	print OUT @buffer;

	$mode = 2;
	next;
    }

    if ($mode == 2) {
	if (!/f\s+(\d+)\s+(\d+)/) {
	    die "Parse error 3: $_";
	}
	print;
	print OUT;
	unitary();
	$mode = 3;
	next;
    }

    if ($mode == 3) {
	if (!/^e/) {
	    die "Parse error 5: $_";
	}
	my $egs = <F>;
	print OUT "e\n$Egs\n";
	$mode = 4;
	next;
    }

    if ($mode == 4) {
       if (/^[zZT]/) { # Discretization tables
	   $mode = 5;
	   next;
       }
       if (/^[spdvt]/) {
	   print;
	   print OUT;
	   unitary();
	   next;
       }
       die "Parse error 6: $_";
    }
}

my $adddisc = 1;

if ($adddisc) {

  print OUT "z\n";

  my $line = `wc xi.dat`;
  $line =~ /^\s*(\d+)/;
  my $nr = $1-1;
    
  print OUT "$nr\n";
  my $xi = `cat xi.dat`;
  print OUT "$xi";

  print OUT "$nr\n";
  my $zeta = `cat zeta.dat`;
  print OUT "$zeta";
}

sub subtract
{
    my $vals = shift;
    
    my @l = split(' ', $vals);
    @l = map { $_-$smallest } @l;
    
    return "@l\n";
}

sub unitary
{
    my $nr;
    chomp($nr = <F>);
    print OUT "$nr\n";
    
    my $i;
    for ($i = 1; $i <= $nr; $i++) {
	print "unitary i=$i\n";
	my $line;
	chomp($line = <F>);
	if (!($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/)) {
	    die "Parse error 4: $_";
	}

	my $qn1spaces = "$1 $2";
	my $qn1 = "$1.$2";
	my $dim1 = $dim{$qn1};
	my $qn2spaces = "$3 $4";
	my $qn2 = "$3.$4";
	my $dim2 = $dim{$qn2};
	
	print "$qn1 $dim1 $qn2 $dim2\n";
	
	open(O, ">mat") or die "Can't open mat for writing.";
	my $n;
	for ($n = 1; $n <= $dim1; $n++) {
	    if (0) {
	    chomp($line = <F>);
	    $line =~ s/Sqrt\[([^]]+)\]/sqrt($1)/g;
	    my @l = split(' ', $line);
	    @l = map(eval, @l);
	    @l = map { sprintf("%.18g", $_) } @l;
	    if ($verbose) {
	      print "$line\n";
	    }
	    if (scalar @l == 0) {
		print "warning: zero length: $line\n";
	    }
	    print O "@l\n";
	    }
	    
	    $line = <F>;
	    print O $line;
	
        } # for $n
        close(O);
    
        my $cmd = "$unitary -l -q -o mat.res vec.$qn1 mat vec.$qn2";
        system($cmd);
        if ($? != 0) {
	    die "unitary call error: $cmd";
	}
        my $res = `cat mat.res`;
	
        print OUT "$qn1spaces $qn2spaces\n$res";
    } # for $i
}

system "rm vec.* vec val mat mat.res ham ham.[-0-9]*";

#!/usr/bin/env perl
# Compare files which contain numeric data with round-uff errors.
# Lines are parsed and numbers are compaired in pairs up to some
# given accuracy.
# Part of "NRG Ljubljana"
# (C) Rok Zitko, rok.zitko@ijs.si, 2006-2011

# CHANGE LOG
# 24. 5. 2007 - EOF bugfix, getline() function
# 17. 6. 2009 - relaxed error triggering
# 9. 3. 2011 - two verbosity levels
# 29 3. 2011 - new ABSACCURACY setting   

use strict;
use warnings;

# Keep ABSACCURACY low to avoid false negative reports. In particular,          
# the round-off errors in saving results might differ from a system             
# to a system.
my $ABSACCURACY = 1e-5;
my $nrdiffshow = 10;

# verbose=0 = no output
# verbose=1 = minimal output
# verbose>1 = lots of output
my $verbose = 2;

(@ARGV == 2) or die "Two arguments required.\n";

my $fn1 = shift;
my $fn2 = shift;

open(my $F1, "<$fn1") or die "Can't open $fn1 for reading.\n";
open(my $F2, "<$fn2") or die "Can't open $fn2 for reading.\n";

if ($verbose >= 2) {
    print ("Comparing $fn1 and $fn2.\n");
}

# Trim leading and trailing whitespace
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Get line, chomp and trim it.
# Differences in CR/LF line ends and leading/trailing whitespace
# are ignored.
# Do so repeatedly, until a non-comment line is read.
sub getline($) {
    my $F = shift;
    my $line;
    do {
        $line = <$F>;
        chomp($line);
        $line = trim($line);
        $_ = $line;
    } until (!/^#/ || eof($F));
    
    if (eof($F)) {
        return "END OF FILE REACHED\n";
    }
    return $line;
}

# Is argument numeric? (C float)
sub isnumeric($) {
    my $string = shift;
    if ($string =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
        return 1;
    }
    return 0;
}  

my $error = 0;
my $cnt = 0; # line counter
my $cntmatch = 0; # identical matching counter
my $cntnumeric = 0; # numerical matching counter

while (!(eof($F1) || eof($F2))) { 
    $cnt++;
    chomp(my $line1 = getline($F1));
    chomp(my $line2 = getline($F2));
    $line1 =~ s/^\s//;
    $line2 =~ s/^\s//;
    $line1 =~ s/\s$//;
    $line2 =~ s/\s$//;
    
    if ($line1 eq $line2) {
        $cntmatch += 1; 
    } else {
        # If lines don't match, try comparing number by number
        my @numbers1 = split(/\s+/,$line1);
        my @numbers2 = split(/\s+/,$line2);
        my $len1 = @numbers1;
        my $len2 = @numbers2;
        if ($verbose >= 3) {
            print "@numbers1 [$len1]\n";
            print "@numbers2 [$len2]\n";
        }
        if ($len1 != $len2) {
            if ($verbose >= 2) {
                print "Lengths don't match, $len1 vs $len2\n";
            }
            $error++;
        } else {
            my $i;
            for ($i = 0; $i < $len1; $i++) {
                my $n1 = pop(@numbers1);
                my $n2 = pop(@numbers2);
                if (isnumeric($n1) && isnumeric($n2)) {
                    my $diff = abs($n1)-abs($n2);
                    if (abs($diff) > $ABSACCURACY) {
                        if ($verbose >= 2) {
                            print "$cnt : $n1 != $n2 (even ignoring signs)\n";
                        }
                        $error++;
                    }
                } else {
                    if ($verbose >= 2) {
                        print "Not numeric: $n1 vs. $n2\n";
                    }
                    if (!isnumeric($n1) && isnumeric($n2)) {
                        # If the second (reference) is numeric, but the
                        # first one isn't, than it's probably really
                        # an error.
                        $error++;
                    }
                }
            }
        }
        if ($error == 0) {
            $cntnumeric += 1;
        }
    } 
}

close($F1);
close($F2);

my $cntdiff = $cnt - ($cntmatch+$cntnumeric);

if ($verbose >= 1 && $error > 0) {
    print ("$cnt lines. $cntmatch identical. $cntnumeric numerically equivalent.");
    print (" $cntdiff differing.\n");
    print ("Error count=$error\n");
}

if ($verbose == 1) {
    if ($error == 0) {
        print "$fn1 and $fn2 match.\n";
    } else {
        print "Error: $fn1 and $fn2 different.\n";
    }
}

exit($error);

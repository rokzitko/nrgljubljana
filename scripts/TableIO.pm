package TableIO;

use strict;
use warnings;

use Exporter qw(import);
use File::Basename qw(dirname);
use File::Temp qw(tempfile);
use POSIX qw(isfinite);

our @EXPORT_OK = qw(atomic_write atomic_write_many data_rows finite_result parse_number
                    read_records scaled_mean scaled_mean_sem);

my $number = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;

sub parse_number {
    my ($text, $where) = @_;
    defined($text) && $text =~ /\A$number\z/
        or die "$where: invalid number '$text'\n";
    my $value = 0 + $text;
    isfinite($value) or die "$where: nonfinite number '$text'\n";
    if ($value == 0 && $text !~ /\A[+-]?(?:0+(?:\.0*)?|\.0+)(?:[eE][+-]?\d+)?\z/) {
        die "$where: number '$text' underflows to zero\n";
    }
    return $value;
}

sub finite_result {
    my ($value, $where) = @_;
    isfinite($value) or die "$where: arithmetic produced a nonfinite result\n";
    return $value;
}

sub scaled_mean {
    my ($values, $where) = @_;
    @$values or die "$where: can't average an empty sample\n";
    my $scale = 0;
    for my $value (@$values) {
        $scale = abs($value) if abs($value) > $scale;
    }
    return 0 if $scale == 0;

    my ($sum, $correction) = (0, 0);
    for my $value (@$values) {
        my $adjusted = $value / $scale - $correction;
        my $next = $sum + $adjusted;
        $correction = ($next - $sum) - $adjusted;
        $sum = $next;
    }
    return finite_result($scale * ($sum / @$values), $where);
}

sub scaled_mean_sem {
    my ($values, $where) = @_;
    @$values >= 2 or die "$where: at least two samples are required\n";
    my $scale = 0;
    for my $value (@$values) {
        $scale = abs($value) if abs($value) > $scale;
    }
    return (0, 0) if $scale == 0;

    my ($mean, $m2, $count) = (0, 0, 0);
    for my $value (@$values) {
        $count++;
        my $normalized = $value / $scale;
        my $delta = $normalized - $mean;
        $mean += $delta / $count;
        $m2 += $delta * ($normalized - $mean);
    }
    my $result = finite_result($scale * $mean, "$where mean");
    my $error = finite_result($scale * sqrt($m2 / ($count * ($count - 1))), "$where SEM");
    return ($result, $error);
}

sub read_records {
    my ($path, %options) = @_;
    my $fh;
    if ($path eq '-') {
        open($fh, '<&', \*STDIN) or die "Can't read standard input: $!\n";
    } else {
        open($fh, '<', $path) or die "Can't open $path for reading: $!\n";
    }

    my @records;
    my $line_number = 0;
    my $record_width;
    while (my $raw = <$fh>) {
        $line_number++;
        if ($raw =~ /^\s*$/) {
            push @records, { type => 'blank', raw => $raw };
            next;
        }
        if ($raw =~ /^\s*#/) {
            push @records, { type => 'comment', raw => $raw };
            next;
        }

        (my $line = $raw) =~ s/[\r\n]+\z//;
        my @tokens = split ' ', $line;
        my $where = "$path:$line_number";
        if (defined($options{columns}) && @tokens != $options{columns}) {
            die "$where: expected $options{columns} columns, found " . scalar(@tokens) . "\n";
        }
        if (defined($options{minimum_columns}) && @tokens < $options{minimum_columns}) {
            die "$where: expected at least $options{minimum_columns} columns, found " . scalar(@tokens) . "\n";
        }
        if ($options{fixed_columns}) {
            $record_width //= scalar(@tokens);
            @tokens == $record_width
                or die "$where: expected $record_width columns, found " . scalar(@tokens) . "\n";
        }

        my $numeric_columns = $options{numeric_columns} // scalar(@tokens);
        my @values;
        for my $index (0 .. $numeric_columns - 1) {
            push @values, parse_number($tokens[$index], "$where: column " . ($index + 1));
        }
        push @records, {
            type => 'data', raw => $raw, tokens => \@tokens, values => \@values,
            path => $path, line => $line_number,
        };
    }
    close($fh) or die "Can't finish reading $path: $!\n";
    return \@records;
}

sub data_rows {
    my ($records) = @_;
    return [ grep { $_->{type} eq 'data' } @$records ];
}

sub _stage {
    my ($path, $content, $mode_source) = @_;
    -l $path and die "Refusing to replace symlink $path\n";
    -e $path && !-f $path and die "Refusing to replace non-file $path\n";
    my ($fh, $temporary) = tempfile('.tmp.XXXXXX', DIR => dirname($path), UNLINK => 0);
    my $open = 1;
    eval {
        my @source_stat = stat(defined($mode_source) ? $mode_source : $path);
        if (@source_stat) {
            chmod($source_stat[2] & 07777, $temporary)
                or die "Can't set permissions on $temporary: $!\n";
        } else {
            my $mask = umask;
            umask($mask);
            chmod(0666 & ~$mask, $temporary)
                or die "Can't set permissions on $temporary: $!\n";
        }
        print {$fh} $content or die "Can't write $temporary: $!\n";
        close($fh) or die "Can't finish writing $temporary: $!\n";
        $open = 0;
        1;
    } or do {
        my $error = $@ || "Can't stage $path\n";
        close($fh) if $open;
        unlink($temporary);
        die $error;
    };
    return $temporary;
}

sub atomic_write {
    my ($path, $content, $mode_source) = @_;
    my $temporary;
    eval {
        $temporary = _stage($path, $content, $mode_source);
        rename($temporary, $path) or die "Can't replace $path: $!\n";
        undef $temporary;
        1;
    } or do {
        my $error = $@ || "Can't replace $path\n";
        unlink($temporary) if defined($temporary);
        die $error;
    };
}

sub atomic_write_many {
    my ($outputs) = @_;
    my (%temporary, %backup, @published);
    eval {
        for my $path (sort keys %$outputs) {
            $temporary{$path} = _stage($path, $outputs->{$path});
        }
        for my $path (sort keys %temporary) {
            next unless -e $path;
            my ($fh, $name) = tempfile('.backup.XXXXXX', DIR => dirname($path), UNLINK => 0);
            close($fh) or die "Can't prepare backup for $path: $!\n";
            unlink($name) or die "Can't prepare backup for $path: $!\n";
            rename($path, $name) or die "Can't back up $path: $!\n";
            $backup{$path} = $name;
        }
        for my $path (sort keys %temporary) {
            rename($temporary{$path}, $path) or die "Can't replace $path: $!\n";
            delete $temporary{$path};
            push @published, $path;
        }
        1;
    } or do {
        my $error = $@ || "Can't publish output files\n";
        unlink($_) for @published;
        for my $path (sort keys %backup) {
            rename($backup{$path}, $path)
                or $error .= "Can't restore $path from $backup{$path}: $!\n";
            delete $backup{$path};
        }
        unlink($_) for values %temporary;
        unlink($_) for values %backup;
        die $error;
    };
    unlink($_) or die "Can't remove output backup $_: $!\n" for values %backup;
}

1;

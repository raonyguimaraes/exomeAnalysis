use warnings;
use strict;

my $file = shift;

my $length;
open(IN, "<$file");
while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/>/) {
	chomp(my $seq = <IN>);
	$length += length($seq);
    }
}
close(IN);

print $file, "\t", $length, "\n";

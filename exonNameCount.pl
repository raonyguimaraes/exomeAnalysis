use warnings;
use strict;

my @a = <*fa>;

my %h;
foreach my $a (@a) {
    open(IN, "<$a");
    my %tmp;
    while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
	    my $id = $1;
	    $id =~ s/_\d+//;
	    $id =~ s/_exon\d+//;
	    $id =~ s/_\d+//;
	    $h{$id}++ unless $tmp{$id};
	    $tmp{$id}++;
	}
    }
    close(IN);
}

foreach my $id (sort {$a cmp $b} keys %h) {
    print $id, "\t", $h{$id}, "\n" if $h{$id} ==  4;
}

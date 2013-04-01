use warnings;
use strict;

my $dir = '/Users/singhal/thesisWork/introgression/';
my @contacts = qw(carlia gillies sjo nBMB);

#full measures
my @measures = qw(divpoly dnds fst);

foreach my $contact (@contacts) {
	my $seq = $dir . 'fullTranscripts/' . $contact . ".fa";
	open(IN, "<$seq");
	my %seq;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*(ENS\S+)/) {
			$seq{$1} = $2;
			}
		}	

	foreach my $measure (@measures) {
	
		my $file = $dir . 'summaryStatistics/' . $contact . '_full_' .  $measure . '.out';
		my $file1 = $dir . 'summaryStatistics/' . $contact . '_full_' .  $measure . '.out2';
		
		open(OUT, ">$file1");
	
		my %tmp;
		open(IN, "<$file");
		my $head = <IN>; chomp($head);
		print OUT "loci\t", $head, "\n";
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
					
			print OUT $seq{$d[0]}, "\t", $line, "\n";
			}
		close(IN);	
		close(OUT);
		
		}					
	}	
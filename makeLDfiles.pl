use warnings;
use strict;

my $dir = '/Users/singhal/thesisWork/introgression/';
my $sdir = $dir . "LD/";
mkdir($sdir) unless(-d $sdir);
my $file = $dir . 'clinesOnly';
my $genome = $dir . 'probeDesign/genome/AnoCar2_70.fa';

my %d; my %loci;
open(IN, "<$file");
while(<IN>) {
	chomp(my $line = $_);
	
	my @d = split(/\t/,$line);
	my $locus = $d[0];
	$locus =~ s/_\S+//;
		
	$d{$d[7]}{$locus}{$d[0]}{$d[1]}{'width'} = $d[3];	
	$d{$d[7]}{$locus}{$d[0]}{$d[1]}{'center'} = $d[2];
	$loci{$d[7]}{$locus}++;		
	}
close(IN);

foreach my $c (keys %d) {
	my $seq = $dir . 'targetSequences/final/' . $c . '_targets.fa.final';
	my $s = parseSeq($seq);
	
	my $tmp = $sdir . "exomeContig_" . $c . ".fa";
	unless(-f $tmp) {
		open(TMP, ">$tmp"); 	
		foreach my $locus (keys %{$d{$c}}) {
			if ($loci{$c}{$locus} > 1) {			
				foreach my $contig (keys %{$d{$c}{$locus}}) {				
					print TMP ">", $contig, "\n", $s->{$contig}, "\n";			
					}						
				}
			}		
		close(TMP);	
		}
	
	my $blatout = $sdir . "blat_" . $c . ".out";
	my $call = system("blat $genome $tmp $blatout -minIdentity=50 -out=blast8") unless (-f $blatout);	
	open(IN, "<$blatout");
	my %match;
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		unless($match{$d[0]}) {
			$match{$d[0]}{'chr'} = $d[1];
			$match{$d[0]}{'start'} = $d[8] - $d[6];
			}
		}
	close(IN);
	
	my $out = $sdir . "LD_" . $c . ".out";
	open(OUT, ">$out");
	foreach my $locus (keys %{$d{$c}}) { 
		foreach my $ec (keys %{$d{$c}{$locus}}) {
			foreach my $pos (keys %{$d{$c}{$locus}{$ec}}) {
				if ($match{$ec}) {	
					print OUT $c, "\t", $locus, "\t", $ec, "\t", $pos, "\t", $d{$c}{$locus}{$ec}{$pos}{'width'}, "\t", $d{$c}{$locus}{$ec}{$pos}{'center'}, "\t"; 				
					my $loc = $match{$ec}{'start'} + $pos;
					print OUT $match{$ec}{'chr'}, "\t", $loc, "\n";
					}				
				else {
					print OUT $c, "\t", $locus, "\t", $ec, "\t", $pos, "\t", $d{$c}{$locus}{$ec}{$pos}{'width'}, "\t", $d{$c}{$locus}{$ec}{$pos}{'center'}, "\t"; 				
					print OUT $ec, "\t", $pos, "\n";
					}
				}
			}
		}					
	close(OUT);	
	}	

sub parseSeq {
	my ($file) = @_;
	open(IN, "<$file");
	my $c = ''; my %seq;	
 	while(<IN>) {
 		chomp(my $line = $_);
 		if ($line =~ m/>(\S+)/) {
 			$c = $1;
 			}
 		else {
 			$seq{$c} .= $line;
 			}
 		}
 	close(IN);
 	return(\%seq);
 	}
use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/activeWork/annotateExome/';
#a TDT file with CONTIGNAME / POSITION / REFERENCE BASE / SNP CALL
my $snpFile = $dir . 'snp.out';
#my annotated exome contig file
my $seq = $dir . 'finalExomeAssembly.fa.annotated';

my %snp;
open(IN, "<$snpFile");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);	
	$snp{$d[0]}{$d[1]}{'1'} = $d[2];
	$snp{$d[0]}{$d[1]}{'2'} = $d[3];
	}
close(IN);

my $s = parseSeq($seq);
	
foreach my $contig (keys %snp) {
	foreach my $pos (keys %{$snp{$contig}}) {
		my $type = 'NA';
		my @gs = @{$s->{$contig}->{'gs'}};
		my @ge = @{$s->{$contig}->{'ge'}};
		
		if ($s->{$contig}->{'utr'}) {
			$type = 'nc';
			}
		elsif (scalar(@gs) > 0) {	
			$type = 'nc';							
			for (my $i = 0; $i < scalar(@gs); $i++) {
				if ($pos >= $gs[$i] && $pos <= $ge[$i]) {
					$type = "coding";
					}
				}
				
			if ($type eq "coding") {			
				my $seq1 = $s->{$contig}->{'seq'}; my $seq2 = $seq1;
				my $s_s1 = substr $seq1, 0, $pos - 1;
				my $e_s1 = substr $seq1, $pos;
				$seq1 = $s_s1 . $snp{$contig}{$pos}{'1'} . $e_s1;
					
				my $s_s2 = substr $seq2, 0, $pos - 1;
				my $e_s2 = substr $seq2, $pos;
				$seq2 = $s_s2 . $snp{$contig}{$pos}{'2'} . $e_s2;					

				my $subseq1;
				for (my $i = 0; $i < scalar(@gs); $i++) {
					my $start = $gs[$i] - 1;
					my $l = $ge[$i] - $gs[$i] + 1;		
					$subseq1 .= substr $seq1, $start, $l;		
					}
				
				my $subseq2;
				for (my $i = 0; $i < scalar(@gs); $i++) {
					my $start = $gs[$i] - 1;
					my $l = $ge[$i] - $gs[$i] + 1;		
					$subseq2 .= substr $seq2, $start, $l;		
					}
					
				my $aa1 = translate($subseq1);
				my $aa2 = translate($subseq2);
			
				if ($aa1 eq $aa2) {
					$type = "syn";
					}
				else {
					$type = "nonsyn";
					}						
				}			
			}
		print $contig, "\t", $pos, "\t", $type, "\n";	
		}	
	}
		
sub parseSeq {
	my ($s) = @_;
	my %s;
	
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			my @gs = $line =~ m/gs(\d+)/g;
			my @ge = $line =~ m/ge(\d+)/g;
			chomp(my $seq = <IN>);
			$s{$id}{'seq'} = $seq;
			$s{$id}{'gs'} = \@gs;
			$s{$id}{'ge'} = \@ge;
			if ($line =~ m/_utr/) {
				$s{$id}{'utr'} = '1';
				}
			}
		}
	return(\%s);	
	}		
	
sub translate {
	my $string = shift;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
#			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}			
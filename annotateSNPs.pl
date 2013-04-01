use warnings;
use strict;

my $dir = '/Users/singhal/thesisWork/introgression/';
my $data = $dir . 'clineAndSummaryStats.out';
my @contacts = qw(carlia gillies nBMB sjo);

my %snp;
open(IN, "<$data");
my $junk = <IN>;
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);	
	$snp{$d[0]}{$d[1]}{$d[2]}{'1'} = 'N';
	$snp{$d[0]}{$d[1]}{$d[2]}{'2'} = 'N';
	}
close(IN);

foreach my $contact (@contacts) {
	my $seq = $dir . "targetSequences/final/" . $contact . "_targets.fa.annotated";
	my $s = parseSeq($seq);
	
	my $vcf1 = $dir . "vcf/" . $contact . ".vcf";
	my $vcf2 = $dir . "vcf/" . $contact . "_HZ.vcf";
	
	my $snp = parseVCF($vcf1,\%snp,$contact);
	$snp = parseVCF($vcf2,\%snp, $contact);
	%snp = %{$snp};
	
	foreach my $contig (keys %{$snp{$contact}}) {
		foreach my $pos (keys %{$snp{$contact}{$contig}}) {
			my $type = 'nc';
			my @gs = @{$s->{$contig}->{'gs'}};
			my @ge = @{$s->{$contig}->{'ge'}};
			
			for (my $i = 0; $i < scalar(@gs); $i++) {
				if ($pos >= $gs[$i] && $pos <= $ge[$i]) {
					$type = "coding";
					}
				}
				
			if ($type eq "coding") {
			
				my $seq1 = $s->{$contig}->{'seq'}; my $seq2 = $seq1;
				my $s_s1 = substr $seq1, 0, $pos - 1;
				my $e_s1 = substr $seq1, $pos;
				$seq1 = $s_s1 . $snp->{$contact}->{$contig}->{$pos}->{'1'} . $e_s1;
					
				my $s_s2 = substr $seq2, 0, $pos - 1;
				my $e_s2 = substr $seq2, $pos;
				$seq2 = $s_s2 . $snp->{$contact}->{$contig}->{$pos}->{'2'} . $e_s2;					

				my $subseq1;
				for (my $i = 0; $i < scalar(@gs); $i++) {
					my $start = $gs[$i] - 1;
					my $l = $ge[$i] - $gs[$i] + 1;		
					$subseq1 .= substr $seq1, $start, $l;		
					}
				if ($s->{$contig}->{'rc'}) {
					$subseq1 = reverse($subseq1);
					$subseq1 =~ tr/ATGC/TACG/;
					}
				
				my $subseq2;
				for (my $i = 0; $i < scalar(@gs); $i++) {
					my $start = $gs[$i] - 1;
					my $l = $ge[$i] - $gs[$i] + 1;		
					$subseq2 .= substr $seq2, $start, $l;		
					}
				if ($s->{$contig}->{'rc'}) {
					$subseq2 = reverse($subseq2);
					$subseq2 =~ tr/ATGC/TACG/;
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
			
			print $contact, "\t", $contig, "\t", $pos, "\t", $type, "\n";
			
			}
		}	
	}
	
sub parseVCF {
	my ($vcf,$snp,$contact) = @_;
	
	my %snp = %{$snp};
	
	open(IN, "<$vcf");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^ENS/) {
			my @d = split(/\t/,$line);
			
			if (exists $snp{$contact}{$d[0]}{$d[1]}) {
				$snp{$contact}{$d[0]}{$d[1]}{'1'} = $d[3];
				$snp{$contact}{$d[0]}{$d[1]}{'2'} = $d[4];
				}
				
			}
		}
	close(IN);	

	return(\%snp);
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
			if ($line =~ m/R_/) {
				$s{$id}{'rc'} = '1';
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
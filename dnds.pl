use strict;
use warnings;

use List::Util qw(min);

my $dir = '/Users/singhal/thesisWork/introgression/';
my @contacts = qw(carlia sjo nBMB gillies);
my $qual = 10;

foreach my $contact (@contacts) {
	my $seq = $dir . 'targetSequences/final/' . $contact . '_targets.fa.annotated';
	
	my $out = $dir . 'summaryStatistics/' . $contact . "_contig_dnds.out";
	open(FINAL, ">$out");
	print FINAL "contig\tdnds\tdn\tds\n";
	
	my $seqInfo = parseSeq($seq);
	my %seq = %{$seqInfo};
	my $var = parseVCF($contact,$dir);
	my %var = %{$var};
		
	foreach my $c (keys %seq) {
		if ($var{$c}) {
			
			#start mutations!
			my $seq1 = $seq{$c}{'seq'}; my $seq2 = $seq1;
			foreach my $pos (keys %{$var{$c}}) {
				if ($var{$c}{$pos} >= $qual) {
					my $s_s1 = substr $seq1, 0, $pos - 1;
					my $e_s1 = substr $seq1, $pos;
					$seq1 = $s_s1 . $var{$c}{$pos}{'1'} . $e_s1;
					
					my $s_s2 = substr $seq2, 0, $pos - 1;
					my $e_s2 = substr $seq2, $pos;
					$seq2 = $s_s2 . $var{$c}{$pos}{'2'} . $e_s2;					
					}
				}	
				
			my @gs = @{$seq{$c}{'gs'}};
			my @ge = @{$seq{$c}{'ge'}};

			my $subseq1;
			for (my $i = 0; $i < scalar(@gs); $i++) {
				my $start = $gs[$i] - 1;
				my $l = $ge[$i] - $gs[$i] + 1;		
				$subseq1 .= substr $seq1, $start, $l;		
				}
			if ($seq{$c}{'rc'}) {
				$subseq1 = reverse($subseq1);
				$subseq1 =~ tr/ATGC/TACG/;
				}
				
			my $subseq2;
			for (my $i = 0; $i < scalar(@gs); $i++) {
				my $start = $gs[$i] - 1;
				my $l = $ge[$i] - $gs[$i] + 1;		
				$subseq2 .= substr $seq2, $start, $l;		
				}
			if ($seq{$c}{'rc'}) {
				$subseq2 = reverse($subseq2);
				$subseq2 =~ tr/ATGC/TACG/;
				}	
				
			#should check lengths here	
			if (length($subseq1) % 3) {
				$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
				$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
				}
			if (length($subseq1) % 3) {
				$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
				$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
				}	
				
			my $aa1 = translate($subseq1);
			my $aa2 = translate($subseq2);
			if ($aa1 =~ m/[A-Z]\*[A-Z]/ || $aa2 =~ m/[A-Z]\*[A-Z]/) {
			
				my @stops = (); my %stops;
				for (my $i = 0; $i < 3; $i++) {
					my $stop = 0;
					my $sseq1 = substr $subseq1, $i;
					my $sseq2 = substr $subseq2, $i;
					my $aa1 = translate($sseq1);
					my $aa2 = translate($sseq2);
					$stop++ if $aa1 =~ m/[A-Z]\*[A-Z]/;
					$stop++ if $aa2 =~ m/[A-Z]\*[A-Z]/;
					my $length1 = length($1) if $aa1 =~ m/([A-Z]+)/;
					my $length2 = length($1) if $aa2 =~ m/([A-Z]+)/;
					push(@stops,$i) if $stop == 0;
					$stops{$i} = min($length1,$length2);
					}
					
				if (scalar(@stops)) {
					if (scalar(@stops) == 1) {
						$subseq1 = substr $subseq1, $stops[0];
						$subseq2 = substr $subseq2, $stops[0];
						
						if (length($subseq1) % 3) {
							$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
							$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
							}
						if (length($subseq1) % 3) {
							$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
							$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
							}	
						
						my ($omega,$dn,$ds) = runPaml_pair($subseq1,$subseq2);
						print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
						}
					#what to do if there is more than one high hit?
					else {
						$subseq1 = substr $subseq1, $stops[0];
						$subseq2 = substr $subseq2, $stops[0];
						
						if (length($subseq1) % 3) {
							$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
							$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
							}
						if (length($subseq1) % 3) {
							$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
							$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
							}
							
						my ($omega,$dn,$ds) = runPaml_pair($subseq1,$subseq2);
						print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
						}
					}
				else {
					my @s = sort {$stops{$b} <=> $stops{$a}} keys %stops;
					#my frame is $s[0]
					$subseq1 = substr $subseq1, $s[0];
					$subseq2 = substr $subseq2, $s[0];
					
					if (length($subseq1) % 3) {
							$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
							$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
							}
					if (length($subseq1) % 3) {
						$subseq1 = substr $subseq1, 0, length($subseq1) - 1;
						$subseq2 = substr $subseq2, 0, length($subseq2) - 1;
						}
							
					my ($omega,$dn,$ds) = runPaml_pair($subseq1,$subseq2);
					print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
					}					
				}
			else {
				my ($omega,$dn,$ds) = runPaml_pair($subseq1,$subseq2);
				print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
				}							
			}
		else {	
			print FINAL "$c\t0\t0\t0\n";
			}
		}		
	}
	
sub runPaml_pair {
	my ($seq1,$seq2) = @_;
	my $out = "/Users/singhal/programs/codingseq.phy";
				
	my $length = length($seq1);			
	open(OUT, ">$out");
	print OUT " 2 $length\n";	
	print OUT "seq1  $seq1\nseq2  $seq2";
	close(OUT);
					
	my $call = system("yn00 /Users/singhal/programs/yn00.ctl");
	open(IN, "<yn");
	my $omega = '0'; my $dn = 0; my $ds = 0;
	while(<IN>) {
		if ($_ =~ m/omega/) {
			my $line = <IN>;
			chomp(my $data = <IN>);
			my @d = ($data =~ /(\S+)/g);
			$omega = $d[6];
			$omega = 0 if $omega == 99;		
			$dn = $d[7];
			$dn = 0 if $dn eq '-0.0000';
			$ds = $d[10];
			$ds = 0 if $ds eq '-0.0000';
			}
		}
	close(IN);	
	
	unlink($out);
	
	return($omega,$dn,$ds);	
	}	

sub parseVCF {
	my ($contact,$dir) = @_;
	
	my $vcf1 = $dir . 'vcf/' . $contact . ".vcf";
	my $vcf2 = $dir . 'vcf/' . $contact . "_HZ.vcf";
	
	my %var;
	
	open(IN, "<$vcf1");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^ENS/) {
			my @d = split(/\t/,$line);
			
			my $count1 = 0; my $count2 = 0;
			for (my $i = 9; $i < 14; $i++) {
				$count1++ if $d[$i] =~  m/1\//;
				$count1++ if $d[$i] =~  m/\/1/;
				}
			for (my $i = 14; $i < 19; $i++) {
				$count2++ if $d[$i] =~  m/1\//;
				$count2++ if $d[$i] =~  m/\/1/;	
				}
								
			if ($count1 >= 5) {
				$var{$d[0]}{$d[1]}{'1'} = $d[4];
				}
			else {
				$var{$d[0]}{$d[1]}{'1'} = $d[3];
				}
			
			if ($count2 >= 5) {
				$var{$d[0]}{$d[1]}{'2'} = $d[4];
				}
			else {
				$var{$d[0]}{$d[1]}{'2'} = $d[3];
				}
				
			$var{$d[0]}{$d[1]}{'qual'} = $d[5];
			}
		}
	close(IN);	
	
	open(IN, "<$vcf2");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^ENS/) {
			my @d = split(/\t/,$line);
			if (exists $var{$d[0]}{$d[1]}) {
				$var{$d[0]}{$d[1]}{'qual'} = $d[5] unless $var{$d[0]}{$d[1]}{'qual'} > $d[5];
				}
			}
		}
	close(IN);
		
	return(\%var);
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
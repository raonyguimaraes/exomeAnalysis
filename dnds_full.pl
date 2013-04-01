use strict;
use warnings;

use List::Util qw(min);

my $dir = '/Users/singhal/thesisWork/introgression/';
my @contacts = qw(carlia sjo nBMB gillies);
my $qual = 10;

foreach my $contact (@contacts) {
	my $seq = $dir . 'fullTranscripts/' . $contact . '.fa';
	
	my $out = $dir . 'summaryStatistics/' . $contact . "_full_dnds.out";
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
				
			$seq1 = substr $seq1, $seq{$c}{'start'}, $seq{$c}{'length'};
			$seq2 = substr $seq2, $seq{$c}{'start'}, $seq{$c}{'length'};
				
			my $aa1 = translate($seq1);
			my $aa2 = translate($seq2);
			if ($aa1 =~ m/[A-Z]\*[A-Z]/ || $aa2 =~ m/[A-Z]\*[A-Z]/) {
			
				my @stops = (); my %stops;
				for (my $i = 0; $i < 3; $i++) {
					my $stop = 0;
					my $subseq1 = substr $seq1, $i;
					my $subseq2 = substr $seq2, $i;
					my $aa1 = translate($subseq1);
					my $aa2 = translate($subseq2);
					$stop++ if $aa1 =~ m/[A-Z]\*[A-Z]/;
					$stop++ if $aa2 =~ m/[A-Z]\*[A-Z]/;
					my $length1 = length($1) if $aa1 =~ m/([A-Z]+)/;
					my $length2 = length($1) if $aa2 =~ m/([A-Z]+)/;
					push(@stops,$i) if $stop == 0;
					$stops{$i} = min($length1,$length2);
					}
					
				if (scalar(@stops)) {
					if (scalar(@stops) == 1) {
						$seq1 = substr $seq1, $stops[0];
						$seq2 = substr $seq2, $stops[0];
						
						if (length($seq1) % 3) {
							$seq1 = substr $seq1, 0, length($seq1) - 1;
							$seq2 = substr $seq2, 0, length($seq2) - 1;
							}
						if (length($seq1) % 3) {
							$seq1 = substr $seq1, 0, length($seq1) - 1;
							$seq2 = substr $seq2, 0, length($seq2) - 1;
							}
						
						my ($omega,$dn,$ds) = runPaml_pair($seq1,$seq2);
						print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
						}
					#what to do if there is more than one high hit?
					else {
						$seq1 = substr $seq1, $stops[0];
						$seq2 = substr $seq2, $stops[0];
						
						if (length($seq1) % 3) {
							$seq1 = substr $seq1, 0, length($seq1) - 1;
							$seq2 = substr $seq2, 0, length($seq2) - 1;
							}
						if (length($seq1) % 3) {
							$seq1 = substr $seq1, 0, length($seq1) - 1;
							$seq2 = substr $seq2, 0, length($seq2) - 1;
							}
							
						my ($omega,$dn,$ds) = runPaml_pair($seq1,$seq2);
						print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
						}
					}
				else {
					my @s = sort {$stops{$b} <=> $stops{$a}} keys %stops;
					#my frame is $s[0]
					$seq1 = substr $seq1, $s[0];
					$seq2 = substr $seq2, $s[0];
					
					if (length($seq1) % 3) {
						$seq1 = substr $seq1, 0, length($seq1) - 1;
						$seq2 = substr $seq2, 0, length($seq2) - 1;
						}
					if (length($seq1) % 3) {
						$seq1 = substr $seq1, 0, length($seq1) - 1;
						$seq2 = substr $seq2, 0, length($seq2) - 1;
						}
							
					my ($omega,$dn,$ds) = runPaml_pair($seq1,$seq2);
					print FINAL $c, "\t", $omega, "\t", $dn, "\t", $ds, "\n";
					}					
				}
			else {
				my ($omega,$dn,$ds) = runPaml_pair($seq1,$seq2);
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
	
	my $vcf1 = $dir . 'vcf/' . $contact . "_full.vcf";
	
	my %var;
	
	open(IN, "<$vcf1");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^contig/) {
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
		
	return(\%var);
	}		
	
sub parseSeq {
	my ($s) = @_;
	my %s;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			my $start = $1 if $line =~ m/gs(\d+)/;
			my $end = $1 if $line =~ m/ge(\d+)/;			
			chomp(my $seq = <IN>);
			
			$s{$c}{'seq'} = $seq;
			$s{$c}{'start'} = $start - 1;
			my $length = $end - $start + 1;
			unless( $length % 3) {
				$length = 3 * int($length / 3);
				}
			$s{$c}{'length'} = $length;
			}
		}
	close(IN);	
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
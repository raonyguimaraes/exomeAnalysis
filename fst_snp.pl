use warnings;
use strict;

my $dir = '/Users/singhal/thesisWork/introgression/';
my @contacts = qw(carlia sjo nBMB gillies);
my $qual = 10;

foreach my $contact (@contacts) {
	my $seq = $dir . 'targetSequences/final/' . $contact . '_targets.fa.final';
	
	my $seqInfo = parseSeq($seq);
	my %seq = %{$seqInfo};
	my $var = parseVCF($contact,$dir);
	my %var = %{$var};
		
	my $out = $contact . '_SNP_fst.out';
	open(OUT, ">$out");
	print OUT "gene\tpos\tFst\n";
			
	foreach my $c (keys %var) {
		my @fst;
		foreach my $pos (keys %{$var{$c}}) {
			my ($a1, $a2);
			if ($var{$c}{$pos}{'1'}) {
				$a1 = $var{$c}{$pos}{'1'};
				}
			else {
				$a1 = 0;
				}
			if ($var{$c}{$pos}{'2'}) {
				$a2 = $var{$c}{$pos}{'2'};
				}
			else {
				$a2 = 0;
				}
			
			my $fis1 = 2 * ( $a1 / 10 ) * ( (10 - $a1) / 10 );	
			my $fis2 = 2 * ( $a2 / 10 ) * ( (10 - $a2) / 10 );	
			my $fis = ($fis1 + $fis2) / 2;
			my $a = $a1 + $a2;
			my $fit  = 2 * ( $a / 20 ) * ( (20 - $a) / 20 );
			my $fst;
			if ($fit == 0) {
				$fst = 0;
				}
			else {	
				$fst = ($fit - $fis) / $fit;
				}	
			print OUT $c, "\t", $pos, "\t", $fst, "\n";
			}		
		}		
	close(OUT);	
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
			for (my $i = 9; $i < 14; $i++) {
				$var{$d[0]}{$d[1]}{'1'}++ if $d[$i] =~  m/1\//;
				$var{$d[0]}{$d[1]}{'1'}++ if $d[$i] =~  m/\/1/;
				}
			for (my $i = 14; $i < 19; $i++) {
				$var{$d[0]}{$d[1]}{'2'}++ if $d[$i] =~  m/1\//;
				$var{$d[0]}{$d[1]}{'2'}++ if $d[$i] =~  m/\/1/;	
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
			my $c = $1;
			chomp(my $seq = <IN>);
			$s{$c} = $seq;
			}
		}
	close(IN);	
	return(\%s);	
	}	
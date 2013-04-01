use warnings;
use strict;

my @contacts = qw(sjo gillies nBMB carlia);
my $dir = '/Users/singhal/thesisWork/introgression/';
my $qual = 10;

foreach my $contact (@contacts) {
	my $seq = $dir . 'fullTranscripts/' . $contact . ".fa";
	my $s = parseSeq($seq);
	my $var = parseVCF($contact,$dir);
	makeHaplo($contact,$s, $var);	
	}
	
sub makeHaplo {
	my ($contact,$s,$var) = @_;
	my %var = %{$var};
	my %seq = %{$s};
	
	open(OUT, ">$contact" . "_full_divpoly.out");
	print OUT "contig\tpoly1\tpoly2\traw\tnet\tlength\n";
	
	foreach my $c (keys %var) {
		my @haplo;
		for (my $i = 9; $i < 19; $i++) {
			for (my $j = 0; $j < 2; $j++) {
				my $haplo;
				foreach my $pos (keys %{$var{$c}}) {
					if ($var{$c}{$pos}{'qual'} >= $qual) {
						my $rand = rand();
						if ($rand > 0.5) {
							$haplo .= $var{$c}{$pos}{$i}[1];
							}
						else {
							$haplo .= $var{$c}{$pos}{$i}[0];
							}
						}	
					}
				push(@haplo,$haplo);
				}
			}
			
		if ($haplo[0]) {	
			my @h1 = @haplo[0..9];	
			my @h2 = @haplo[10..19];
			#now time to calculate polymorphism and divergence
			my $within1 =  sprintf("%.4f",varwithin(\@h1)/$seq{$c});
			my $within2 = sprintf("%.4f",varwithin(\@h2)/$seq{$c});
			my $raw = sprintf("%.4f",varbtn(\@h1,\@h2)/$seq{$c});
			my $net = $raw - ($within1 + $within2) / 2;
		
			print OUT $c, "\t", $within1, "\t", $within2, "\t", $raw, "\t", $net, "\t$seq{$c}\n";
			delete($seq{$c});
			}
		else {
			print OUT $c, "\t", "0", "\t", "0", "\t", "0", "\t", "0", "\t$seq{$c}\n";
			delete($seq{$c});
			}
		}
		
	foreach my $c (keys %seq) {
		print OUT $c, "\t", "0", "\t", "0", "\t", "0", "\t", "0", "\t$seq{$c}\n";
		}
		
	close(OUT);	
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
			for (my $i = 9; $i < 19; $i++) {
				my $geno1 = $1 if $d[$i] =~ m/([0|1])\//;
				my $geno2 = $1 if $d[$i] =~ m/\/([0|1])/;
				$var{$d[0]}{$d[1]}{$i} = [$geno1,$geno2];
				$var{$d[0]}{$d[1]}{'qual'} = $d[5];
				}
			}
		}
	close(IN);	
		
	return(\%var);
	}
	
	
sub varwithin {
	my ($haplo) = @_;
	my @h = @{$haplo};
	my $sum; my $compare;
	for (my $i = 0; $i  < scalar(@h); $i++) {
		for (my $j = 0; $j  < scalar(@h); $j++) {
			my $diff = 0;
			my $mask = $h[$i] ^ $h[$j];
			$diff++ while ($mask =~ /[^\0]/g);
			$sum += $diff;
			$compare++;
			}
		}
	$sum = $sum/$compare;	
	return($sum);
	}


sub varbtn {
	my ($haplo1,$haplo2) = @_;
	my @h1 = @{$haplo1};
	my @h2 = @{$haplo2};
	my $sum; my $compare;
	for (my $i = 0; $i  < scalar(@h1); $i++) {
		for (my $j = 0; $j  < scalar(@h2); $j++) {
			my $diff = 0;
			my $mask = $h1[$i] ^ $h2[$j];
			$diff++ while ($mask =~ /[^\0]/g);
			$sum += $diff;
			$compare++;
			}
		}
	$sum = $sum/$compare;	
	return($sum);
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
			$s{$c} = length($seq);
			}
		}
	close(IN);	
	return(\%s);	
	}		
	
	

	
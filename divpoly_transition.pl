use warnings;
use strict;

my @contacts = qw(sjo gillies nBMB carlia);
my $dir = '/Users/singhal/thesisWork/introgression/';

foreach my $contact (@contacts) {
	my $seq = $dir . 'targetSequences/final/' . $contact . "_targets.fa.final";
	my $s = parseSeq($seq);
	my $var = parseVCF($contact,$dir);
	makeHaplo($contact,$s,$var);	
	}
	
sub makeHaplo {
	my ($contact,$s,$var) = @_;
	my %var = %{$var};
	my %seq = %{$s};
	
	open(OUT, ">$contact" . "_transition_divpoly.out");
	print OUT "contig\tpoly1\tpoly2\traw\tnet\tlength\n";
	
	foreach my $c (keys %var) {
		my @h1; my @h2;
	
		for (my $i = 0; $i < 10; $i++) {
			my $haplo;
			foreach my $pos (keys %{$var{$c}}) {
				if (exists $var{$c}{$pos}{'10kN'} && exists $var{$c}{$pos}{'10kS'}) {
					my $rand = rand();
					if ($rand <= $var{$c}{$pos}{'10kN'}) {
						$haplo .= '0';
						}
					else {
						$haplo .= '1';
						}
					}	
				}
			push(@h1,$haplo);	
			}	
			
		for (my $i = 0; $i < 10; $i++) {
			my $haplo;
			foreach my $pos (keys %{$var{$c}}) {
				if (exists $var{$c}{$pos}{'10kN'} && exists $var{$c}{$pos}{'10kS'}) {
					my $rand = rand();
					if ($rand <= $var{$c}{$pos}{'10kS'}) {
						$haplo .= '0';
						}
					else {
						$haplo .= '1';
						}
					}	
				}
			push(@h2,$haplo);	
			}	
				
		if ($h1[0]) {	
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
	
	my $vcf1 = $dir . 'clineAF/' . $contact . ".cline.out";
	
	my %var;
	
	open(IN, "<$vcf1");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($d[2] =~ m/10k/) {
			$var{$d[0]}{$d[1]}{$d[2]} = $d[3] unless ($d[3] eq 'NA');
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
	
	

	
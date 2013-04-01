use warnings;
use strict;

my $dir = '/Users/singhal/thesisWork/introgression/';
my $data =  $dir . 'clineData';
my @contacts = qw(carlia gillies sjo nBMB);

my %d;
my %head;

my @index = qw(2 3 4 5 6);
open(IN, "<$data");
my $head = <IN>; chomp($head);
my @head = split(/\t/,$head);
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	
	for (my $i = 0; $i < scalar(@index); $i++) {
		$d{$d[7]}{$d[0]}{$d[1]}{$head[$index[$i]]} = $d[$index[$i]];
		$head{$head[$index[$i]]}++;
		}	
	}
close(IN);

#snp measures
foreach my $contact (@contacts) {
	my $snp = $dir . 'summaryStatistics/' . $contact . "_SNP_fst.out";
	
	open(IN, "<$snp");
	my $junk = <IN>;
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		
		if (exists $d{$contact}{$d[0]}{$d[1]}) {
			$d{$contact}{$d[0]}{$d[1]}{'Fst_snp'} = $d[2];
			}
		}
	close(IN);	
	}	
$head{'Fst_snp'}++;	

#contig measures
my @measures = qw(divpoly dnds fst);
foreach my $contact (@contacts) {
	foreach my $measure (@measures) {
	
		my $file = $dir . 'summaryStatistics/' . $contact . '_contig_' .  $measure . '.out';
	
		my %tmp;
		open(IN, "<$file");
		my $head = <IN>; chomp($head);
		my @head = split(/\t/,$head);
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
					
			for (my $i = 1; $i < scalar(@d); $i++) {
				$tmp{$d[0]}{'contig_' . $head[$i]} = $d[$i];
				}
			}
		close(IN);	
		
		foreach my $c (keys %{$d{$contact}}) {
			if ($tmp{$c}) {
				foreach my $key (keys %{$tmp{$c}}) {
					foreach my $pos (keys %{$d{$contact}{$c}}) {
						$d{$contact}{$c}{$pos}{$key} = $tmp{$c}{$key};
						$head{$key}++;
						}
					}
				}
			}
		}					
	}	

#transition measures
foreach my $contact (@contacts) {	
	my $file = $dir . 'summaryStatistics/' . $contact . '_transition_divpoly.out';
	
	my %tmp;
	open(IN, "<$file");
	my $head = <IN>; chomp($head);
	my @head = split(/\t/,$head);
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
					
		for (my $i = 1; $i < scalar(@d); $i++) {
			$tmp{$d[0]}{'trans_' . $head[$i]} = $d[$i];
			}
		}
	close(IN);	
		
	foreach my $c (keys %{$d{$contact}}) {
		if ($tmp{$c}) {
			foreach my $key (keys %{$tmp{$c}}) {
				foreach my $pos (keys %{$d{$contact}{$c}}) {
					$d{$contact}{$c}{$pos}{$key} = $tmp{$c}{$key};
					$head{$key}++;
					}
				}
			}
		}					
	}

#full measures
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
	
		my %tmp;
		open(IN, "<$file");
		my $head = <IN>; chomp($head);
		my @head = split(/\t/,$head);
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
					
			for (my $i = 1; $i < scalar(@d); $i++) {
				$tmp{$seq{$d[0]}}{'full_' . $head[$i]} = $d[$i];
				}
			}
		close(IN);	
		
		foreach my $c (keys %{$d{$contact}}) {
			my $new_c = $c;
			$new_c =~ s/_\d+//;
			$new_c =~ s/_exon\d+//;
			
			if ($tmp{$new_c}) {
				foreach my $key (keys %{$tmp{$new_c}}) {
					foreach my $pos (keys %{$d{$contact}{$c}}) {
						$d{$contact}{$c}{$pos}{$key} = $tmp{$new_c}{$key};
						$head{$key}++;
						}
					}
				}
			}
		}					
	}	

#do outlier business
foreach my $contact (@contacts) {
	my @width;
	foreach my $contig (keys %{$d{$contact}}) {
		foreach my $pos (keys %{$d{$contact}{$contig}}) {
			if ($d{$contact}{$contig}{$pos}{'width'} ne 'NA') {
				push(@width,$d{$contact}{$contig}{$pos}{'width'}) if $d{$contact}{$contig}{$pos}{'width'} > 0;
				}
			}
		}
		
	@width = sort {$a <=> $b} @width;	

	my $wmin = $width[int(scalar(@width) * 0.05)];
	my $dump = $width[int(scalar(@width) * 0.95)];
	my $wmax = $width[int(scalar(@width) * 0.90)];
	
	foreach my $contig (keys %{$d{$contact}}) {
		foreach my $pos (keys %{$d{$contact}{$contig}}) {
			my $width = $d{$contact}{$contig}{$pos}{'width'};
			if ($width ne 'NA') {
				if ($width < 0) {
					#do nothing
					}
				else {
					if ($width <= $wmin) {
						 $d{$contact}{$contig}{$pos}{'outlier'} = 'outlier';
						 $d{$contact}{$contig}{$pos}{'outtype'} = 'narrow';
						}
					else {
						if ($width >= $wmax) {
							if ($width >= $dump) {
								}
							else {
								$d{$contact}{$contig}{$pos}{'outlier'} = 'outlier';
								$d{$contact}{$contig}{$pos}{'outtype'} = 'wide';
								}						
							}
						else {
							$d{$contact}{$contig}{$pos}{'outlier'} = 'normal';
							$d{$contact}{$contig}{$pos}{'outtype'} = 'normal';
							}
						}
					}
				}
			}
		}
	}
$head{'outlier'}++;
$head{'outtype'}++;

#print this nonsense
print "contact\tcontig\tpos\t";
print join("\t", sort {$a cmp $b} keys %head);
print "\n";

foreach my $contact (@contacts) {
	foreach my $contig (keys %{$d{$contact}}) {
		foreach my $pos (sort {$a <=> $b} keys %{$d{$contact}{$contig}}) {
			print $contact, "\t", $contig, "\t", $pos, "\t";
			foreach my $key (sort {$a cmp $b} keys %head) {
				if ($d{$contact}{$contig}{$pos}{$key}) {
					print $d{$contact}{$contig}{$pos}{$key}, "\t";
					}
				else {
					print "NA\t";
					}
				}
			print "\n";
			}
		}
	}
	

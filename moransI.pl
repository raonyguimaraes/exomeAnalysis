use warnings;
use strict;

my @contact = qw(sjo nBMB gillies carlia);
#my @contact = qw(carlia);
my @dist = qw(0 100 200 300 400 500 600 700 800 900 1000);

foreach my $contact (@contact) {
	my $file = "/Users/singhal/thesisWork/introgression/LD/LD_" . $contact . ".out";

	my %d;
	open(IN,"<$file");
	while(<IN>) {
		chomp(my $line = $_);
		
		my @d = split(/\t/,$line);
		push(@{$d{$d[2]}},\@d); 
	
		}
	close(IN);	

	for (my $n = 0; $n < scalar(@dist) - 1; $n++) {
	
		my $nLoci = 0;
		my $sumNumerator = 0;
		my $nComparisons = 0;
		my $sumDenominator = 0;	
	
		my $dist1 = $dist[$n];
		my $dist2 = $dist[$n + 1];
	
		foreach my $c (keys %d) {
			my @tmp = @{$d{$c}};
			if (scalar(@tmp) > 1) {
				$nLoci += scalar(@tmp);  
			
				for (my $i = 0; $i < scalar(@tmp); $i++) {
					$sumDenominator += $tmp[$i][4] ** 2;
					}
			
				for (my $j = 0; $j < scalar(@tmp); $j++) {
					for (my $k = $j + 1; $k < scalar(@tmp); $k++) {
						if ( abs($tmp[$j][3] - $tmp[$k][3]) >= $dist1) {
							if ( abs($tmp[$j][3] - $tmp[$k][3]) <= $dist2) {
								$nComparisons++;
								$sumNumerator += ($tmp[$j][4] * $tmp[$k][4]);
								}
							}	
						}			
					}	
				}
			}	
		my $moranI = ($nLoci * $sumNumerator) / ($nComparisons * $sumDenominator);
		print $contact, "\t", $dist1, "\t", $moranI, "\t$nComparisons\n";
		}
	}	
	
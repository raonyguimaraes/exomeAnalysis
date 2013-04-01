use warnings;
use strict;

my @names = qw(carlia);
my @pop = qw(10kN 2.5kN 1kN northTail center southTail 1kS 2.5kS 10kS);

foreach my $name (@names) {
	my $geno = '/Users/singhal/thesisWork/introgression/benchtruth/geno/' . $name . 'Genotypes';
	my $out = '/Users/singhal/thesisWork/introgression/benchtruth/af/' . $name . 'AF';
	my $pop = '/Users/singhal/thesisWork/introgression/pooledLibraries/pooledSamples/' . $name . 'Amounts.txt';
	
	my %pop;
	open(IN, "<$pop");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$pop{$d[1]}{$d[0]}++ unless $d[0] eq 'water';
		}
	close(IN);
	
	my %snp;
	open(IN, "<$geno");
	chomp(my $head = <IN>); my @head = split(/\t/,$head);
	my @loci = @head[1..$#head];
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		for (my $i = 1; $i < scalar(@head); $i++) {
			$snp{$d[0]}{$head[$i]} = $d[$i];
			}
		}
	close(IN);	
	
	open(OUT, ">$out");
	foreach my $loci (@loci) {
		foreach my $pop (@pop) {
			my $popsize = 0;
			my $actualP = scalar(keys %{$pop{$pop}});
			my $af = 0;
			foreach my $i (keys %{$pop{$pop}}) {
				my $geno = $snp{$i}{$loci};
				if ($geno) {	
					if ($geno eq 'N') {
						$popsize++;
						}
					elsif ($geno eq 'H') {
						$af++;
						$popsize++;
						}
					elsif ($geno eq 'C' || $geno eq 'S') {
						$af += 2;
						$popsize++;
						}	
					}	
				else {
					#print $loci, "\t", $pop, "\t", $i, "\n";
					}
				}
			my $alleleFreq;	
			if ($popsize) {	
				$alleleFreq = sprintf("%.03f",$af / (2 * $popsize));	
				}
			else {
				$alleleFreq = 'NA';
				}
			print OUT $name, "\t", $loci, "\t", $pop, "\t", $popsize, "\t", $actualP, "\t", $alleleFreq, "\n";   	
			#print $name, "\t", $loci, "\t", $pop, "\t", $popsize, "\t", $actualP, "\t", $alleleFreq, "\n" if $actualP ne $popsize; 	
			}
		}
	close(OUT);	
	}	

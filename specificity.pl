use warnings;
use strict;

my @contact = qw(gillies nBMB sjo carlia);
my @pops = qw(10kN 2kN 1kN nTail center sTail 1kS 2kS 10kS);

my $i = 0; my $j = 0;

my $file = shift;

open(IN, "<$file");
while(<IN>) {
	chomp(my $line = $_);
			if ($line =~ m/reads/) {
				my $reads1 = $1 if $line =~ m/^(\d+)/;
				
				chomp(my $junk1 = <IN>);
				chomp(my $junk2 = <IN>);
				chomp(my $junk3 = <IN>);
				chomp(my $junk4 = <IN>);
				
				chomp($line = <IN>);
				my $percent1 = $1 if $line =~ m/^(\S+)\%/;
				$percent1 = $percent1 / 100;
				
				chomp($line = <IN>);
				my $reads2 = $1 if $line =~ m/^(\d+)/;
				
				chomp(my $junk5 = <IN>);
				chomp(my $junk6 = <IN>);
				chomp(my $junk7 = <IN>);
				chomp(my $junk8 = <IN>);
				
				chomp($line = <IN>);
				my $percent2 = $1 if $line =~ m/^(\S+)\%/;
				$percent2 = $percent2 / 100;
				
				my $per =  $reads1 * $percent1 + $reads2 * $percent2;
				$per = $per / $reads1;
				
				print $contact[$i], "\t", $pops[$j], "\t", $per, "\n";
				
				$j++;
				if ($j % 9 == 0) {
					$j = 0; $i++;
					}
				}
	}	
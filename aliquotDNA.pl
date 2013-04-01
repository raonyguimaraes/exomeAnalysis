use warnings;
use strict;

my $target = 5000;
my $totVol = ($target/10) * 0.9;
my $intro = shift;

my %count;
my %all;

open(IN, "<$intro");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	$count{$d[2]}++;
	$all{$d[2]}{$d[0]} = $d[4];
	}
close(IN);

foreach my $pop (sort {$a cmp $b} keys %count) {
	my $vol = $totVol;
	my $each = $target/$count{$pop};
	foreach my $indiv (sort {$a cmp $b} keys %{$all{$pop}}) {
		my $ali = sprintf("%.1f",$each / $all{$pop}{$indiv});
		$vol = $vol - $ali;
		print $indiv, "\t", $pop, "\t", $all{$pop}{$indiv}, "\t", $ali, "\n";
		}
	print "water", "\t", $pop, "\t", "NA", "\t", $vol, "\n";
	}

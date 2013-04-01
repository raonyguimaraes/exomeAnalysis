use warnings;
use strict;

my $file = shift;
my %total;
my %min;

my $vol = 75;

open(IN, "<$file");
my $header = <IN>;
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	my $contact = $1 if $d[1] =~ m/(\S+)_/;
	$total{$contact}{$d[0]}{'qbit'} = $d[3];
	$total{$contact}{$d[0]}{'nano'} = $d[2];
	if ($min{$contact}) {
		$min{$contact} = $d[3] if $d[3] < $min{$contact};
		}
	else {
		$min{$contact} = $d[3];
		}
	}
close(IN);

foreach my $c (keys %total) {
	my $target = $min{$c} * $vol;
	my $nanoNG; my $totVol;
	foreach my $id (sort {$a cmp $b} keys %{$total{$c}}) {
		my $idVol = sprintf("%.1f",$target/$total{$c}{$id}{'qbit'});
		my $nano = $idVol * $total{$c}{$id}{'nano'};
		$nanoNG += $nano;
		$totVol += $idVol;
		print $c, "\t", $id, "\t", $idVol, "\n";
		}
	print $c, "\t", $nanoNG, "\t", $totVol, "\n";
	}	
	
	
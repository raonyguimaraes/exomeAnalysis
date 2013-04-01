use strict;
use warnings;

my $dir = '/Users/singhal/thesisWork/introgression/probeDesign/probeDesign/';
my @cz = qw(Lampro1 Lampro2 Sapro Carlia);
#my @cz = qw(Lampro2);

foreach my $cz (@cz) {
	my $file = $dir . $cz . '_probes.txt';
	open(IN,"<$file");
	my %d;
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		my $id = $d[0];
		$id =~ s/chr\d+_//;
		$id =~ s/^L2_//;
		$id =~ s/\:.*//;
		$id =~ s/_[0|1]$//;
		$id =~ s/_\d+//;
		$d{$id}++;
		}
	close(IN);	
		
	my %seq;
	my $seq = $dir . $cz . '/' . $cz . '_preMasked.fa';
	open(SEQ, "<$seq");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/>chr\d+_(\S+)_[0|1]$/) {
			my $id = $1;
			$id =~ s/_\d+//;
			chomp(my $s = <SEQ>);
			if ($d{$id}) {
				push(@{$seq{$id}}, $s);
				}
			}
		}	
	close(SEQ);		

	my $error = $cz . '_error.out';
	my $out = $cz . '_targets.fa';
	open(OUT, ">$out");
	open(ERROR, ">$error");
		
	foreach my $id (keys %seq) {
		open(TMP, ">tmp.fa");
		my $tracker = 1;
		foreach my $s (@{$seq{$id}}) {
			print TMP ">", $id, "_", $tracker, "\n", $s, "\n";
			$tracker++;
			}
		close(TMP);
		my $call = system("cap3 tmp.fa -z 1 -o 16 -e 11");
		my $assembled = 'tmp.fa.cap.contigs';
		my $singlets =  'tmp.fa.cap.singlets';
		
		my %c;
		my $c = readFile($assembled,\%c,'contig');
		$c = readFile($singlets,$c,'single');
		%c = %{$c};
		
		if (scalar(keys %c) == 1) {
			foreach my $i (keys %c) {
				print OUT ">", $id, "\n";
				print OUT $c{$i}, "\n";
				}
			}
		elsif (scalar(keys %c) > 1) {
			print ERROR $id, "\tTooMany\n";
			my $j = 1;
			foreach my $i (keys %c) {
				print OUT ">", $id, "_", $j, "\n";
				$j++;
				print OUT $c{$i}, "\n";
				}
			}
		else {
			print ERROR $id, "\tTooFew\n";
			}	
		}		
	close(OUT); close(ERROR);	
	}
	
sub readFile {
	my ($file,$hash,$base) = @_;
	if (-f $file) {
		open(IN, "<$file");	
		my $id; my $tracker = 1;
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $base . $tracker;
				$tracker++;
				}
			else {
				$hash->{$id} .= $line;
				}
			}
		close(IN);	
		}	
	return($hash);
	}
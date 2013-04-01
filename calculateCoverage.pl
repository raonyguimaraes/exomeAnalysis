use warnings;
use strict;

my @contacts = qw(carlia sjo gillies nBMB);
my @pops = qw(10kN 2kN 1kN nTail center sTail 1kS 2kS 10kS);
my $dir = '/media/Elements/introgression/ancMapping/';
my $seqdir = '/media/Elements/introgression/targetSequences/final/';

makeDepth();
parseDepth();

sub parseDepth {
    foreach my $c (@contacts) {
	
	my $seq = $seqdir . $c . "_targets.fa.final";
	my %seq;
	open(SEQ, "<$seq");
	while(<SEQ>) {
	    chomp(my $line = $_);
	    if ($line =~ m/>(\S+)/) {
		my $id = $1;
		chomp(my $s = <SEQ>);
		$seq{$id} = length($s);
	    }
	}
	close(SEQ);

	my %depth;
	my $depth = $dir . $c . ".depth.out";
	open(IN, "<$depth");
	while(<IN>) {
	    chomp(my $line = $_);
	    my @d = split(/\t/,$line);
	    for (my $i = 2; $i < scalar(@d); $i++) {
		$depth{$d[0]}{$pops[$i - 2]} += $d[$i];
	    }
	}
	close(IN);

	my $out = $dir . $c . ".depth.summarized";
	open(OUT,">$out");
	foreach my $s (keys %seq) {
	    print OUT $s, "\t";
	    if ($depth{$s}) {
		for (my $i = 0; $i < scalar(@pops); $i++) {
		    my $cov = sprintf("%.1f",$depth{$s}{$pops[$i]}/$seq{$s});
		    print OUT $cov, "\t";
		}
	    }
	    else {
		for (my $i = 0; $i < scalar(@pops); $i++) {
		    print OUT "NA\t";
		}
	    }
	    print OUT "\n";
	}
	close(OUT);
    }
}
	
sub makeDepth {
    foreach my $c (@contacts) {   
	my @files;
	foreach my $pop (@pops) {
	    my $file = $dir . $c . '/' . $pop . ".sorted.bam";
	    push(@files,$file);
	}

	my $files = join("\t",@files);
	my $out = $dir . $c . ".depth.out";
	my $call = system("samtools depth $files > $out") unless (-f $out);
    }
}

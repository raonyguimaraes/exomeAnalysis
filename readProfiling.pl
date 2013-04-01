use strict;
use warnings;

my $dir = '/media/Elements/introgression/';
my @c = qw(carlia sjo nBMB gillies);
my @pop = qw(10kN 10kS 2kN 2kS 1kN 1kS nTail sTail center);
#my @c = qw(gillies);
#my @pop = qw(test);

my %data;

foreach my $c (@c) {
    foreach my $pop (@pop) {

	my $file1 = $dir . $c . '/' . $pop . '/' . $pop . '_1.fastq.gz';
	my $file2 = $file1; $file2 =~ s/_1/_2/;
	my $count1 = count($file1);
	my $count2 = count($file2);
	my $totraw = ($count1 + $count2 )/4 * 100; 
	$data{$c}{$pop}{'1raw'} = $totraw;

	my $final1 = $dir . $c . '/' . $pop . '/' . $pop . '_1_final.fastq.gz';
	my $final2 = $final1; $final2 =~ s/_1/_2/;
	my $finalu = $final1; $finalu =~ s/_1/_u/;
	my $fcount1 = fcount($final1);
	my $fcount2 = fcount($final2);
	my ($fcount3,$qual) = fcountu($finalu);

	my $totfin = $fcount1 + $fcount2 + $fcount3;
	$data{$c}{$pop}{'2fin'} = $totfin;
	$data{$c}{$pop}{'3qual'} = $qual;
    }
}

foreach my $c (keys %data) {
    foreach my $pop (keys %{$data{$c}}) {
	print $c, "\t", $pop, "\t";
	foreach my $values (sort {$a cmp $b} keys %{$data{$c}{$pop}}) {
	    print $data{$c}{$pop}{$values}, "\t";
	}
	print "\n";
    }
}

sub fcountu {
    my ($file) = @_;

    my $call1 = system("gunzip $file");
    $file =~ s/\.gz//;
    my $count; my $totqual; my $num;
    open(IN, "<$file");
    while(<IN>) {
	chomp(my $line = $_);
	chomp(my $seq =<IN>);
	chomp(my $id = <IN>);
	chomp(my $qual = <IN>);
	
	my $tmpqual;
	my @q = split(//,$qual);
	foreach my $q (@q) {
	    $tmpqual += ord($q);
	}
	$tmpqual = int($tmpqual/length($qual));
	$totqual += $tmpqual;
	$num++;

	$count += length($seq);
    }   
    close(IN);

    my $call3 = system("gzip $file");

    $totqual = int($totqual/$num) - 33;
    
    return($count,$totqual);
}




sub fcount {
    my ($file) = @_;

    my $call1 = system("gunzip $file");
    $file =~ s/\.gz//;
    my $count;
    open(IN, "<$file");
    while(<IN>) {
	chomp(my $line = $_);
	chomp(my $seq = <IN>);
	chomp(my $id = <IN>);
	chomp(my $qual = <IN>);
	
	$count += length($seq);
    }
    close(IN);

    my $call3 = system("gzip $file");

    return($count);
}

sub count {
    my ($file) = @_;
    
    my $call1 = system("gunzip $file");
    $file =~ s/\.gz//;
    my @call2 = `wc $file`;
    my $count = $1 if $call2[0] =~ /^\s+(\S+)/;
    my $call3 = system("gzip $file");
    
    return($count);
}

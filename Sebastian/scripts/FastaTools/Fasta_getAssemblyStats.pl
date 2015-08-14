#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

#use GD::Graph::bar;

if (!defined($ARGV[0]) ) { die "Missing arguments, please provided a fasta file" }

open(FASTAINPUT, $ARGV[0]) or die "Can't open file $ARGV[0]\n";
open(OUTPUT, ">$ARGV[0]_stats_report.txt");

my $sequenceSize = 0;

my $aCount = 0;
my $tCount = 0;
my $gCount = 0;
my $cCount = 0;
my $nCount = 0;
my $otherCount = 0;

my @sequenceSizes;
my $numSequences;

my $firstLine = <FASTAINPUT>;


while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\r//;

	if( $line =~ /^>/) {

		push(@sequenceSizes, $sequenceSize);
		++$numSequences;

		$sequenceSize = 0;

	} else {

		$sequenceSize += length($line);

		$aCount += ($line =~ tr/aA/aA/);
		$tCount += ($line =~ tr/tT/tT/);
		$gCount += ($line =~ tr/gG/gG/);
		$cCount += ($line =~ tr/cC/cC/);
		$nCount += ($line =~ tr/nN/nN/);

		$line =~ s/a|A|t|T|g|G|c|C|n|N//g;

		$otherCount = length($line);

	} #end else

} #end while

push(@sequenceSizes, $sequenceSize);
++$numSequences;

close(FASTAINPUT);

my $totalSequence = $aCount + $tCount + $gCount + $cCount + $nCount +$otherCount;

my $average = $totalSequence/ ($numSequences-1);

my $GCcontent = ( ($gCount + $cCount) / $totalSequence ) * 100;

my $nContent = ( $nCount / $totalSequence ) * 100;

my @sortedSizes = sort {$b <=> $a} @sequenceSizes;

my $N75_threshold = $totalSequence * 0.25;
my $N50_threshold = $totalSequence * 0.50;
my $N25_threshold = $totalSequence * 0.75;

my $N75;
my $N50;
my $N25;

my %sequenceDistr;
my %amountSequenceDistr;

my @SizeRanges = ('100', '200', '300', '400', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '6000', '7000', '8000', '9000', '10000', '15000', '20000','50000','100000');

my $cumulativeSize = 0;

foreach my $size (@sortedSizes) {

	$cumulativeSize += $size;

	if ( ($cumulativeSize > $N25_threshold) && (!defined($N75)) ) { $N75 = $size }
	if ( ($cumulativeSize > $N50_threshold) && (!defined($N50)) ) { $N50 = $size }
	if ( ($cumulativeSize > $N75_threshold) && (!defined($N25)) ) { $N25 = $size }

	foreach my $sizeLimit (@SizeRanges) {
		if ($sizeLimit > $size) {
			++$sequenceDistr{$sizeLimit};
			$amountSequenceDistr{$sizeLimit} += $size;
			last;
		} #end if range found

	} #end foreach ranges

	if ($size > 20000) {
		++$sequenceDistr{'100001'}; 
		$amountSequenceDistr{'100001'} += $size;
	} #end if big sequence

} #end foreach size

print OUTPUT "File analyzed:	$ARGV[0]\n\n";

print OUTPUT "N75 :	$N75\n";
print OUTPUT "N50 :	$N50\n";
print OUTPUT "N25 :	$N25\n";

print OUTPUT "Number of sequences :	$numSequences\n";

print OUTPUT "Average sequence length:	", sprintf("%.3f", $average), "\n";

my $maxSize = shift(@sortedSizes);
my $minSize = pop(@sortedSizes);

print OUTPUT "Min sequence Length :	$minSize\n";
print OUTPUT "Max sequence Length :	$maxSize\n\n";

print OUTPUT "Nucleotide composition\n";
print OUTPUT "Number of A's:	", $aCount, "\n";
print OUTPUT "Number of T's:	", $tCount, "\n";
print OUTPUT "Number of G's:	", $gCount, "\n";
print OUTPUT "Number of C's:	", $cCount, "\n";
print OUTPUT "Number of N's:	", $nCount, "\n";
print OUTPUT "Number of Non-ATGCN characters: ", $otherCount, "\n\n";

print OUTPUT "GC content :	", sprintf("%.2f", $GCcontent), "%\n";
print OUTPUT "N content :	", sprintf("%.2f", $nContent), "%\n";

print OUTPUT "Total number of bases: ", $totalSequence, "\n\n";



print OUTPUT "Sequence size and bases distribution per bins\n";
print OUTPUT "Bin\tNumber of Sequences\tNumber of bases (in the bin)\n";
my $lowerLimit = 0;
foreach my $upperLimit (@SizeRanges) {
	if(defined($sequenceDistr{$upperLimit})) {
		print OUTPUT $lowerLimit+1, "-", $upperLimit, "\t", $sequenceDistr{$upperLimit}, "\t", $amountSequenceDistr{$upperLimit}, "\n";
	} else {
		print OUTPUT $lowerLimit+1, "-", $upperLimit, "\t0\t0\n";
	} #end else

	$lowerLimit = $upperLimit;

} #end foreach ranges

if(defined($sequenceDistr{'20001'})) {
		print OUTPUT "+100001\t", $sequenceDistr{'100001'}, "\t", $amountSequenceDistr{'100001'}, "\n";
} #end if big sequences found

close(OUTPUT);

exit;


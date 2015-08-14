#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV < 1) {

	print "Two input parameters are needed, please input them.\n";
	print "Usage: Coverage_calculatorr.pl <Output from SAM_Ranges_Summarizer.pl> <Output prefix>\n";
	exit;

} #end if

open(TABLE, $ARGV[0]);
open(COVERAGES, ">$ARGV[1].coveragePerPositions.txt") or die "Can't open file $ARGV[1].coveragePerPositions.txt\n";
print COVERAGES "#Sequence\tPosition\tCoverage\n";

my %scaffolds;

# Load read data from table
while (my $line = <TABLE>) {

	chomp($line);
	my @data = split("\t", $line);

	push(@{$scaffolds{$data[1]}},"$data[3]-$data[4]");

} #end while
close(TABLE);

print "Alignment data loaded\n";

# Get number of scaffolds with reads
my $numScaffolds = keys %scaffolds;
print "Datapoints for $numScaffolds scaffolds were loaded\n";

# Initialize counting varaibles
my $countScaffolds = 0;
my $m = $numScaffolds/100;

# Initialize global variables to store information
my %coverages;
my %sequenceStats;

my $totalCount;
my $numPositions;

my $max = 1;

foreach my $scaffold (keys %scaffolds) {

	# Initialize local variables for counting
	my %coverage;
	my $totalCov;
	my $totalCovSequence;
	my $readsSequence;

	# Print progress update
	print STDERR int(100 * $countScaffolds / $numScaffolds), " " if $countScaffolds++ % $m == 0;

	# For each read increase the coverage count along the length of the read
	foreach my $read (@{$scaffolds{$scaffold}}) {

		my @coordinates = split("-", $read);

		++$readsSequence;

		for (my $i = $coordinates[0]; $i <= $coordinates[1]; ++$i) { ++$coverage{$i} }

	} #end foreach

	# For all the positions found with reads print and the coverage and add to global variables
	foreach my $pos (sort { $a <=> $b}  keys %coverage) {

		print COVERAGES $scaffold, "\t", $pos, "\t", $coverage{$pos}, "\n";

		if($coverage{$pos} > $max) {$max = $coverage{$pos}}

		push(@{$coverages{$scaffold}}, $coverage{$pos});

		$totalCov += $coverage{$pos};

		++$totalCovSequence;

	} #end foreach

	my $averageSeq = $totalCov/ ($totalCovSequence-1);

	$totalCount += $totalCov;
	$numPositions += $totalCovSequence;

	$sequenceStats{$scaffold} = [$averageSeq, $totalCovSequence, $readsSequence];

} #end foreach

print STDERR "\n";

print STDERR "Done calculating coverages, generating statistics\n";

my $average = $totalCount/ ($numPositions-1);

open(GENERAL, ">$ARGV[1].generalStats.txt") or die "Can't open file $ARGV[1].generalStats.txt\n";

open(SEQUENCE, ">$ARGV[1].statsPerSequence.txt") or die "Can't open file $ARGV[1].statsPerSequence.txt\n";
print SEQUENCE "#Sequence\tCoveredPositions\tAverageCoverage\tNumMappedReads\tStandardDeviation\n";

print GENERAL "Total number of positions with reads: $numPositions\n";
print GENERAL "Average number of reads per positions: ", sprintf("%.3f", $average), "\n";
print GENERAL "Maximun number of reads for a single positions: $max\n";

my $totalSqDev = 0;

foreach my $scaffold (keys %coverages) {

	my $totalSqDevScaffold = 0;

	foreach my $coverage (@{$coverages{$scaffold}} ) {

		$totalSqDev += ($average - $coverage) ** 2;

		$totalSqDevScaffold += ($sequenceStats{$scaffold}[0] - $coverage) ** 2;

	} #end foreach

	my $stDevScaffold = ($totalSqDevScaffold / ($sequenceStats{$scaffold}[1] - 1) ) ** 0.5;

	print SEQUENCE $scaffold, "\t", $sequenceStats{$scaffold}[1], "\t", sprintf("%.3f", $sequenceStats{$scaffold}[0]), "\t", $sequenceStats{$scaffold}[2], "\t", sprintf("%.4f", $stDevScaffold), "\n";

} #end foreach

my $stDev = ($totalSqDev / ($numPositions-1) ) ** 0.5;

print GENERAL "Standard deviation of the coverage across all the scaffolds: ", sprintf("%.4f", $stDev), "\n";

print STDERR "Completed!!\n\n";

close(COVERAGES);

exit;
	




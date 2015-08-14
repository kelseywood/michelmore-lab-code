#!/usr/bin/perl

use strict;
use warnings;

if( !defined($ARGV[1]) ) {
	print "FastQ_QualityCalculator.pl <Fastq file> <Quality threshold>\n";
	print "Script will calculate percentage of bases above the inputed threshold and mean quality for the entire read\n";
	print "Threshold should be numeric and within the ranges of the quality values (dependant of the quality format)\n";
	print "To find out which values correspond to each format please see http://en.wikipedia.org/wiki/FASTQ_format\n";
	print "OUTPUT:\n";
	print "Will generate a new file fastq.quality_scores.txt, with the calculated qualities\n";
	die "missing arguments\n"
}

my $fastq = $ARGV[0];
my $threshold = $ARGV[1];

# Getting distribution of reads across GC bins
open(FASTQ, $fastq) or die "Can't open file $fastq\n";
open(RESULTS, ">$fastq.quality_scores.txt");

print RESULTS "ReadName\tPercentage_bases_above_threshold\tMean_quality\n";

while(my $header = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	chomp($header);
	chomp($qual);

	if ( !($header =~ m/^@/) || !($header2 =~ /^\+/) ) { print "Issue with fastq format, please check your file formatting\n"; die}

	my @qualList = split("", $qual);
	my $qualSum = 0;
	my $goodBases = 0;

	foreach my $value (@qualList) {

		if(ord($value) >= $threshold) {++$goodBases}
		$qualSum += ord($value);

	} #end foreach

	print RESULTS $header, "\t", sprintf("%.2f", ($goodBases/length($qual)*100)), "\t", sprintf("%.2f", ($qualSum/length($qual))), "\n";

} #end while

close(FASTQ);
close(RESULTS);

exit;



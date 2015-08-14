#!/usr/bin/perl
#######################################################
#                                                     #
#                 Pile up to GFF3                     #
#                                                     #
#              Sebastian Reyes-Chin-Wo                #
#               sreyesch@ucdavis.edu                  #
#                                                     #
#######################################################
use strict;
use warnings;

if( !(defined($ARGV[0])) || (defined($ARGV[1]))) {

	print "Missing Arguments or extra arguments.\n";
	print "Please input the pileup file (only one file).\n";
	die "Missing Arguments";

} #end if

open(PILEUP, $ARGV[0]);
open(GFF, ">$ARGV[0].gff3");

my $windownSize = 1000;

my $firsLine = <PILEUP>;
chomp($firsLine);

my @firstData = split($firsLine);

my $prevScaffold = $firstData[0];

my $firstPos = $firstData[1];
my $prevPos = $firstData[1];

my $totalReads = 1 ;
my $totalPositions = $firstData[6];

my $readsCounts = 1;
my $numPositions = $firstData[6];

my @posDepth = [$firstData[6]];

my $numBases = 1;

while (my $line = <PILEUP>) {

	chomp($line);

	my @data = split($line);

	if($data[0] ne $prevScaffold) {

		my $averDepth = $readsCounts/$numPositions;

		my $dif = 0;

		foreach my $read (@posDepth) {

			$dif += abs($read - $averDepth);

		} #end foreach

		my $stdev = $dif/$numPositions;

		my @scaffoldName = split("_", $prevScaffold);

		my $scaffoldID = shift(@scaffoldName);

		print GFF $prevScaffold, "\t", "GenomicReads", "\t", "BGI", "\t", $firstPos, "\t", $prevPos, "\t", $averDepth, "\t", ".", "\t", ".", "\t",  "ID=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";Name=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";StDev=",$stdev,";NumPos=",$numPositions,";";

		$prevScaffold = $data[0];

		$firstPos = $data[1];

		$prevPos = $data[1];

		++$totalPositions;
		$totalReads += $data[6];

		$readsCounts = 1;
		$numPositions = $data[6];

		@posDepth = [$data[6]];

		$numBases = 1;

	} elsif ($numBases > $windownSize) {

		my $averDepth = $readsCounts/$numPositions;

		my $dif = 0;

		foreach my $read (@posDepth) {

			$dif += abs($read - $averDepth);

		} #end foreach

		my $stdev = $dif/$numPositions;

		my @scaffoldName = split("_", $prevScaffold);

		my $scaffoldID = shift(@scaffoldName);

		print GFF $prevScaffold, "\t", "GenomicReads", "\t", "BGI", "\t", $firstPos, "\t", $prevPos, "\t", $averDepth, "\t", ".", "\t", ".", "\t",  "ID=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";Name=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";StDev=",$stdev,";NumPos=",$numPositions,";";


		$firstPos = $data[1];

		$prevPos = $data[1];

		++$totalPositions;
		$totalReads += $data[6];

		$readsCounts = 1;
		$numPositions = $data[6];

		@posDepth = [$data[6]];

		$numBases = 1;

	} else {

		$prevPos = $data[1];

		++$totalPositions;
		$totalReads += $data[6];

		++$readsCounts;
		$numPositions += $data[6];

		 push(@posDepth, $data[6]);

		$numBases = 1;

	} #end else

} #end while

my $averDepth = $readsCounts/$numPositions;

my $dif = 0;

foreach my $read (@posDepth) {

	$dif += abs($read - $averDepth);

} #end foreach

my $stdev = $dif/$numPositions;

my @scaffoldName = split("_", $prevScaffold);

my $scaffoldID = shift(@scaffoldName);

print GFF $prevScaffold, "\t", "GenomicReads", "\t", "BGI", "\t", $firstPos, "\t", $prevPos, "\t", $averDepth, "\t", ".", "\t", ".", "\t",  "ID=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";Name=Depth-",$scaffoldID,"_",$firstPos,"-",$prevPos,";StDev=",$stdev,";NumPos=",$numPositions,";";

print "Main Stats for $ARGV[0]\n";
print "Total Number of positions = $totalPositions\n";
print "Average depth = ,", $totalReads/$totalPositions, "\n\n";

exit;






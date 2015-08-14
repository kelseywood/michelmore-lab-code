#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

if (!defined($ARGV[0])) {

	print "One argument is needed, please input it.\n";
	print "GeVo links file from SynMap\n";
        exit 0;

} #end if



open(MATCHES, $ARGV[0]);

my $inputFile = $ARGV[0];
$inputFile =~ s/\.[^.]+$//;

open(BLOCKS, ">$inputFile.links.txt");

print BLOCKS join("\t", ("sequence1", "start1", "end1", "sequence2", "start2", "end2", "numMatches", "meanIdentity", "sumLengths") );
print BLOCKS "\n";

my @blocks;

my $numMatches = 0;
my $sumIdentities = 0;
my $sumLengths = 0;

my $strand;

my $sequence1;
my $start1 = 0;
my $end1 = 0;

my $sequence2;
my $start2 = 0;
my $end2 = 0;

# Look for the first block
while (my $firstline = <MATCHES>) {

	if($firstline =~ /^#\d/) {

		chomp($firstline);

		my @blockInfo = split("\t", $firstline);

		$sumIdentities = $blockInfo[1];

		$sequence1 = $blockInfo[2];
		$sequence2 = $blockInfo[3];

		$sequence1 =~ s/^[a-z][0-9]+_//;
		$sequence2 =~ s/^[a-z][0-9]+_//;

		$strand = $blockInfo[4];

		$numMatches = $blockInfo[5];

		last;

	} #end if first block found

} #end of looking for first block


# read the whole file
while (my $line = <MATCHES>) {

	chomp($line);

	if ($line =~ /^#\d/) {

		my $meanIdentity = $sumIdentities/$numMatches;

		print BLOCKS join("\t", ($sequence1, $start1, $end1, $sequence2, $start2, $end2, $strand, $numMatches, sprintf("%.2f", $meanIdentity), $sumLengths) );
		print BLOCKS "\n";
		
		my @blockInfo = split("\t", $line);

		$numMatches = 0;
		$sumIdentities = 0;
		$sumLengths = 0;

		$sumIdentities = $blockInfo[1];

		$sequence1 = $blockInfo[2];
		$sequence2 = $blockInfo[3];

		$sequence1 =~ s/^[a-z][0-9]+_//;
		$sequence2 =~ s/^[a-z][0-9]+_//;

		$strand = $blockInfo[4];

		$numMatches = $blockInfo[5];

		$start1 = 0;
		$end1 = 0;
		$start2 = 0;
		$end2 = 0;

	} else {

		my @syntenicHit = split("\t", $line);

		$sumLengths += abs($syntenicHit[2] - $syntenicHit[3]);

		if($start1 == 0) { $start1 = $syntenicHit[2] }
		$end1 = $syntenicHit[3];

		if($strand eq "f") {
			if($start2 == 0) { $start2 = $syntenicHit[6] }
			$end2 = $syntenicHit[7];
		} elsif($strand eq "r") {
			if($end2 == 0) { $end2 = $syntenicHit[7] }
			$start2 = $syntenicHit[6];
		} #end strand
			

	} #end 

} #end while

close(MATCHES);
close(BLOCKS);

exit;

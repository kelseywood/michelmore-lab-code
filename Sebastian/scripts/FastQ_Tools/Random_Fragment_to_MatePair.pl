#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

############################################
#                                          #
#             #
#        Sebastian Reyes-Chin-Wo           #
#                                          #
#         sreyesch@ucdavis.edu             #
#                                          #
#                                          #
############################################

# This script was designed to generate mate pairs reads out of a fragments file

if (!defined($ARGV[3]) ) {
	print "Missing arguments\n";
	print "Usage Random_PairEnd_read_generator.pl <Fragment Fasta File> <Output File prefix> <Read length> <Quality score to asign to all the bases>\n";
	print "Maximum quality values:\n";
	print "  Sanger Phred+33 : I\n";
	print "  Illumina Phred+33 : J\n";
	print "  Illumina Phred+64 : h\n";
	print "NOTE:\n";
	print "Fasta Fragment file should contain sequence in one line\n";
	die;

} #end if insufficient arguments

my $fragments = $ARGV[0];
my $output = $ARGV[1];
my $readLength = $ARGV[2];
my $qualityScore = $ARGV[3];
my $variance = 0.10;

# Counting variables
my $readsGen;
my $numFragments;

# Generating the quality
my $qual;
for(my $j = 0; $j < $readLength,; ++$j) {$qual .= $qualityScore};

# Opening output files
open(FASTQ1, ">$output.1.fastq");
open(FASTQ2, ">$output.2.fastq");

## making the reads
print "Will start generating the reads\n";

# Opening input files
open(FRAGMENTS, $fragments);

while (my $fastaHeader = <FRAGMENTS>) {

	++$numFragments;

	chomp($fastaHeader);

	$fastaHeader =~ s/>/@/;

	my $seq = <FRAGMENTS>;
	chomp($seq);

	# Get forward read
	my $forward = reverse(substr($seq, 0, $readLength));

	# Get reverse read
	my $reverse = substr($seq, -$readLength, $readLength);

	#Verifying number of N's
	my $numN = ($forward =~ tr/Nn/Nn/) + ($reverse =~ tr/Nn/Nn/) ;
	if ($numN != 0) { next }

	$forward =~ tr/ACGTacgt/TGCAtgca/;

	# Print reads to fastq files
	print FASTQ1 $fastaHeader, "/1\n";
	print FASTQ1 $forward, "\n";
	print FASTQ1 "+\n";
	print FASTQ1 $qual, "\n";

	print FASTQ2 $fastaHeader, "/2\n";
	print FASTQ2 $reverse, "\n";
	print FASTQ2 "+\n";
	print FASTQ2 $qual, "\n";

	++$readsGen;

} #end while


print "Read generation concluded\n";
print $readsGen, " reads were generated\n";
print "Out of ", $numFragments," fragments\n";

close(FRAGMENTS);
close(FASTQ1);
close(FASTQ2);

exit;

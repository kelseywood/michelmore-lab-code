#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Bio::Seq;
use Bio::Index::Fasta;
use Math::Random qw(:all);

############################################
#                                          #
#      Random PairEnd read generator       #
#        Sebastian Reyes-Chin-Wo           #
#                                          #
#         sreyesch@ucdavis.edu             #
#                                          #
#                                          #
############################################

### Update from version 1
# Change sequence read to indexed format (usage of big genomes without need of a lot of memory
# Added function to add variance to insertsizes
# Added variance for the quality values (deprecated)


# This script was designed to generate random mate pairs reads out of a fasta file

if (!defined($ARGV[6]) ) {
	print "Missing arguments\n";
	print "Usage Random_PairEnd_read_generator.pl <Fasta File> <Output File prefix> <Fragment size> <Read length> <Number of reads> <Quality score to asign to all the bases> <Read name prefix>\n";
	print "Maximum quality values:\n";
	print "  Sanger Phred+33 : I\n";
	print "  Illumina Phred+33 : J\n";
	print "  Illumina Phred+64 : h\n";
	die;

} #end if insufficient arguments

my $fasta = $ARGV[0];
my $output = $ARGV[1];
my $fragmentSize = $ARGV[2];
my $readLength = $ARGV[3];
my $numberReads = $ARGV[4];
my $qualityScore = $ARGV[5];
my $readPrefix = $ARGV[6];
my $variance = 0.10;

### Format Fasta index

my $idx = Bio::Index::Fasta->new(
	                         '-filename' => "$fasta.idx",
	                         '-write_flag' => 1
	                        );

$idx->make_index("$fasta");

### Load fasta sequence

my @sequenceHeaders;

open(FASTA, $fasta);

while (my $fastaLine = <FASTA>) {

	chomp($fastaLine);

	if( $fastaLine =~ /^>/ ) {

		$fastaLine =~ s/^>//;

		my @header = split(" ", $fastaLine);

		push(@sequenceHeaders, $header[0]);

	} #end loading info

} #end while

close(FASTA);

## making the reads

print "Will start generating the reads\n";

# Opening files
open(FASTQ1, ">$output.1.fastq");
open(FASTQ2, ">$output.2.fastq");
open(FRAGMENTS, ">$output.0.fasta");

# Getting number of sequences
my $numSequences = scalar(@sequenceHeaders);

# Generating array of insert sizes
my @insertSizes = random_normal($numberReads, 0.01, 0.1);

# Generating the quality
my $qual;
for(my $j = 0; $j <= $readLength,; ++$j) {$qual .= $qualityScore};

for (my $i = 0; $i < $numberReads; ++$i) {

	# Calculate insert size
	my $insert = floor($fragmentSize + ($fragmentSize*$insertSizes[$i]));

	my $status = "unprinted";

	while($status eq "unprinted") {

		# Fetch sequence from index
		my $seqIndex = random_uniform(1, 0, $numSequences);
#		print $sequenceHeaders[$seqIndex], "\n";
		my $seq = $idx->fetch($sequenceHeaders[$seqIndex]);

		# Verify that sequence it's big enough
		if( $seq->length <= $insert ) {next}

		# Get a random start position for the fragment
		my $randPos = floor(rand(($seq->length)-$insert));

		# Verifying proper start
		if ($randPos <= 0) {print $randPos; next}

		# Get forward read
		my $forward = $seq->subseq($randPos, ($randPos+$readLength) );

		# Get reverse read, "$randPos+$fragmentLength-$readLength" it's the start point for the reverse read
		my $reverse = $seq->subseq(($randPos+$insert-$readLength), ($randPos+$insert) );

		#Verifying number of N's
		my $numN = ($forward =~ tr/Nn/Nn/) + ($reverse =~ tr/Nn/Nn/) ;
		if ($numN != 0) { next }

		$reverse =~ tr/ACGTacgt/TGCAtgca/;

		# Print reads to fastq files
		print FASTQ1 "@", $readPrefix, "_", $i, "/1\n";
		print FASTQ1 $forward, "\n";
		print FASTQ1 "+\n";
		print FASTQ1 $qual, "\n";

		print FASTQ2 "@", $readPrefix, "_", $i, "/2\n";
		print FASTQ2 $reverse, "\n";
		print FASTQ2 "+\n";
		print FASTQ2 $qual, "\n";

		print FRAGMENTS ">", $readPrefix, "_", $i, "\n";
		print FRAGMENTS $seq->subseq($randPos, ($randPos+$insert) ), "\n";

		$status = "printed";

	} #end while generating reads

} #end selecting reads

print "Reads succesfully generated\n";

close(FASTQ1);
close(FASTQ2);

exit;

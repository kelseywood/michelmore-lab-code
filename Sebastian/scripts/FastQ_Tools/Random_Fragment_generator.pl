#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Bio::Seq;
use Bio::Index::Fasta;
use Math::Random qw(:all);

############################################
#                                          #
#      r       #
#        Sebastian Reyes-Chin-Wo           #
#                                          #
#         sreyesch@ucdavis.edu             #
#                                          #
#                                          #
############################################

# This script was designed to generate random fragmetns out of a fasta file

if (!defined($ARGV[4]) ) {
	print "Missing arguments\n";
	print "Usage Random_PairEnd_read_generator.pl <Fasta File> <Output File prefix> <Fragment size> <Number of fragmetns> <Fragment name prefix>\n";
	die;

} #end if insufficient arguments

my $fasta = $ARGV[0];
my $output = $ARGV[1];
my $fragmentSize = $ARGV[2];
my $numberFragments = $ARGV[3];
my $fragmentPrefix = $ARGV[4];
my $variance = 0.10;

### Format Fasta index

my $idx = Bio::Index::Fasta->new(
	                         '-filename' => "$fasta.idx",
	                         '-write_flag' => 1
	                        );

$idx->make_index("$fasta");

### Load fasta sequence

print STDERR "Loading fasta headers\n";

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

print "Will start generating the fragments\n";

# Opening files
open(FRAGMENTS, ">$output.0.fasta");

# Getting number of sequences
my $numSequences = scalar(@sequenceHeaders);

# Generating array of insert sizes
my @insertSizes = random_normal($numberFragments, 0.01, 0.1);

my $generatedRead = 0;
my $m = int $numberFragments/100;
#When less than 100 reads are generated m need to be corrected
if($m == 0) {$m = 1};

for (my $i = 0; $i < $numberFragments; ++$i) {

	#print advance status
	print STDERR int(100 * $generatedRead / $numberFragments), " " if $generatedRead++ % $m == 0;	

	# Calculate insert size
	my $insert = floor($fragmentSize + ($fragmentSize*$insertSizes[$i]));

	my $status = "unprinted";

	while($status eq "unprinted") {

		# Fetch sequence from index
		my $seqIndex = random_uniform(1, 0, $numSequences);
		my $seq = $idx->fetch($sequenceHeaders[$seqIndex]);

		# Verify that sequence it's big enough
		if( $seq->length <= $insert ) {next}

		# Get a random start position for the fragment
		my $randPos = floor(rand(($seq->length)-$insert));

		# Verifying proper start
		if ($randPos <= 0) {print $randPos; next}

		print FRAGMENTS ">", $fragmentPrefix, "_", $i, " FragmentLength ", $insert, "\n";
		print FRAGMENTS $seq->subseq($randPos, ($randPos+$insert) ), "\n";

		$status = "printed";

	} #end while generating reads

} #end selecting reads

print "\nFragments succesfully generated\n";

close(FRAGMENTS);

exit;

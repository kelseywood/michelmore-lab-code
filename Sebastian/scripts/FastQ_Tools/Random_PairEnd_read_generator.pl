#!/usr/bin/perl

use strict;
use warnings;
use POSIX;


############################################
#                                          #
#      Random PairEnd read generator       #
#        Sebastian Reyes-Chin-Wo           #
#                                          #
#         sreyesch@ucdavis.edu             #
#                                          #
#                                          #
############################################

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

### Load fasta sequence

my @sequence;
my $possiblePairs;
my $totalSequence;

open(FASTA, $fasta);

while (my $fastaLine = <FASTA>) {

  chomp($fastaLine);

  if( !($fastaLine =~ /^>/) && ($fastaLine =~ /a|c|g|t|A|G|C|T/) ) {

    push(@sequence, $fastaLine);

    my $length = length($fastaLine);

    $possiblePairs += ($length - $fragmentSize)/2;

    $totalSequence += $length;

  } #end if

} #end while

close(FASTA);

print "Total amount of sequence in the fasta file ", $totalSequence, "\n";
print "Total number of 2 bases stack pairs ", $possiblePairs, "\n";


## making the reads

print "Will start generating the reads\n";

open(FASTQ1, ">$output.1.fastq");
open(FASTQ2, ">$output.2.fastq");
open(FRAGMENTS, ">$output.0.fasta");

my @quals;

for(my $j = 0; $j < $readLength,; ++$j) {push(@quals, $qualityScore)};

#my $qual =~ s/^/'$qualityScore' x $readLength/e;;

my $numSequences = scalar(@sequence);

for (my $i = 1; $i <= $numberReads; ++$i) {

  # Get a random index for the sequence
  my $randSeq = floor(rand($numSequences));

  # Get a random start position for the fragment
  my $randPos = floor(rand(length($sequence[$randSeq])-$fragmentSize));

  # Get forward read
  my $forward = substr($sequence[$randSeq], $randPos, $readLength);

  # Get reverse read, "$randPos+$fragmentLength-$readLength" it's the start point for the reverse read
  my $reverse = reverse(substr($sequence[$randSeq], ($randPos+$fragmentSize-$readLength), $readLength));

  $reverse =~ tr/ACGTacgt/TGCAtgca/;

#  print $randPos, "\t", $randPos+$fragmentSize-$readLength, "\t", $randPos+$fragmentSize, "\n";

  # Print reads to fastq files
  print FASTQ1 "@", $readPrefix, "_", $i, "/1\n";
  print FASTQ1 $forward, "\n";
  print FASTQ1 "+\n";
  print FASTQ1 @quals, "\n";

  print FASTQ2 "@", $readPrefix, "_", $i, "/2\n";
  print FASTQ2 $reverse, "\n";
  print FASTQ2 "+\n";
  print FASTQ2 @quals, "\n";

  print FRAGMENTS substr($sequence[$randSeq], $randPos, $fragmentSize), "\n";

} #end selecting reads

print "Reads succesfully generated\n";

close(FASTQ1);
close(FASTQ2);

exit;

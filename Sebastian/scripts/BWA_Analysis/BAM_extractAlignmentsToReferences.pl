#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[2])) {

	print "Error, please inpute a text file with the references to extract, a BAM file (sorted or unsorted) and an output file prefix)\n";
	die "Missing arguments\n";

} #end if

if ( !(-f "$ARGV[1].bai") ) {

	print "Will generate bam index\n";

	if( system("/home/sreyesch/MichelmoreBin/bin/samtools.1.18 index $ARGV[1]") != 0 ) {

		print STDERR "ERROR: Couldn't generate BAM indes.\n";
		die;

	} #end if

} #end generating bam index

if( system("/home/sreyesch/MichelmoreBin/bin/samtools.1.18 view -H $ARGV[1] > $ARGV[2].sam") != 0 ) {

	print STDERR "ERROR: Couldn't retrieve SAM header from $ARGV[1].\n";

} #end if

open(LIST, $ARGV[0]) or die "Can't open file $ARGV[0]\n";
while (my $reference = <LIST>) {

	chomp($reference);

	if( system("/home/sreyesch/MichelmoreBin/bin/samtools.1.18 view $ARGV[1] $reference >> $ARGV[2].sam") != 0 ) {

		print STDERR "ERROR: Couldn't retrieve alignments for $reference.\n";

	} #end if



} #end while
close(LIST);

if( system("/home/sreyesch/MichelmoreBin/bin/samtools.1.18 view -b -S -o $ARGV[2].bam $ARGV[2].sam") != 0 ) {

	print "ERROR: Couldn't transform SAM to BAM.\n";

} #end if

exit;














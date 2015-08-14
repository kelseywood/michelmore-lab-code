#!/usr/bin/perl
use strict; use warnings;

if (scalar(@ARGV) <= 2) { die "usage:  <fastq> <kmer start> <kmer length>\n"}
my $fastq = $ARGV[0];
my $kmerStart = $ARGV[1];
my $kmerLength = $ARGV[2];

--$kmerStart;

my %barcodes;

open(FASTQ, $fastq) or die "Can't open $fastq\n";

while (my $fastqHeader = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $secondeHeader = <FASTQ>;
	my $qual = <FASTQ>;

	my $barcode = substr($sequence, $kmerStart, $kmerLength);

	++$barcodes{$barcode};	

} #end while

print "Statistics of the barcodes\n";

print length(keys %barcodes), " were found\n";

foreach my $barcode ( sort hashValueDescendingNum ( keys (%barcodes) ) ) {

	print $barcode, "\t", $barcodes{$barcode}, "\n";

} #end foreach

exit;

sub hashValueDescendingNum {
   $barcodes{$b} <=> $barcodes{$a};
} #end of hashValueDescendingNum sub


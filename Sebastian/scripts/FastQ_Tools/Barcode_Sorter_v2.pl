#!/usr/bin/perl
#use strict;
use warnings;

### Diff from version 0, add capability to use a different length of barcode


if (scalar(@ARGV) <= 3) { die "usage:  <barcodes names> <length of barcode> <third read> <fastq reverse> <fastq forward>\n"}
my $barcodeFile = $ARGV[0]; #Tab delimited file with the barcodes and sample names (two columns, samplename and barcode)
my $barcodeLength = $ARGV[1];
my $third = $ARGV[2];
my $fastq1 = $ARGV[3];
my $fastq2 = $ARGV[4];

my %barcodes;
my $countBarcodes;

open(BARCODES, $barcodeFile);# or die "Can't open file $barcodeFile";

while (my $barcodeLine = <BARCODES>) {

	chomp($barcodeLine);

        my @barcode = split("\t", $barcodeLine);

        my $code = substr($barcode[1],0,$barcodeLength);

        $barcodes{$code}{'reference'} = $barcode[0];

        my $outputFile1 = $barcode[0]."_1.fastq";
        my $outputFile2 = $barcode[0]."_2.fastq";

        my $filehandle1 = $code."_1";
        my $filehandle2 = $code."_2";

        open( $filehandle1, ">$outputFile1");
        open( $filehandle2, ">$outputFile2");

        ++$countBarcodes;

} #end while for barcodes

open(THIRD, $third);
open(FASTQ1, $fastq1);
open(FASTQ2, $fastq2);

while (my $thirdHeader = <THIRD>) {

      	#Get all the lines from the third read file

	my $thirdSeq = <THIRD>;
	my $thirdHeader2 = <THIRD>;
	my $thirdQual = <THIRD>;

	chomp($thirdSeq);

        my $code = substr($thirdSeq,0,$barcodeLength);

        my $filehandle1 = $code."_1";
        my $filehandle2 = $code."_2";

	++$barcodes{$code}{'counts'};

	#Check if barcode as a reference
	if(defined($barcodes{$code}{'reference'}) ) {

		#Get and print lines from forward and reverse

		for(my $x = 0; $x < 4; ++$x) {

			my $fastq1_line = <FASTQ1>;
		        print $filehandle1 $fastq1_line;

			my $fastq2_line = <FASTQ2>;
		        print $filehandle2 $fastq2_line;

	      	} #end for fastq printing

	} #end if

} #end while

print "Statistics of the barcodes\n";

print $countBarcodes, " were found\n";

foreach my $barcode (sort { $barcodes{$b}->{'counts'} <=> $barcodes{$a}->{'counts'} } keys %barcodes ) {

	if(defined($barcodes{$barcode}{'reference'}) ) {

		print $barcode, "\t", $barcodes{$barcode}{'counts'}, "\t", $barcodes{$barcode}{'reference'}, "\n";

	} else {

		print $barcode, "\t", $barcodes{$barcode}{'counts'}, "\tNo_sample\n";

	} #end else

} #end foreach

exit;

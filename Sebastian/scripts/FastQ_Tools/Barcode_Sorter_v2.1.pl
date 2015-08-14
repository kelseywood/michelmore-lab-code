#!/usr/bin/perl
#use strict;
use warnings;

### Diff from version 0, add capability to use a different length of barcode


if (scalar(@ARGV) <= 4) { die "usage:  <barcodes names> <length of barcode> <third read> <fastq reverse> <fastq forward>\n"}
my $barcodeFile = $ARGV[0]; #Tab delimited file with the barcodes and sample names (two columns, samplename and barcode)
my $barcodeLength = $ARGV[1];
my $third = $ARGV[2];
my $fastq1 = $ARGV[3];
my $fastq2 = $ARGV[4];

my %barcodes;
my $countBarcodes;
my @perfectBarcodes;
my @allBarcodes;

open(LOG, ">$ARGV[0].log");

print LOG "Analysis started at ", getTime(), "\n\n";

open(BARCODES, $barcodeFile) or die "Can't open file $barcodeFile";

while (my $barcodeLine = <BARCODES>) {

	chomp($barcodeLine);

        my @barcode = split("\t", $barcodeLine);

        my $origCode = substr($barcode[1],0,$barcodeLength);

        my $outputFile1 = $barcode[0]."_1.fastq";
        my $outputFile2 = $barcode[0]."_2.fastq";

        my $filehandle1 = $origCode."_1";
        my $filehandle2 = $origCode."_2";

        open( $filehandle1, ">$outputFile1");
        open( $filehandle2, ">$outputFile2");

	push(@perfectBarcodes, $origCode);

	my @bases = ("A", "C", "G", "T", "N");

	for( my $i = 0; $i < $barcodeLength; ++$i) {

		foreach my $base (@bases) {

			my $varCode = $origCode;

			substr($varCode, $i, 1, $base);

			push(@allBarcodes, $varCode);

		        $barcodes{$varCode}{'reference'} = $barcode[0];

		        $barcodes{$varCode}{'barcode'} = $origCode;

			$barcodes{$varCode}{'counts'} = 0;

		} #end foreach bases

	} #end for length barcode

        ++$countBarcodes;

} #end while for barcodes

print  "Barcodes loaded\n\n";

open(THIRD, $third);
open(FASTQ1, $fastq1);
open(FASTQ2, $fastq2);

my $totalReads = 0;

while (my $thirdHeader = <THIRD>) {

	++$totalReads;

      	#Get all the lines from the third read file

	my $thirdSeq = <THIRD>;
	my $thirdHeader2 = <THIRD>;
	my $thirdQual = <THIRD>;

	chomp($thirdSeq);

        my $code = substr($thirdSeq,0,$barcodeLength);

	++$barcodes{$code}{'counts'};

	if(defined($barcodes{$code}{'barcode'}) ) {

	        my $filehandle1 = $barcodes{$code}{'barcode'}."_1";
	        my $filehandle2 = $barcodes{$code}{'barcode'}."_2";

		 #Get and print lines from forward and reverse
	        for(my $x = 0; $x < 4; ++$x) {

		        my $fastq1_line = <FASTQ1>;
	               print $filehandle1 $fastq1_line;

	       	 my $fastq2_line = <FASTQ2>;
       	        print $filehandle2 $fastq2_line;

	      	} #end for fastq printing

	} #end if barcode present

} #end while

print "DeMultiplexing completted, preparing statistics\n\n";

my $totalMatchCodes = 0;
foreach my $code (@allBarcodes) { if( defined($barcodes{$code}{'counts'}) ) { $totalMatchCodes += $barcodes{$code}{'counts'} } } #end foreach

my $perfectMatchCodes = 0;
foreach my $code (@perfectBarcodes) { if( defined($barcodes{$code}{'counts'}) ) { $perfectMatchCodes += $barcodes{$code}{'counts'} } } #end foreach


print LOG "Statistics of the barcodes\n";

print LOG $countBarcodes, " are been used\n";
print LOG scalar(@allBarcodes), " variable codes were generated\n\n";
print LOG "Total number of reads: ", $totalReads, "\n";
print LOG "Total number of deMultiplex reads: ", $totalMatchCodes, "\n";
print LOG "Total number of reads with perfect barcode: ", $perfectMatchCodes, "\n\n";

print LOG "Percentage of succesfully demultiplex reads: ", ($totalMatchCodes/$totalReads)*100, "\n";
print LOG "Percentage of reads with perfect matched barcode: ", ($perfectMatchCodes/$totalReads)*100, "\n";

#print LOG join(",", @perfectBarcodes), "\n\n";
#print LOG join(",", @allBarcodes), "\n\n";


foreach my $barcode (sort { $barcodes{$b}->{'counts'} <=> $barcodes{$a}->{'counts'} } keys %barcodes ) { 

	if(defined($barcodes{$barcode}{'reference'}) ) {

		print LOG $barcode, "\t", $barcodes{$barcode}{'counts'}, "\t", $barcodes{$barcode}{'reference'}, "\n";

	} else {

		print LOG $barcode, "\t", $barcodes{$barcode}{'counts'}, "\tNo_sample\n";

	} #end else
} #end foreach

print LOG "Analysis finished at ", getTime(), "\n";

print "Analysis completed\n\n";

exit;

####################################
#                                  #
#          Get time sub            #
#                                  #
####################################

sub getTime {

	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	return "$hour:$minute:$second, $months[$month] $dayOfMonth, $year";

} #end of getTime sub


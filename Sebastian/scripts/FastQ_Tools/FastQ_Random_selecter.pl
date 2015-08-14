#!/usr/bin/perl

use strict; use warnings;
use POSIX;

if( !defined($ARGV[2]) ) {
	print "Please provide a fastq file and a maximum number of N's allowed per read\n";
	print "FastQ_N_filterer <Fastq forward> <Fastq reverse> <Fraction to select>\n";
	die "missing argumetns\n"
}

my $fractionOfReadsNeedIt = $ARGV[2];

if ( (0 > $fractionOfReadsNeedIt) && ( $fractionOfReadsNeedIt >= 1 ) ) {

	print "Fraction of reads should be between 0 and 1\n";
	print "Please input a fraction\n";

} #end checking fraction its a fraction


# Reopen fastq file and create output files
open(FASTQ1, $ARGV[0]);
my $file1 = $ARGV[0];
$file1 =~ s{\.[^.]+$}{};
open(SELECTEDFASTQ1, ">$file1.selected.$fractionOfReadsNeedIt.fastq");

open(FASTQ2, $ARGV[1]);
my $file2 = $ARGV[1];
$file2 =~ s{\.[^.]+$}{};
open(SELECTEDFASTQ2, ">$file2.selected.$fractionOfReadsNeedIt.fastq");

print "Distribution set and fractions calculated, will start selecting the reads, please wait\n\n";

# Select and print reads
while(my $header1 = <FASTQ1>) {

	my $sequence1 = <FASTQ1>;
	my $qualheader1 = <FASTQ1>;
	my $qual1 = <FASTQ1>;

	my $header2 = <FASTQ2>;
	my $sequence2 = <FASTQ2>;
	my $qualheader2 = <FASTQ2>;
	my $qual2 = <FASTQ2>;

	my $probability = getRandomUnif($fractionOfReadsNeedIt);

	if($probability == 1) {

		print SELECTEDFASTQ1 $header1;
		print SELECTEDFASTQ1 $sequence1;
		print SELECTEDFASTQ1 $qualheader1;
		print SELECTEDFASTQ1 $qual1;

		print SELECTEDFASTQ2 $header2;
		print SELECTEDFASTQ2 $sequence2;
		print SELECTEDFASTQ2 $qualheader2;
		print SELECTEDFASTQ2 $qual2;

	} #end of selecting

} #end while

close(FASTQ1);
close(FASTQ2);
close(SELECTEDFASTQ1);
close(SELECTEDFASTQ2);

print "Read selection completed\n\n";

exit;

###############################
#
# Subroutines
#
###############################

# This function selects a random number between 0 and 1 and return '1' if it falls between the fraction that need to be selected, considering that the fraction that needs to be selected its between 0 and $prob
sub getRandomUnif {

	my ($prob) = @_;

	my $point = rand(1);

	if($point <= $prob) {return 1 }

	else {return 0}

} #end getRandomUnif sub


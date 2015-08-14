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

	print "Fraction of reads need should be between 0 and 1\n";
	print "Please input a fraction\n";

} #end checking fraction its a fraction


#Initialize read counts

my $totalReads;

my $GC10_30;
my $GC30_40;
my $GC40_45;
my $GC45_50;
my $GC50_60;
my $GC60_80;

print "Determining the distribution of the reads\n\n";

# Getting distribution of reads across GC bins
open(FASTQ1, $ARGV[0]);
while(my $header = <FASTQ1>) {

	my $sequence = <FASTQ1>;
	my $header2 = <FASTQ1>;
	my $qual = <FASTQ1>;

	chomp($sequence);

	my $length = length($sequence);

	my $GCcount = $sequence =~ tr/(g|G|c|C)//;

	my $GCcontent = ($GCcount/$length)*100;

	if ( (10 <= $GCcontent) && ($GCcontent < 30) ) { ++$GC10_30; ++$totalReads }

	elsif ( (30 <= $GCcontent) && ($GCcontent < 40) ) { ++$GC30_40; ++$totalReads }

	elsif ( (40 <= $GCcontent) && ($GCcontent < 45) ) { ++$GC40_45; ++$totalReads }

	elsif ( (45 <= $GCcontent) && ($GCcontent < 50) ) { ++$GC45_50; ++$totalReads }

	elsif ( (50 <= $GCcontent) && ($GCcontent < 60) ) { ++$GC50_60; ++$totalReads }

	elsif ( (60 <= $GCcontent) && ($GCcontent <= 80) ) { ++$GC60_80; ++$totalReads } #end elsif GC content

} #end while

close(FASTQ);

my $readsNeedPerOneGCpercent = ($totalReads*$fractionOfReadsNeedIt)/100;

# Calculate how many reads you need per bin (amount to be extracted it's normalized to a 100% among the 71percentiles)
my $readsNeedIt10_30 = ($readsNeedPerOneGCpercent * 28);
my $readsNeedIt30_40 = ($readsNeedPerOneGCpercent * 14);
my $readsNeedIt40_45 = ($readsNeedPerOneGCpercent * 7.5);
my $readsNeedIt45_50 = ($readsNeedPerOneGCpercent * 7.5);
my $readsNeedIt50_60 = ($readsNeedPerOneGCpercent * 14);
my $readsNeedIt60_80 = ($readsNeedPerOneGCpercent * 29);

# Calculate fraction of each bin it's need it
my $fractionGC10_30 = $readsNeedIt10_30/$GC10_30;
my $fractionGC30_40 = $readsNeedIt30_40/$GC30_40;
my $fractionGC40_45 = $readsNeedIt40_45/$GC40_45;
my $fractionGC45_50 = $readsNeedIt45_50/$GC45_50;
my $fractionGC50_60 = $readsNeedIt50_60/$GC50_60;
my $fractionGC60_80 = $readsNeedIt60_80/$GC60_80;


# adjusting the fractions in case of decrease number of reads
my $remainingReads;

my $binLessReads = 0;

if($fractionGC10_30 > 1) { 
	$fractionGC10_30 = 1;
	$remainingReads += $readsNeedIt10_30-$GC10_30;
	++$binLessReads;
} #end

if($fractionGC30_40 > 1) { 
	$fractionGC30_40 = 1;
	$remainingReads += $readsNeedIt30_40-$GC30_40;
	++$binLessReads;
} #end

if($fractionGC50_60 > 1) { 
	$fractionGC50_60 = 1;
	$remainingReads += $readsNeedIt50_60-$GC50_60;
	++$binLessReads;
} #end

if($fractionGC60_80 > 1) { 
	$fractionGC60_80 = 1;
	$remainingReads += $readsNeedIt60_80-$GC60_80;
	++$binLessReads;
} #end


# in case that the tail don't have enough reads, add the extra percentages to the center peak
if ($binLessReads > 0) {

	if($fractionGC30_40 < 1) { $fractionGC30_40 += ceil($remainingReads/(6-$binLessReads))/$GC30_40 }
	if($fractionGC40_45 < 1) { $fractionGC40_45 += ceil($remainingReads/(6-$binLessReads))/$GC40_45 }
	if($fractionGC45_50 < 1) { $fractionGC45_50 += ceil($remainingReads/(6-$binLessReads))/$GC45_50 }
	if($fractionGC50_60 < 1) { $fractionGC50_60 += ceil($remainingReads/(6-$binLessReads))/$GC50_60 }

} #end if lacking reads

# Counters of selected reads
my $totalSelected;

my $selectedGC10_30;
my $selectedGC30_40;
my $selectedGC40_45;
my $selectedGC45_50;
my $selectedGC50_60;
my $selectedGC60_80;

# Reopen fastq file and create output files
open(FASTQ1, $ARGV[0]);
my $file1 = $ARGV[0];
$file1 =~ s{\.[^.]+$}{};
open(SELECTEDFASTQ1, ">$file1.selected.fastq");

open(FASTQ2, $ARGV[1]);
my $file2 = $ARGV[1];
$file2 =~ s{\.[^.]+$}{};
open(SELECTEDFASTQ2, ">$file2.selected.fastq");

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

	chomp($sequence1);

	my $length = length($sequence1);

	my $GCcount = $sequence1 =~ tr/(g|G|c|C)//;

	my $GCcontent = ($GCcount/$length)*100;

	if ( (10 <= $GCcontent) && ($GCcontent < 30) ) {

		my $probability = getRandomUnif($fractionGC10_30);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC10_30;
			++$totalSelected;

		} #end if to print

	} elsif ( (30 <= $GCcontent) && ($GCcontent < 40) ) {

		my $probability = getRandomUnif($fractionGC30_40);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC30_40;
			++$totalSelected;

		} #end if to print

	} elsif ( (40 <= $GCcontent) && ($GCcontent < 45) ) {

		my $probability = getRandomUnif($fractionGC40_45);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC40_45;
			++$totalSelected;

		} #end if to print

	} elsif ( (45 <= $GCcontent) && ($GCcontent < 50) ) {

		my $probability = getRandomUnif($fractionGC45_50);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC45_50;
			++$totalSelected;

		} #end if to print

	} elsif ( (50 <= $GCcontent) && ($GCcontent < 60) ) {

		my $probability = getRandomUnif($fractionGC50_60);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC50_60;
			++$totalSelected;

		} #end if to print

	} elsif ( (60 <= $GCcontent) && ($GCcontent <= 80) ) {

		my $probability = getRandomUnif($fractionGC60_80);

		if($probability == 1) {

			printFastq($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2);

			++$selectedGC60_80;
			++$totalSelected;

		} #end if to print

	} #end elsif GC content

} #end while

close(FASTQ);
close(SELECTEDFASTQ1);
close(SELECTEDFASTQ2);

print "Read selection completed, preparing stats file\n\n";

open(STATS, ">$file1.selectionStats.txt");

print STATS "TOTAL NUMBERS\n\n";

print STATS "Total number of reads :						", $totalReads, "\n";
print STATS "Number of reads between 10-30 GC content :			", $GC10_30, "\n";
print STATS "Number of reads between 30-40 GC content :			", $GC30_40, "\n";
print STATS "Number of reads between 40-45 GC content :			", $GC40_45, "\n";
print STATS "Number of reads between 45-50 GC content :			", $GC45_50, "\n";
print STATS "Number of reads between 50-60 GC content :			", $GC50_60, "\n";
print STATS "Number of reads between 60-80 GC content :			", $GC60_80, "\n\n";

print STATS "Fraction of reads between 10-30 GC content :			", sprintf("%.3f", $GC10_30/$totalReads), "\n";
print STATS "Fraction of reads between 30-40 GC content :			", sprintf("%.3f", $GC30_40/$totalReads), "\n";
print STATS "Fraction of reads between 40-45 GC content :			", sprintf("%.3f", $GC40_45/$totalReads), "\n";
print STATS "Fraction of reads between 45-50 GC content :			", sprintf("%.3f", $GC45_50/$totalReads), "\n";
print STATS "Fraction of reads between 50-60 GC content :			", sprintf("%.3f", $GC50_60/$totalReads), "\n";
print STATS "Fraction of reads between 60-80 GC content :			", sprintf("%.3f", $GC60_80/$totalReads), "\n\n\n";


print STATS "ESTIMATED NUMBERS OF SELECTION\n\n";

print STATS "Number of reads to be selected :				", sprintf("%.1f", $fractionOfReadsNeedIt*$totalReads), "\n";
print STATS "Number of reads to be selected between 10-30 GC content :	", sprintf("%.1f", $readsNeedIt10_30), "\n";
print STATS "Number of reads to be selected between 30-40 GC content :	", sprintf("%.1f", $readsNeedIt30_40), "\n";
print STATS "Number of reads to be selected between 40-45 GC content :	", sprintf("%.1f", $readsNeedIt40_45), "\n";
print STATS "Number of reads to be selected between 45-50 GC content :	", sprintf("%.1f", $readsNeedIt45_50), "\n";
print STATS "Number of reads to be selected between 50-60 GC content :	", sprintf("%.1f", $readsNeedIt50_60), "\n";
print STATS "Number of reads to be selected between 60-80 GC content :	", sprintf("%.1f", $readsNeedIt60_80), "\n\n";

print STATS "Fraction of reads to be selected between 10-30 GC content :	", sprintf("%.3f", $fractionGC10_30), "\n";
print STATS "Fraction of reads to be selected between 30-40 GC content :	", sprintf("%.3f", $fractionGC30_40), "\n";
print STATS "Fraction of reads to be selected between 40-45 GC content :	", sprintf("%.3f", $fractionGC40_45), "\n";
print STATS "Fraction of reads to be selected between 45-50 GC content :	", sprintf("%.3f", $fractionGC45_50), "\n";
print STATS "Fraction of reads to be selected between 50-60 GC content :	", sprintf("%.3f", $fractionGC50_60), "\n";
print STATS "Fraction of reads to be selected between 60-80 GC content :	", sprintf("%.3f", $fractionGC60_80), "\n\n\n";


print STATS "NUMBERS OF SELECTED READS\n\n";

print STATS "Number of selected reads :					", $totalSelected, "\n";
print STATS "Number of selected reads between 10-30 GC content :		", $selectedGC10_30, "\n";
print STATS "Number of selected reads between 30-40 GC content :		", $selectedGC30_40, "\n";
print STATS "Number of selected reads between 40-45 GC content :		", $selectedGC40_45, "\n";
print STATS "Number of selected reads between 45-50 GC content :		", $selectedGC45_50, "\n";
print STATS "Number of selected reads between 50-60 GC content :		", $selectedGC50_60, "\n";
print STATS "Number of selected reads between 60-80 GC content :		", $selectedGC60_80, "\n\n";

print STATS "Fraction of selected reads between 10-30 GC content :		", sprintf("%.3f", $selectedGC10_30/$totalSelected), "\n";
print STATS "Fraction of selected reads between 30-40 GC content :		", sprintf("%.3f", $selectedGC30_40/$totalSelected), "\n";
print STATS "Fraction of selected reads between 40-45 GC content :		", sprintf("%.3f", $selectedGC40_45/$totalSelected), "\n";
print STATS "Fraction of selected reads between 45-50 GC content :		", sprintf("%.3f", $selectedGC45_50/$totalSelected), "\n";
print STATS "Fraction of selected reads between 50-60 GC content :		", sprintf("%.3f", $selectedGC50_60/$totalSelected), "\n";
print STATS "Fraction of selected reads between 60-80 GC content :		", sprintf("%.3f", $selectedGC60_80/$totalSelected), "\n\n";

close(STATS);

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
	
sub printFastq {

	my ($header1,$sequence1,$qualheader1,$qual1,$header2,$sequence2,$qualheader2,$qual2) = @_;

	print SELECTEDFASTQ1 $header1;
	print SELECTEDFASTQ1 $sequence1, "\n";
	print SELECTEDFASTQ1 $qualheader1;
	print SELECTEDFASTQ1 $qual1;

	print SELECTEDFASTQ2 $header2;
	print SELECTEDFASTQ2 $sequence2, "\n";
	print SELECTEDFASTQ2 $qualheader2;
	print SELECTEDFASTQ2 $qual2;

} #end printFastq sub




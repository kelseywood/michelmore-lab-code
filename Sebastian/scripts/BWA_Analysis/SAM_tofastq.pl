#!/usr/bin/perl
use strict;
use warnings;

my $sam = $ARGV[0];

if ($sam eq "") {

	print "Error, not sam input was entered. Please enter a SAM file for analysis\n";
	die "Missing arguments\n";

} #end if

#open(FILTER, ">$sam.filter.sam");
open(FASTQUNSIN, ">$sam.unMapped.single.fastq");
open(FASTQUNPA1, ">$sam.unMapped.paired.1.fastq");
open(FASTQUNPA2, ">$sam.unMapped.paired.2.fastq");

open(FASTQMASIN, ">$sam.Mapped.single.fastq");
open(FASTQMAPA1, ">$sam.Mapped.paired.1.fastq");
open(FASTQMAPA2, ">$sam.Mapped.paired.2.fastq");


my %reads;

my $numMapReads = 0;

open(SAM, $sam);
while (my $samLine = <SAM>) {

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if ($alignment[2] ne "\*" && defined($alignment[9])) {

		++$reads{$alignment[0]};

		++$numMapReads;

#		print FILTER $samLine, "\n";

	} #end else map or unmap

} #end while
close(SAM);

my $numTotalReads = 0;

my $numUnPair = 0;
my $numUnSing = 0;
my $numMaPair = 0;
my $numMaSing = 0;

print STDERR "Finish loading reads names\n";

open(SAM, $sam);
while (my $samLine = <SAM>) {

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if( !($samLine =~ /^@/) ) {

		++$numTotalReads;

		if( defined($reads{$alignment[0]}) ) {

			if($reads{$alignment[0]} == 3 ) {

				print FASTQMAPA2 "@", $alignment[0],"/2\n";
				print FASTQMAPA2 $alignment[9],"\n";
				print FASTQMAPA2 "+\n";
				print FASTQMAPA2 $alignment[10],"\n";

				++$numMaPair;

			} elsif($reads{$alignment[0]} == 2 ) {

				print FASTQMAPA1 "@", $alignment[0],"/1\n";
				print FASTQMAPA1 $alignment[9],"\n";
				print FASTQMAPA1 "+\n";
				print FASTQMAPA1 $alignment[10],"\n";

				++$reads{$alignment[0]};

				++$numMaPair;

			} elsif($reads{$alignment[0]} == 1) {

				if($alignment[2] eq "\*") {

					print FASTQUNSIN "@", $alignment[0],"/1\n";
					print FASTQUNSIN $alignment[9],"\n";
					print FASTQUNSIN "+\n";
					print FASTQUNSIN $alignment[10],"\n";

					++$numUnSing

				} else {

					print FASTQMASIN "@", $alignment[0],"/1\n";
					print FASTQMASIN $alignment[9],"\n";
					print FASTQMASIN "+\n";
					print FASTQMASIN $alignment[10],"\n";

					++$numMaSing

				} #end else

			} #end if single unMapping

		} else {

			if($alignment[1] == 141 ) { 

				print FASTQUNPA2 "@", $alignment[0],"/2\n";
				print FASTQUNPA2 $alignment[9],"\n";
				print FASTQUNPA2 "+\n";
				print FASTQUNPA2 $alignment[10],"\n";

			} else {

				print FASTQUNPA1 "@", $alignment[0],"/1\n";
				print FASTQUNPA1 $alignment[9],"\n";
				print FASTQUNPA1 "+\n";
				print FASTQUNPA1 $alignment[10],"\n";

			} #end else

			++$numUnPair;

		} # end if both map

	} #end if its an alignment line

} #end while
close(SAM);

print "\nAlignment Stats\n\n";

print "Totals\n";
print "Total number of reads		", $numTotalReads, "\n";
print "Total number of map reads	", $numMapReads, "\n";
print "Total number of un-map reads	", $numTotalReads-$numMapReads, "\n\n";

print "Read Counts\n";
print "Total number of map reads in pairs		", $numMaPair, "\n";
print "Total number of map reads as broken pairs	", $numMaSing, "\n";
print "Total number of un-map reads in pairs		", $numUnPair, "\n";
print "Total number of un-map reads  as broken pairs	", $numUnSing, "\n\n";

print "Read Percents\n";
print "Percentage of map reads in pairs		", sprintf("%.2f", ($numMaPair/$numTotalReads)*100 ), "\n";
print "Percentage of map reads as broken pairs		", sprintf("%.2f", ($numMaSing/$numTotalReads)*100 ), "\n";
print "Percentage of un-map reads in pairs		", sprintf("%.2f", ($numUnPair/$numTotalReads)*100 ), "\n";
print "Percentage of un-map reads  as broken pairs	", sprintf("%.2f", ($numUnSing/$numTotalReads)*100 ), "\n\n";

exit;














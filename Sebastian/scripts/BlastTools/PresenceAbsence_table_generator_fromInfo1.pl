#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 4) {

	print "Usage: Lida_Families_Finder.pl <Identity Threeshold> <Size Ratio Threshold> <Output Prefix> <info1 file 1> <info1 file 2> ... <info1 file n>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $identityThreshold = shift(@ARGV);
my $sizeRatioThreshold = shift(@ARGV);
my $outputfile = shift(@ARGV);

my @info1Files = @ARGV;

print "Identity threshold : ", $identityThreshold, "\n";
print "Size ration threshold : ", $sizeRatioThreshold, "\n";

my %subjects;

foreach my $info1 (@info1Files) {

	open(INFO1, $info1) or die "Can't open $info1\n";

	while(my $hitLine = <INFO1>) {

		chomp($hitLine);

		my @fields = split("\t", $hitLine);

		if ( !defined($fields[5]) ) { next }

		# Get information out of the fields
		my $subject = $fields[1];

		my $identity = $fields[4];
		my $lengthAlignment = $fields[6];
		my $lengthQuery = $fields[8];
		my $lengthSubject = $fields[9];

		#Calculate the size ration
		my $sizeRatio;
		if($lengthSubject > $lengthQuery) { $sizeRatio = $lengthQuery/$lengthSubject }
		else { $sizeRatio = $lengthSubject/$lengthQuery }

		#Check if alignment pass the filters
		if( ($identity >= $identityThreshold) && ($sizeRatio >= $sizeRatioThreshold) ) {

			$subjects{$subject}{$info1} = 1;

#		} else {
#			if(!defined($hits{$query}) ) { print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n" }
#			print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n"

		} #end if alignment pass filters

	} #end of reading info1

} #end of while for multiple info1's

close(INFO1);

print "All hits have been gathered\n";
print "Number of sequences with hits ", scalar(keys %subjects), "\n\n";

open(PRESENCEABSENCE, ">$outputfile");

# print header for table
foreach my $info2 (@info2Files) { print PRESENCEABSENCE "\t", $info2}

print PRESENCEABSENCE "\n";

my @listSubjects = (keys %subjects);

foreach my $subject ( @listSubjects ) {

	print PRESENCEABSENCE $subject;

	foreach my $info1 (@info1Files) {

		if( defined($subjects{$subject}{$info1}) ) {
			print PRESENCEABSENCE "\t1";
		} else {
			print PRESENCEABSENCE "\t0";
		} #end if present or absent

	} #end readind subjects

	print PRESENCEABSENCE "\n";

} #end of printing subjects

exit;

### END OF SCRIPT


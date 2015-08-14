#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 4) {

	print "Usage: Lida_Families_Finder.pl <Identity Threeshold> <Perc coverage query> <Output Prefix> <info2 file 1> <info2 file 2> ... <info2 file n>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $identityThreshold = shift(@ARGV);
my $percQueryThreshold = shift(@ARGV);
my $outputfile = shift(@ARGV);

my @info2Files = @ARGV;

print "Identity threshold : ", $identityThreshold, "\n";
print "Percentage query aligned threshold : ", $percQueryThreshold, "\n";

my %queries;

foreach my $info2 (@info2Files) {

	open(INFO2, $info2) or die "Can't open $info2\n";

	while(my $hitLine = <INFO2>) {

		chomp($hitLine);

		my @fields = split("\t", $hitLine);

		if ( !defined($fields[5]) ) { next }

		# Get information out of the fields
		my $query = $fields[0];

		my $identity = $fields[4];
		my $lengthAlignment = $fields[6];
		my @lengthSequences = split("/", $fields[13]);
		my $lengthQuery = $lengthSequences[0];

		#Calculate the percentage alignment

		my $percAlignment = $lengthQuery/$lengthAlignment;

		#Check if alignment pass the filters
		if( ($identity >= $identityThreshold) && ($percAlignment >= $percQueryThreshold) ) {

			$queries{$query}{$info2} = 1;

#		} else {
#			if(!defined($hits{$query}) ) { print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n" }
#			print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n"

		} #end if alignment pass filters

	} #end of reading info2

} #end of while for multiple info2's

close(INFO2);

print "All hits have been gathered\n";
print "Number of sequences with hits ", scalar(keys %queries), "\n\n";

open(PRESENCEABSENCE, ">$outputfile");

# print header for table
foreach my $info2 (@info2Files) { print PRESENCEABSENCE "\t", $info2}

print PRESENCEABSENCE "\n";

my @listQueries= (keys %queries);

foreach my $query ( @listQueries ) {

	print PRESENCEABSENCE $query;

	foreach my $info2 (@info2Files) {

		if( defined($queries{$query}{$info2}) ) {
			print PRESENCEABSENCE "\t1";
		} else {
			print PRESENCEABSENCE "\t0";
		} #end if present or absent

	} #end readind subjects

	print PRESENCEABSENCE "\n";

} #end of printing subjects

exit;

### END OF SCRIPT


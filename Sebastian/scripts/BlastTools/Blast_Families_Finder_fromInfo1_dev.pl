#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 4) {

	print "Usage: Lida_Families_Finder.pl <info1 file> <Identity Threeshold> <Size Ratio Threshold> <Output Prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $info1 = $ARGV[0];
my $identityThreshold = $ARGV[1];
my $sizeRatioThreshold = $ARGV[2];
my $outputPrefix = $ARGV[3];

print "Alignments file : ", $info1, "\n";
print "Identity threshold : ", $identityThreshold, "\n";
print "Size ration threshold : ", $sizeRatioThreshold, "\n";

my %hits;
my %hitMatrix;

my %badHits;

open(INFO1, $info1);

#sopen(TEST, ">test.txt");

while(my $hitLine = <INFO1>) {

	chomp($hitLine);

	my @fields = split("\t", $hitLine);

	if ( !defined($fields[7]) ) { next }

	# Get information out of the fields
	my $query = $fields[0];
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
		$hits{$query}{$subject} = [$identity, $lengthSubject];
		$hitMatrix{$query}{$subject} = [$identity, $lengthSubject];
	} else {

		$badHits{$query} = "";
#		if(!defined($hits{$query}) ) { print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n" }
#		print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n"

	} #end if alignment pass filters

} #end while for reach info2 file

close(INFO1);

print "Matrix fully loaded\n";
print "Number of sequences with hits ", scalar(keys %hits), "\n\n";

# Creating output files
open(REPRESENTATIVE, ">$outputPrefix.id_$identityThreshold.sizeRatio_$sizeRatioThreshold.representative.txt");
open(FAMILIES, ">$outputPrefix.id_$identityThreshold.sizeRatio_$sizeRatioThreshold.families.txt");
open(GROUPED, ">$outputPrefix.id_$identityThreshold.sizeRatio_$sizeRatioThreshold.grouped.txt");

# Storage variables
my @representatives;
my %grouped;

my %clusters;
my %origin;

my $clusterCount = 1;

# Generating the families
foreach my $query (keys %hits) {

	# If grouped defined skip the iteration
	if( defined($grouped{$query}) ) { next } #end 

	my $clusterID = $clusterCount;

	# Get all the members of the family out of the hits connections
	my @familyMembers = check_relatives($query);

	push(@familyMembers, check_relatives($hitMatrix{$query}) );

	# Check if any of the nodes it's already in a cluster
	foreach my $line (@familyMembers) {
		if(defined($origin{$line})) { $clusterID = $origin{$line}; last } #end if
	} #end foreach

	# Getting the mean identity and size of each member 
	foreach my $line (@familyMembers) {

		# Initialize variables
		my $sumIdentity = 0;
		my $size = 0;
		my $numHits = 0;

		foreach my $hit (keys %{$hitMatrix{$line}}) {

			$sumIdentity += $hitMatrix{$line}{$hit}[0];
			$size = $hitMatrix{$line}{$hit}[1];
			++$numHits;
		
		} #end foreach

		my $meanIdentity;

		if($numHits == 0 ) {

			$meanIdentity = 0;

		} else {
		
			$meanIdentity = $sumIdentity/$numHits;

		} #end else

		$clusters{$clusterCount}{$line}{'meanIdentity'} = $meanIdentity;
		$clusters{$clusterCount}{$line}{'size'} = $size;

		$origin{$line} = $clusterCount;

	} #end foreach

	if ($clusterID == $clusterCount) { ++$clusterCount }

#	++$clusterCount;

} #end foreach queries

print "Finished constructing clusters\n";

print scalar(keys %clusters), " independent groups were found\n";

foreach my $cluster (keys %clusters) {

	# Hash to store the family members
	my %nodes = %{$clusters{$cluster}};

	my @nodesNames = (sort { $nodes{$b}->{'meanIdentity'} <=> $nodes{$a}->{'meanIdentity'} } keys %nodes );

	# Get the first element of the family and call it representative 
	my $representative = shift(@nodesNames);
	push (@representatives, $representative);
	$grouped{$representative} = "";

	# Print representative to all the files
	print REPRESENTATIVE $representative, "\t", scalar(@nodesNames) + 1, "\t", sprintf("%.2f", $nodes{$representative}{'meanIdentity'}), "\t", $nodes{$representative}{'size'}, "\n";
	print FAMILIES $representative, "\t", $representative, "\n";
	print GROUPED $representative, "\n";

	# Print the rest of the family to the files
	foreach my $node (@nodesNames) {

		delete($badHits{$node});

		print FAMILIES $representative, "\t", $node, "\n";

#		print FAMILIES "\t", $node, "\t", $hitMatrix{$node}{$representative}[0], "\n";
		print GROUPED $node, "\n";
		$grouped{$node} = "";

	} #end foreach
	
} #end foreach

print scalar(keys (%grouped) ), " sequences were grouped\n";
print scalar(keys (%badHits) ), " sequences are unique\n";

close(REPRESENTATIVE);
close(FAMILIES);

exit;

### END OF SCRIPT

sub check_relatives {

	my $value = shift;
	my @familyMembers;

	if( defined($hits{$value} ) ) {

		push(@familyMembers, $value);
		foreach my $hit (keys %{$hits{$value}}) {

			undef %{$hits{$value}};
			push(@familyMembers, $hit);

			if( defined($hits{$hit} ) ) { 

				push(@familyMembers, check_relatives($hit));

			} #end if hit defined
			
		} #end foreach

	} #end if

	return(@familyMembers);

} #end check relatives sub









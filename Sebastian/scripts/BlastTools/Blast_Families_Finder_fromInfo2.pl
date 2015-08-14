#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 4) {

	print "Usage: Lida_Families_Finder.pl <info2 file> <Identity Threeshold> <Size Ratio Threshold> <Percent Alignment Threeshold> <Output Prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $info2 = $ARGV[0];
my $identityThreshold = $ARGV[1];
my $sizeRatioThreshold = $ARGV[2];
my $percentAlignThreeshold = $ARGV[3];
my $outputPrefix = $ARGV[4];

print "Alignments file : ", $info2, "\n";
print "Identity threshold : ", $identityThreshold, "\n";
print "Size ration threshold : ", $sizeRatioThreshold, "\n";
print "Percent aligment threshold : ", $percentAlignThreeshold, "\n\n";

my %hits;
my %hitMatrix;

open(INFO2, $info2);

#sopen(TEST, ">test.txt");

while(my $hitLine = <INFO2>) {

	chomp($hitLine);

	my @fields = split("\t", $hitLine);

	if ( !defined($fields[5]) ) { next }

	# Get information out of the fields
	my $query = $fields[0];
	my $subject = $fields[1];

	my $identity = $fields[4];
	my $lengthAlignment = $fields[6];
	my @lengthSequences = split("/", $fields[12]);
	my $gapStat = $fields[13];

	my $percentAlign = ($lengthSequences[0]/$lengthAlignment) * 100 ;

	#Calculate the size ration
	my $sizeRatio;
	if($lengthSequences[0] > $lengthSequences[1]) { $sizeRatio = $lengthSequences[1]/$lengthSequences[0] }
	else { $sizeRatio = $lengthSequences[0]/$lengthSequences[1] }

	#Check if alignment pass the filters
	if( ($identity >= $identityThreshold) && ($sizeRatio >= $sizeRatioThreshold) && ($percentAlign >= $percentAlignThreeshold) && ($gapStat eq "NO_GAPS") ) {
		$hits{$query}{$subject} = [$identity, $lengthSequences[0]];
		$hitMatrix{$query}{$subject} = [$identity, $lengthSequences[0]];
	} else {

#		if(!defined($hits{$query}) ) { print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n" }
#		print TEST $query, "\t", $subject, "\t", $identity, "\t", $sizeRatio,"\t", $percentAlign, "\t", $gapStat, "\n"

	} #end if alignment pass filters

} #end while for reach info2 file

close(INFO2);

print "Matrix fully loaded\n";
print "Number of sequences with hits ", scalar(keys %hits), "\n\n";

# Creating output files
open(REPRESENTATIVE, ">$outputPrefix.id_$identityThreshold.alg_$percentAlignThreeshold.representative.txt");
open(FAMILIES, ">$outputPrefix.id_$identityThreshold.alg_$percentAlignThreeshold.families.txt");
open(GROUPED, ">$outputPrefix.id_$identityThreshold.alg_$percentAlignThreeshold.grouped.txt");

# Storage variables
my @representatives;
my %grouped;

# Generating the families
foreach my $query (keys %hits) {

	# If grouped defined skip the iteration
	if( defined($grouped{$query}) ) { next } #end 

	# Get all the members of the family out of the hits connections
	my @familyMembers = check_relatives($query);

	# Hash to store the family members
	my %nodes;

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
		
		my $meanIdentity = $sumIdentity/$numHits;
		$nodes{$line}{'meanIdentity'} = $meanIdentity;
		$nodes{$line}{'size'} = $size;

	} #end foreach

	my @nodesNames = (sort { $nodes{$b}->{'meanIdentity'} <=> $nodes{$a}->{'meanIdentity'} } keys %nodes );

	# Get the first element of the family and call it representative 
	my $representative = shift(@nodesNames);
	push (@representatives, $representative);
	$grouped{$representative} = "";

	# Print representative to all the files
	print REPRESENTATIVE $representative, "\t", scalar(@nodesNames) + 1, "\t", sprintf("%.2f", $nodes{$representative}{'meanIdentity'}), "\t", $nodes{$representative}{'size'}, "\n";
#	print FAMILIES $representative;
	print GROUPED $representative, "\n";

	# Print the rest of the family to the files
	foreach my $node (@nodesNames) {

		print FAMILIES $representative, "\t", $node, "\n";

#		print FAMILIES "\t", $node, "\t", $hitMatrix{$node}{$representative}[0], "\n";
		print GROUPED $node, "\n";
		$grouped{$node} = "";

	} #end foreach
	
} #end foreach

print scalar(@representatives), " independent groups were found\n";
print "Out of ", scalar(keys (%grouped) ), " sequences \n\n";

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









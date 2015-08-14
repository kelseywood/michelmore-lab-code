#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 3) {

	print "Usage: Lida_Families_Finder.pl <Evidence Overlaping file> <Sources File> <Output Prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $overlapping = $ARGV[0];
my $sources = $ARGV[1];
my $outputPrefix = $ARGV[2];

print "Overlapping file : ", $overlapping, "\n";
print "Sources file : ", $sources, "\n";



# open source related variables
my %geneModelsSources;

open(SOURCES, $sources);

print "\nReading sources and types from inputed file.\n";

#Load sources information into hash
while(my $source = <SOURCES>) {

	chomp($source);

	my @valuePairs = split("\t", $source);

	if(!(defined($valuePairs[1]))) {next}

	if($valuePairs[2] eq "GeneModels") {$geneModelsSources{$valuePairs[0]} = 0}

} #end while

my %perfectHits;
my %perfectHitMatrix;

open(OVERLAPPING, $overlapping);

while(my $hitLine = <OVERLAPPING>) {

	chomp($hitLine);

	my @fields = split("\t", $hitLine);

	#Check if alignment pass the filters
	if( defined($geneModelsSources{$fields[1]}) && ($fields[3] eq "match") && ($fields[4] == $fields[5]) && ($fields[5] == $fields[6]) ) {
		$perfectHits{$fields[0]}{$fields[2]} = $fields[4];
	} #end if alignment pass filters

} #end while for reach info2 file

close(OVERLAPPING);

my $numPerfectHits = scalar(keys %perfectHits);

print "Matrix fully loaded\n";

# Creating output files
open(GENEMODELS, ">$outputPrefix.geneModels.txt");

# Storage variables
my %redundants;
my %collapsed;

my $skipped;

# Generating the families
foreach my $geneModel (keys %perfectHits) {

	# If grouped defined skip the iteration
	if( defined($redundants{$geneModel}) ) { ++$skipped; next } #end 

	# Get all the members of the family out of the hits connections
	my %family;
	my %parents;

	my %familyMembers = check_relatives($geneModel, %family, %parents);

	my @uniqueMembers = sort keys %familyMembers;

	# Create the representative variable
	my $representative;

	# Verify if members belong to other groups
	foreach my $member (@uniqueMembers) { if(defined($redundants{$member})) { $representative = $redundants{$member}; last } }

	# If not, defined new group
	if(!defined($representative) ) { $representative = shift(@uniqueMembers) }

	undef %{$perfectHits{$representative}};

	$redundants{$representative} = $representative;

	# Print the rest of the family to the files
	foreach my $node (@uniqueMembers) {

		$redundants{$node} = $representative;

		$collapsed{$representative}{$node} = "";

		undef %{$perfectHits{$node}};


	} #end foreach
	
} #end foreach

foreach my $represativeNode (sort keys %collapsed) {

	my @nodes = keys %{$collapsed{$represativeNode}};

	print GENEMODELS $represativeNode, "\t", scalar(@nodes);
	foreach my $node (@nodes) {print GENEMODELS "\t", $node}

	print GENEMODELS "\n";

} #end foreach

print scalar(keys %collapsed), " independent groups were found\n";
print "With ", scalar(keys %redundants), " gene models\n";

close(GENEMODELS);

open(REDUNDANT, ">$outputPrefix.redundants.txt");

foreach my $redundant (sort keys %redundants) { if($redundant ne $redundants{$redundant}) {  print REDUNDANT $redundant, "\t", $redundants{$redundant}, "\n" } }


exit;

### END OF SCRIPT

sub check_relatives {

	my ($genemodel, %family, %parents) = @_;

#	my %family = %{$family_hashreference};

	if( defined($perfectHits{$genemodel}) && !defined($parents{$genemodel}) ) {

		$family{$genemodel} = "";

		my @hits = keys %{$perfectHits{$genemodel}};
		
		undef $perfectHits{$genemodel};

		foreach my $hit (@hits) {

			$family{$hit} = "";

			if( defined($perfectHits{$hit} ) ) { 

				%family = check_relatives($hit, %family, %parents);

			} #end if hit defined
			
		} #end foreach

	} #end if

	return(%family);


} #end check relatives sub


#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 3) {

	print "Usage: Lida_Families_Finder.pl <Evidence Overlaping file> <Sources File> <Scores File> <Output Prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $overlapping = $ARGV[0];
my $sources = $ARGV[1];
my $scores = $ARGV[2];
my $outputPrefix = $ARGV[3];

print "Overlapping file : ", $overlapping, "\n";
print "Sources file : ", $sources, "\n";
print "Scores file : ", $scores, "\n\n";



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

my %scores;

open(SCORES,  $scores);

while (my $scoreLine = <SCORES>) {

	chomp($scoreLine);

	my @score = split("\t", $scoreLine);

	$scores{$score[0]} = $score[1];

} #end of loading scores


my %Hits;
#my %perfectHitMatrix;

open(OVERLAPPING, $overlapping);

while(my $hitLine = <OVERLAPPING>) {

	chomp($hitLine);

	my @fields = split("\t", $hitLine);

	#Check if alignment pass the filters
	if( defined($geneModelsSources{$fields[1]}) ) {
		$Hits{$fields[0]}{$fields[2]} = $fields[4];
	} #end if alignment pass filters

} #end while for reach info2 file

close(OVERLAPPING);

my $numHits = scalar(keys %Hits);

print "Matrix fully loaded\n";

# Storage variables
my %redundants;
my %collapsed;

# Generating the families
foreach my $geneModel (keys %Hits) {

	# If grouped defined skip the iteration
	if( defined($redundants{$geneModel}) ) { next } #end 

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

	undef %{$Hits{$representative}};

	$redundants{$representative} = $representative;

	# Print the rest of the family to the files
	foreach my $node (@uniqueMembers) {

		$redundants{$node} = $representative;

		$collapsed{$representative}{$node} = "";

		undef %{$Hits{$node}};


	} #end foreach
	
} #end foreach

# Creating output files
open(REPRESENTATIVE, ">$outputPrefix.representative.txt");

foreach my $represativeNode (sort keys %collapsed) {

	my @nodes = keys %{$collapsed{$represativeNode}};

	print REPRESENTATIVE $represativeNode, "\t", scalar(@nodes);
	foreach my $node (@nodes) {print REPRESENTATIVE "\t", $node}

	print REPRESENTATIVE "\n";

} #end foreach

print scalar(keys %collapsed), " independent groups were found\n";
print "With ", scalar(keys %redundants), " gene models\n";

close(REPRESENTATIVE);

open(REDUNDANT, ">$outputPrefix.locis.txt");

foreach my $redundant (sort keys %redundants) { print REDUNDANT $redundant, "\t", $redundants{$redundant}, "\n" }

close(REDUNDANT);

open(UNGROUPED, ">$outputPrefix.ungrouped.txt");

foreach my $ungrouped (sort keys %scores) { if (!defined($redundants{$ungrouped})) { print UNGROUPED $ungrouped, "\n" } }

close(UNGROUPED);



exit;

### END OF SCRIPT

sub check_relatives {

	my ($genemodel, %family, %parents) = @_;

#	my %family = %{$family_hashreference};

	if( defined($Hits{$genemodel}) && !defined($parents{$genemodel}) ) {

		$family{$genemodel} = "";

		my @hits = keys %{$Hits{$genemodel}};
		
		undef $Hits{$genemodel};

		foreach my $hit (@hits) {

			$family{$hit} = "";

			if( defined($Hits{$hit} ) ) { 

				%family = check_relatives($hit, %family, %parents);

			} #end if hit defined
			
		} #end foreach

	} #end if

	return(%family);


} #end check relatives sub


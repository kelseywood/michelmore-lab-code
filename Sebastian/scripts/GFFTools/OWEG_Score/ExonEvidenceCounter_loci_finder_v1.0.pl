#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 4) {

	print "Usage: Lida_Families_Finder.pl <Evidence Overlaping file> <Sources File> <Scores file> <Output Prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $overlapping = $ARGV[0];
my $sources = $ARGV[1];
my $scores = $ARGV[2];
my $outputPrefix = $ARGV[3];

print "Overlapping file : ", $overlapping, "\n";
print "Sources file : ", $sources, "\n";
print "Scores file : ", $scores, "\n";



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

my %hits;
my %hitMatrix;

open(OVERLAPPING, $overlapping);

while(my $hitLine = <OVERLAPPING>) {

	chomp($hitLine);

	my @fields = split("\t", $hitLine);

	#Check if alignment pass the filters
	if( defined($geneModelsSources{$fields[1]}) ) {
		$hits{$fields[0]}{$fields[2]} = $fields[4];
		$hitMatrix{$fields[0]}{$fields[2]} = $scores{$fields[2]};
	} #end if alignment pass filters

} #end while for reach info2 file

close(OVERLAPPING);

print "Matrix fully loaded\n";

# Creating output files
open(LOCI, ">$outputPrefix.loci.txt");
open(GROUPS, ">$outputPrefix.groups.txt");

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
#		my $size = $scores{$line};

#		foreach my $hit (keys %{$hitMatrix{$line}}) {
#
#			$size = $hitMatrix{$line}{$hit};
#		
#		} #end foreach

		$nodes{$line}{'size'} = $scores{$line};

	} #end foreach

	my @nodesNames = (sort { $nodes{$b}->{'size'} <=> $nodes{$a}->{'size'} } keys %nodes );

	# Get the first element of the family and call it representative 
	my $representative = shift(@nodesNames);
	push (@representatives, $representative);
	$grouped{$representative} = "";

	# Print representative to all the files
	print LOCI $representative, "\t", scalar(@nodesNames) + 1, "\t", $nodes{$representative}{'size'}, "\n";
	print GROUPS $representative, "\t", $representative, "\n";

	# Print the rest of the family to the files
	foreach my $node (@nodesNames) {

		if(defined($grouped{$node})) {next}

		print GROUPS $representative, "\t", $node, "\n";

		$grouped{$node} = "";

	} #end foreach
	
} #end foreach

print scalar(@representatives), " independent groups were found\n";
print "Out of ", scalar(keys (%grouped) ), " sequences \n\n";

close(LOCI);
close(GROUPS);

open(UNGROUPED, ">$outputPrefix.ungrouped.txt");

foreach my $geneModel (keys %scores) {

	if(defined($grouped{$geneModel})) {next}

	print UNGROUPED $geneModel, "\n";

} #end foreach

close(UNGROUPED);

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









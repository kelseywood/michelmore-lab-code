#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

########################################################
#                                                      #
#                                                      #
#                 MST Mapper V1                        #
#                                                      #
#             Sebastian Reyes-Chin-Wo                  #
#                 Genome Center                        #
#        University of California at Davis             #
#                                                      #
#                                                      #
########################################################

# This script was designed to generate a genetic map from a loc file. The Loci file can have headers and annotation commented with ";" or "#".
# Also the MSTMap command line was build to map one single LG, it won't split data into indepedent LG and tails would be cut of.
#

#### Generating Global Variables 

# Running variables
my $loci_file;
my $tail_threshold = 20;
my $lg;

my $missingCharacter;
my $missingThreshold;

my $DRiterations = 0;

my $outDir;

# This variable should be updated to the location of the MadMapper script
my $madMapperPath = "~/scripts/Mapping_Tools/";
my $ChkMtxPath = "~/scripts/Mapping_Tools/";

# MSTMap variables
my $population_name;

my $population_type = "RIL6";
my $distance_function = "kosambi";
my $cut_off_p_value = 1;
my $no_map_dist = 15;
my $no_map_size = 2;
my $missing_threshold = 0.25;
my $estimation_before_clustering = "no";
my $detect_bad_data = "yes";
my $objective_function = "ML";


# Retrieve values from user input
GetOptions (

'loci_file=s' => \$loci_file,

'tail_threshold=i' => \$tail_threshold,
'double_recombinant_iterations=i' => \$DRiterations,

'missing_data_threshold=f' => \$missingThreshold,
'missing_character=s' => \$missingCharacter,

'lg=s' => \$lg,
'outDir=s' => \$outDir,

'population_name=s' => \$population_name,

'population_type=s' => \$population_type,
'distance_function=s' => \$distance_function,
'no_map_dist=i' => \$no_map_dist,
'no_map_size=i' => \$no_map_size,
'missing_threshold=i' => \$missing_threshold,
'estimation_before_clustering=s' => \$estimation_before_clustering,
'detect_bad_data=s' => \$detect_bad_data,
'objective_function=s' => \$objective_function

);

# Check if variables were defined
if (!defined($population_name) || !defined($loci_file) || !defined($lg) || !defined($outDir) || !defined($missingThreshold) || !defined($missingCharacter) ) {

	print "Please provided the required argumentsn";

	die "Missing arguments";

}


if($DRiterations > 3) {print STDERR  "Max allow of iterations of double recombinant masking it's three\n"; $DRiterations = 3}


#### Filter high rate missing data markers

open(LOCI, $loci_file) or die "can't opent $loci_file\n";
open(CLEANLOCI, ">$loci_file.missings$missingThreshold.loc") or die "can't opent $loci_file\n";

while(my $loci = <LOCI>) {
	if( $loci =~ /^;/ || $loci =~ /^#/) {print CLEANLOCI $loci}

	chomp($loci);

	my @genotypes = split("\t", $loci);
	my $markerName = shift(@genotypes);

	my $genotype = join("\t", @genotypes);

	my $missing = @{[$genotype =~ /$missingCharacter/g]};

	if( ( $missing / scalar(@genotypes) ) <= $missingThreshold) {

		print CLEANLOCI $loci, "\n";

	} #end if less missing

} #end getting numbers
close(LOCI);

#### Reformat Loci file

print STDERR "Generating MSTMap input file\n\n";

#open new file for MSTInput
open(MSTINPUT, ">$loci_file.MSTInput.loc") or die "Can't open $loci_file.MSTInput.loc\n";

# Create counts
my $nLoci = 0;
my $nInd;

# Calculating number of loci and individuals
open(LOCI, $loci_file) or die "can't opent $loci_file\n";
while(my $loci = <LOCI>) {
	if( $loci =~ /^;/ || $loci =~ /^#/) {next}
	my @fields = split("\t", $loci);
	$nInd = scalar(@fields);
	++$nLoci;
} #end getting numbers
close(LOCI);
--$nInd;

# Print the header sections to the file
print MSTINPUT "population_type $population_type\n";
print MSTINPUT " population_name $population_name\n";
print MSTINPUT " distance_function $distance_function\n";
print MSTINPUT " cut_off_p_value $cut_off_p_value\n";
print MSTINPUT " no_map_dist $no_map_dist\n";
print MSTINPUT " no_map_size $no_map_size\n";
print MSTINPUT " missing_threshold $missing_threshold\n";
print MSTINPUT " estimation_before_clustering $estimation_before_clustering\n";
print MSTINPUT " detect_bad_data $detect_bad_data\n";
print MSTINPUT " objective_function $objective_function\n";
print MSTINPUT " number_of_loci $nLoci\n";
print MSTINPUT " number_of_individual $nInd\n\n";

print MSTINPUT " locus_name";
for (my $i = 1; $i <= $nInd; ++$i) { print MSTINPUT "\ti", $i }
print MSTINPUT "\n";

# print loci information
open(LOCI, $loci_file) or die "can't opent $loci_file\n";
while(my $loci = <LOCI>) {
	if( $loci =~ /^;/ || $loci =~ /^#/) {next}

	chomp($loci);

	my @genotypes = split("\t", $loci);
	my $markerName = shift(@genotypes);

	my $genotype = join("\t", @genotypes);

	$genotype =~ s/$missingCharacter/\-/g;

	print MSTINPUT $markerName, "\t", $genotype, "\n";
} #end getting numbers
close(LOCI);



#### Run MSTMap

print STDERR "Running MSTMap\n\n";

if(system("MSTMap $loci_file.MSTInput.loc $loci_file.MSTOutput > $loci_file.MSTOutput.log 2> $loci_file.MSTOutput.error.log")) {
	die "Can't run MSTMap on $loci_file.MSTInput.loc\n";
} #end running MSTmap




#### Reformat MSTMap output

print STDERR "Extracting map information from MSTMap output\n\n";

open(MSTOUTPUT, "$loci_file.MSTOutput");
open(MAP, ">$loci_file.mstmap.map");

my $group = 1;

while(my $line = <MSTOUTPUT>) {
	if ($line =~ /;BEGINOFGROUP/) {
		while(my $mapLine = <MSTOUTPUT>) {

			if ($mapLine =~ /;ENDOFGROUP/) {last}
			chomp($mapLine);
			print MAP $lg, "\t", $mapLine, "\tframework\n";

		} #end while map
	} #end if its the beggining if the map
} #end while



#### Remove tails if asked (ready)

print STDERR "Removing tails\n\n";

if(system("Map_tails_remover.pl $loci_file.mstmap.map $tail_threshold > De-tailer.log") ) {
	die "Couldn't run de-tailer in \n";
} #end running de-tailer



#### Select markers and sort them (ready)

print STDERR "Selecting and sorting markers from loci file \n\n";

if(system("Mapped_Markers_selecter.pl $loci_file $loci_file.mstmap.map.withoutTails.map > MapMarker_selecter.log") ) {
	die "Couldn't run marker selecter and sorter\n";
} #end running de-tailer



### Get stats

print STDERR "Calculating general stats of the map\n\n";

# Get number of double recombinants
if(system("Get_mapStats.pl  $loci_file.mstmap.map.withoutTails.map.loc $loci_file.mstmap.map.withoutTails.map > GetStats.log") ) {
	die "Couldn't run get map stats\n";
} #end of double recombinant counter

# Get number of double recombinants
if(system("Double_recombinants_counter.pl $loci_file.mstmap.map.withoutTails.map.loc > Dbl_Recom_counter.log") ) {
	die "Couldn't run marker double recombinant counter\n";
} #end of double recombinant counter



### Run CheckMatrix

print STDERR "Running MadMapper\n\n";

#Running MadMapper

# Preparing MadMapper input
open(SORTEDLOCIFILE, "$loci_file.mstmap.map.withoutTails.map.loc");
open(MADMAPPERINPUT, ">$loci_file.MadMapperInput.loc");

print MADMAPPERINPUT ";";

for(my $i = 1; $i <= $nInd; ++$i) { print MADMAPPERINPUT "\t", $i}

print MADMAPPERINPUT "\n";

while(my $locLine = <SORTEDLOCIFILE>) {

	if( $locLine =~ /^;/ || $locLine =~ /^#/) {next}

	chomp($locLine);
	my @genotypes = split("\t", $locLine);
	my $markerName = shift(@genotypes);

	my $genotype = join("\t", @genotypes);

	$genotype =~ s/$missingCharacter/\-/g;

	$genotype =~ tr/ab/AB/;

	print MADMAPPERINPUT $markerName, "\t", $genotype, "\n";

} # End reformatting loci file

#if(system("python $madMapperPath/Python_MadMapper_V248_RECBIT_016NR.py $loci_file.MadMapperInput.loc MadMapper.temp 0.2 100 25 X 0.33 50 NOTRIO 3 > MadMapper.log 2> MadMapper.error.log")) {
#	die "Couldn't run MadMapper\n";
#} #end MadMapper

#system("mv MadMapper.temp.pairs_all $loci_file.mstmap.map.withoutTails.map.MadMapper.pairs_all");

#system("rm MadMapper.temp*");

#Running ChkMtx

print STDERR "Generating heatplots\n\n";

my $popType;

if( $population_type =~ /(R|r)(I|i)(L|l)/ ) {$popType = "RIL"} else { $popType = "F2"}

if(system("python $ChkMtxPath/py_matrix_2D_V254_RECBIT_V090710.py $loci_file.mstmap.map.withoutTails.map.MadMapper.pairs_all $loci_file.mstmap.map.withoutTails.map $loci_file.mstmap.map.withoutTails.map.ChkMtx X Y $loci_file.MadMapperInput.loc REC NOGRAPH 0.9 SMALL $popType > ChkMtx.log 2> ChkMtx.error.log")) {
	die "Couldn't run ChkMtx\n";
} #end running ChkMtx


### Loop through masking stages

### Remove double recombinants (ready)
my $workingLociFile = "$loci_file.MadMapperInput.loc";

print $DRiterations, "\n";

for (my $i = 1; $i <= $DRiterations; ++$i) {

	print "Running iteration $i of double recombinant masker\n";

	if(system("Double_recombinants_masker.pl $workingLociFile") ) {
		die "Couldn't run double recombinants masker\n";
	} #end running de-tailer

} #end loop iterations

exit;







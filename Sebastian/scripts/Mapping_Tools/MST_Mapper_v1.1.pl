#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

########################################################
#                                                      #
#                                                      #
#                 MST Mapper V1.1                      #
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


##### UPDATES LOG

### Difference from V1
# Added functionality to split modules
# Major re-estructuration of the script

#### Generating Global Variables

# Running variables
my $input_loci;
my $list_file;
my $tail_threshold;
my $lg;

my $missingCharacter;
my $missingThreshold;
my $nInd;

my $DRiterations = 0;

my $outDir;

my $help;

# This variable should be updated to the location of the MadMapper script
my $madMapperPath = "/home/sreyesch/scripts/Mapping_Tools/";
my $ChkMtxPath = "/home/sreyesch/scripts/Mapping_Tools/";

# MSTMap variables
my $population_name;

my $population_type = "RIL6";
my $distance_function = "kosambi";
my $cut_off_p_value = 1;
my $no_map_dist = 15;
my $no_map_size = 2;
my $estimation_before_clustering = "no";
my $detect_bad_data = "yes";
my $objective_function = "ML";


# Retrieve values from user input
GetOptions (

'loci_file=s' => \$input_loci,
'list_file=s' => \$list_file,

'tail_threshold=i' => \$tail_threshold,
'double_recombinant_iterations=i' => \$DRiterations,

'missing_data_threshold=f' => \$missingThreshold,
'missing_data_character=s' => \$missingCharacter,
'number_individuals=i' => $nInd,

'lg=s' => \$lg,
'outDir=s' => \$outDir,

'help' => \$help,

'population_name=s' => \$population_name,

'population_type=s' => \$population_type,
'distance_function=s' => \$distance_function,
'no_map_dist=i' => \$no_map_dist,
'no_map_size=i' => \$no_map_size,
'estimation_before_clustering=s' => \$estimation_before_clustering,
'detect_bad_data=s' => \$detect_bad_data,
'objective_function=s' => \$objective_function

);

my @errors;

if(!defined($input_loci) && !defined($list_file) ) { push(@errors,"Missing input_loci or list_file\n")}
if(!defined($population_name)) { push(@errors,"Missing population_name\n")}
if(!defined($lg)) { push(@errors,"Missing lg\n")}
if(!defined($outDir)) { push(@errors,"Missing outDir\n")}
if(!defined($missingThreshold)) { push(@errors,"Missing missing_data_threshold\n")}
if(!defined($missingCharacter)) { push(@errors,"Missing missing_data_character\n")}


# Check if variables were defined
if (defined($help)) {

	extendedHelp();	

	exit;	

} elsif (defined($errors[0]) ) {

	printHelp();

	foreach my $indError (@errors) {print $indError}

	print "\n";

	die "Missing arguments\n";

} # end if help it's need it

#!( defined($input_loci) || defined($list_file) ) || 

if($DRiterations > 3) {print STDERR  "Max allow of iterations of double recombinant masking it's three\n"; $DRiterations = 3}

print STDERR "Determining input format\n";

my @lociFiles;
if( defined($input_loci) && !defined($list_file) ) {

	print STDERR "A single file was inputed, mapping will run without batch\n\n";

	push(@lociFiles, $input_loci);

} elsif ( !defined($input_loci) && defined($list_file) ) {

	open(LOCIFILES, $list_file) or die "Can't open $list_file\n\n";

	while(my $loci = <LOCIFILES>) {

		chomp($loci);
		my @fileLine = split("\t", $loci);
		push(@lociFiles, $fileLine[0]);

	} #end while loading loci files

	print STDERR "File list was provided, mapping will run in batch on ", scalar(@lociFiles), " files\n\n";

} else {

	print STDERR "Both input_loci and list_file were define, please provided only one of the two\n";

} # end determining input format

system("mkdir $outDir") if !(-d $outDir);

foreach my $curr_loci_file (@lociFiles) {

	my $loci_file = $curr_loci_file;

	$loci_file =~ s{.*/}{};

	my $workingLociFile = $curr_loci_file;
	my $workingMapFile;

	system("mkdir $outDir/$loci_file") if !(-d "$outDir/$loci_file");

	for (my $i = 0; $i <= $DRiterations; ++$i) {

		#### Iteration variables
	 
		my $workingDir;
		my $iterationNum = "iteration$i";

		# Create counts
		my $nLoci = 0;


		#### Setting and creating directories if necessary

		if($i == 0 ) { $workingDir = "$outDir/$loci_file" }
		else {
			$workingDir = "$outDir/$loci_file/DRiteration$i";
		        print "\nRunning double recombinant masker iteration $i\n\n";
		}

		system("mkdir $workingDir") if !(-d $workingDir);
		system("mkdir $workingDir/logs") if !(-d "$workingDir/logs");


		#### Mask double recombinants

		# Perform double recombinant masking if it's in an iteration

		if($i != 0 ) { 

		        if(system("Double_recombinants_masker.pl $workingLociFile $workingDir/$loci_file.$iterationNum.DblRecMasked.loc") ) {
		        	die "Couldn't run double recombinants masker\n";
		        } #end running de-tailer

		        $workingLociFile = "$workingDir/$loci_file.$iterationNum.DblRecMasked.loc";

		} #end masking double recombinants if it's iterating


		#### Filter high rate missing data markers and generating MST input file

		print STDERR "Filtering markers with high rate of missing datapoints, coverting missing characters and generating MST input file\n\n";

		#open new file for MSTInput
		open(MSTINPUT, ">$workingDir/$loci_file.$iterationNum.MSTInput.loc") or die "Can't open $workingDir/$loci_file.$iterationNum.MSTInput.loc\n";
		# Print the header sections to the file
		print MSTINPUT <<mstinput
population_type $population_type
 population_name $population_name
 distance_function $distance_function
 cut_off_p_value $cut_off_p_value
 no_map_dist $no_map_dist
 no_map_size $no_map_size
 missing_threshold $missingThreshold
 estimation_before_clustering $estimation_before_clustering
 detect_bad_data $detect_bad_data
 objective_function $objective_function
 number_of_loci nLoci
 number_of_individual nInd

mstinput
;
		print MSTINPUT " locus_name";
		for (my $i = 1; $i <= $nInd; ++$i) { print MSTINPUT "\ti", $i }
		print MSTINPUT "\n";

		# open working loci file and clean file
		open(LOCI, $workingLociFile) or die "can't opent $workingLociFile\n";
		open(CLEANLOCI, ">$workingDir/$loci_file.$iterationNum.missings$missingThreshold.loc") or die "can't opent $loci_file\n";

		while(my $loci = <LOCI>) {

			if( $loci =~ /^;/ || $loci =~ /^#/) {
			        print CLEANLOCI $loci;
			        next;
			} # end if header

			chomp($loci);

			my @genotypes = split("\t", $loci);
			my $markerName = shift(@genotypes);

			my $genotype = join("\t", @genotypes);

			my $missingCount = ($genotype =~ s/$missingCharacter/-/g);

			$genotype =~ tr/abh/ABH/;

			my $nCalls = scalar(@genotypes);

			if($nCalls < 1) {next}

			if( ( $missingCount / $nCalls ) <= $missingThreshold) {

			        print CLEANLOCI $markerName, "\t", $genotype, "\n";
				++$nLoci;

				$genotype =~ s/(H|h)/U/g;
				print MSTINPUT $markerName, "\t", $genotype, "\n";

			} #end if less missing

		} #end getting numbers
		close(LOCI);
		close(CLEANLOCI);

		$workingLociFile = "$workingDir/$loci_file.$iterationNum.missings$missingThreshold.loc";

		system("perl -p -i -e 's/nLoci/$nLoci/' $workingDir/$loci_file.$iterationNum.MSTInput.loc");
		system("perl -p -i -e 's/nInd/$nInd/' $workingDir/$loci_file.$iterationNum.MSTInput.loc");


		#### Run MSTMap

		print STDERR "Running MSTMap\n\n";

		if(system("MSTMap $workingDir/$loci_file.$iterationNum.MSTInput.loc $workingDir/$loci_file.$iterationNum.MSTOutput > $workingDir/logs/$loci_file.MSTOutput.log 2> $workingDir/logs/$loci_file.MSTOutput.error.log")) {
			die "Can't run MSTMap on $workingDir/$loci_file.$iterationNum.MSTInput.loc\n";
		} #end running MSTmap


		#### Reformat MSTMap output

		print STDERR "Extracting map information from MSTMap output\n\n";

		open(MSTOUTPUT, "$workingDir/$loci_file.$iterationNum.MSTOutput");
		open(MAP, ">$workingDir/$loci_file.$iterationNum.mstmap.map");

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

		$workingMapFile = "$workingDir/$loci_file.$iterationNum.mstmap.map";


		#### Remove tails if asked

		if(defined($tail_threshold)) {
			print STDERR "Removing tails\n\n";
			if(system("Map_tails_remover.pl $workingMapFile $tail_threshold $workingMapFile.withoutTails.map > $workingDir/logs/$loci_file.De-tailer.log") ) {
			        die "Couldn't run de-tailer in $workingMapFile\n";
			} #end running de-tailer

			$workingMapFile = "$workingMapFile.withoutTails.map";

		} #end if remove tails


		#### Select markers and sort them (ready)

		print STDERR "Selecting and sorting markers from loci file \n\n";

		if(system("Mapped_Markers_selecter.pl $workingLociFile $workingMapFile $workingMapFile.loc > $workingDir/logs/$loci_file.MapMarker_selecter.log") ) {
			die "Couldn't run marker selecter and sorter\n";
		} #end running de-tailer

		$workingLociFile = "$workingMapFile.loc";


		### Get stats

		print STDERR "Calculating general stats of the map\n\n";

		# Get number of double recombinants
		if(system("Get_mapStats.pl  $workingLociFile $workingMapFile $workingMapFile.stats > $workingDir/logs/$loci_file.GetStats.log") ) {
			die "Couldn't run get map stats\n";
		} #end of double recombinant counter


		### Run CheckMatrix

		print STDERR "Running MadMapper\n\n";

		#Running MadMapper

		# Preparing MadMapper input
		open(SORTEDLOCIFILE, "$workingLociFile");
		open(MADMAPPERINPUT, ">$workingLociFile.MadMapperInput.loc");

		print MADMAPPERINPUT ";";

		for(my $i = 1; $i <= $nInd; ++$i) { print MADMAPPERINPUT "\t", $i}

		print MADMAPPERINPUT "\n";

		while(my $locLine = <SORTEDLOCIFILE>) {

			if( $locLine =~ /^;/ || $locLine =~ /^#/) {next}

			chomp($locLine);
			my @genotypes = split("\t", $locLine);
			my $markerName = shift(@genotypes);

			my $genotype = join("\t", @genotypes);

			print MADMAPPERINPUT $markerName, "\t", $genotype, "\n";

		} # End reformatting loci file

		$workingLociFile = "$workingLociFile.MadMapperInput.loc";

		if(system("python $madMapperPath/Python_MadMapper_V248_RECBIT_016NR.py $workingLociFile $workingDir/MadMapper.temp 0.2 100 25 X 0.33 50 NOTRIO 3 > $workingDir/logs/$loci_file.MadMapper.log 2> $workingDir/logs/$loci_file.MadMapper.error.log")) {
			die "Couldn't run MadMapper\n";
		} #end MadMapper

		system("mv $workingDir/MadMapper.temp.pairs_all $workingLociFile.MadMapper.pairs_all");
		system("mv $workingDir/MadMapper.temp.x_tree_clust $workingLociFile.MadMapper.tree_clust");

		system("rm $workingDir/MadMapper.temp*");

		#Running ChkMtx

		print STDERR "Generating heatplots\n\n";

		my $popType;

		if( $population_type =~ /(R|r)(I|i)(L|l)/ ) {$popType = "RIL"} else { $popType = "F2"}

		if(system("python $ChkMtxPath/py_matrix_2D_V254_RECBIT_V090710.py $workingLociFile.MadMapper.pairs_all $workingMapFile $workingMapFile.ChkMtx X Y $workingLociFile REC NOGRAPH 0.9 SMALL $popType > $workingDir/logs/$loci_file.ChkMtx.log 2> $workingDir/logs/$loci_file.ChkMtx.error.log")) {
			die "Couldn't run ChkMtx\n";
		} #end running ChkMtx

	} #end loop iterations

	system("mv $workingLociFile $loci_file.final.loc");

	system("mv $workingMapFile $loci_file.final.map");

} #end of running mappings

exit;

### SUBS

sub printHelp{

print "Usage: MST_Mapper_v1.pl -loci_file <loc file> -population_name <pop identifier> -lg <LG> -outDir <output dir> -missing_data_threshold <float> -missing_data_character <character> ..options\n\n";

print "Require arguments:\n";
print "loci_file = filename : loc file with genotype information (all non-markers lines should be commented with ';' or '#'.\n";
print "population_name = string : population identifier to be use in MSTmap.\n";
print "lg = string : Linkage group which to assign the map.\n";
print "missing_data_threshold = float : ratio from 0 to 1 of which to set the maximun allow missing data points per marker.\n";
print "missing_data_character = character : character that represents missing datapoints in the loc file.\n\n";

print "Optional arguments:\n";
print "double_recombinant_iterations = int[0:3] : if requested run iteration of double recombinants masker (default = 0) (max = 3).\n";
print "tail_threshold = int : cut off value in cM to remove tails from the end of the map (default = 20).\n";
print "estimation_before_clustering = yes/no : activate missing data estimation before clustering\n\n";

print "For extended help use -help option\n";

} #end printHelp sub

sub extendedHelp{

print "Usage: MST_Mapper_v1.pl -loci_file <loc file> -population_name <pop identifier> -lg <LG> -outDir <output dir> -missing_data_threshold <float> -missing_data_character <character> ..options\n\n";

print "Require arguments:\n";
print "loci_file = filename : loc file with genotype information (all non-markers lines should be commented with ';' or '#'.\n";
print "population_name = string : population identifier to be use in MSTmap.\n";
print "lg = string : Linkage group which to assign the map.\n";
print "missing_data_threshold = float : ratio from 0 to 1 of which to set the maximun allow missing data points per marker.\n";
print "missing_data_character = character : character that represents missing datapoints in the loc file.\n\n";

print "Optional arguments:\n";
print "double_recombinant_iterations = int[0:3] : if requested run iteration of double recombinants masker (default = 0) (max = 3).\n";
print "tail_threshold = int : cut off value in cM to remove tails from the end of the map (default = 20).\n";
print "estimation_before_clustering = yes/no : activate missing data estimation before clustering\n\n";

print "Please add the next lines in your .bash_profile or .bashrc file\n";

print "export PATH=\$PATH:/home/sreyesch/MichelmoreBin/MSTMap\n";
print "export PATH=\$PATH:/home/sreyesch/scripts/Mapping_Tools\n";
print "export PATH=/home/sreyesch/MichelmoreBin/anaconda/bin:\$PATH\n\n";

print "OUTPUT\n";
print "Final files will be loci_file.loc and loci_file.map containing the information from the markers used in the map and the map itself\n";
print "To be extended ...\n";

} # end extendedHelp sub


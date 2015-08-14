#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

########################################################
#                                                      #
#                                                      #
#                 MST Mapper V2                        #
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

#### Differences to version 1.
## Add section to bin markers before mapping


#### Generating Global Variables

# Running variables
my $input_loci;
my $tail_threshold = 20;
my $lg;

my $missingCharacter;
my $missingThreshold;

my $DRiterations = 0;

my $binning_before_mapping = "no";

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

'tail_threshold=i' => \$tail_threshold,
'double_recombinant_iterations=i' => \$DRiterations,

'missing_data_threshold=f' => \$missingThreshold,
'missing_data_character=s' => \$missingCharacter,

'lg=s' => \$lg,
'outDir=s' => \$outDir,

'binning_before_mapping=s' => \$binning_before_mapping,

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

# Check if variables were defined
if (defined($help)) {

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

	print "binning_before_mapping = yes/no : bin the markers (with MadMapper) before mapping them\n\n";

	print "Please add the next lines in your .bash_profile or .bashrc file\n";

	print "export PATH=\$PATH:/home/sreyesch/MichelmoreBin/MSTMap\n";
	print "export PATH=\$PATH:/home/sreyesch/scripts/Mapping_Tools\n";
	print "export PATH=/home/sreyesch/MichelmoreBin/anaconda/bin:\$PATH\n\n";

	print "OUTPUT\n";
	print "Final files will be loci_file.loc and loci_file.map containing the information from the markers used in the map and the map itself\n";
	print "To be extended ...\n";

	exit;	

} elsif (!defined($population_name) || !defined($input_loci) || !defined($lg) || !defined($outDir) || !defined($missingThreshold) || !defined($missingCharacter) ) {

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
	print "estimation_before_clustering = yes/no : activate missing data estimation before clustering\n";
	print "binning_before_mapping = yes/no : bin the markers (with MadMapper) before mapping them\n\n";

	print "For extended help use -help option\n";

	die "Missing arguments\n";

} # end if help it's need it


if($DRiterations > 3) {print STDERR  "Max allow of iterations of double recombinant masking it's three\n"; $DRiterations = 3}

system("mkdir $outDir") if !(-d $outDir);
system("mkdir $outDir/logs") if !(-d "$outDir/logs");

my $loci_file = $input_loci;

$loci_file =~ s{.*/}{};

my $workingLociFile = $input_loci;

# Create counts
my $nLoci = 0;
my $nInd;

if($binning_before_mapping =~ m/(Y|y)(E|e)(|S|s)/) {

	print STDERR "Binnig markers\n\n";

	open(LOCI, $workingLociFile) or die "can't opent $loci_file\n";
	open(CLEANLOCI, ">$outDir/$loci_file.missings$missingThreshold.loc") or die "can't opent $loci_file\n";

	while(my $loci = <LOCI>) {
	        if( $loci =~ /^;/ || $loci =~ /^#/) {
	                print CLEANLOCI $loci;
	                next;
	        }

	        chomp($loci);

	        my @genotypes = split("\t", $loci);
	        my $markerName = shift(@genotypes);

	        my $genotype = join("\t", @genotypes);

	        my $missingCount = ($genotype =~ s/$missingCharacter/$missingCharacter/g);

	        $genotype =~ s/$missingCharacter/\-/g;

	        $genotype =~ tr/abh/ABH/;

	        $nInd = scalar(@genotypes);

		if(scalar(@genotypes) < 1) {next}

	        if( ( $missingCount / scalar(@genotypes) ) <= $missingThreshold) {

	                print CLEANLOCI $loci, "\n";
		        ++$nLoci;

	        } #end if less missing

	} #end getting numbers
	close(LOCI);
	close(CLEANLOCI);


	# Preparing MadMapper input
	open(CLEANLOCI, "$outDir/$loci_file.missings$missingThreshold.loc");
	open(MADMAPPERINPUT, ">$outDir/$loci_file.missings$missingThreshold.MadMapperInput.loc");

	print MADMAPPERINPUT ";";

	for(my $i = 1; $i <= $nInd; ++$i) { print MADMAPPERINPUT "\t", $i}

	print MADMAPPERINPUT "\n";

	while(my $locLine = <CLEANLOCI>) {

	        if( $locLine =~ /^;/ || $locLine =~ /^#/) {next}

	        chomp($locLine);
	        my @genotypes = split("\t", $locLine);
	        my $markerName = shift(@genotypes);

	        my $genotype = join("\t", @genotypes);

	        print MADMAPPERINPUT $markerName, "\t", $genotype, "\n";

	} # End reformatting loci file

	if(system("python $madMapperPath/Python_MadMapper_V248_RECBIT_016NR.py $outDir/$loci_file.missings$missingThreshold.loc $outDir/MadMapper.temp 0.2 100 25 X 0.33 50 NOTRIO 3 > $outDir/logs/$loci_file.MadMapper.log 2> $outDir/logs/$loci_file.MadMapper.error.log")) {
	        die "Couldn't run MadMapper\n";
	} #end MadMapper

	system("mv $outDir/MadMapper.temp.x_tree_clust  $outDir/$loci_file.MadMapper.tree_clust");

	system("rm $outDir/MadMapper.temp*");

	if(system("Representative_bin_genotype_extractor.pl $outDir/$loci_file.missings$missingThreshold.loc $outDir/$loci_file.MadMapper.tree_clust $outDir/$loci_file.MadMapperBins > $outDir/logs/$loci_file.RepBin.log 2> $outDir/logs/$loci_file.RepBin.error.log") ) {
		die "Could't run bin extractor\n";
	} #end if representative extractor

	$workingLociFile = "$outDir/$loci_file.MadMapperBins.representatives.loc";

} #end of binning if asked


my $lastMap;
my $lastLoc;

for (my $i = 0; $i <= $DRiterations; ++$i) {

	my $workingDir;

        my $iterationNum = "iteration$i";

	$nLoci = 0;

	if($i == 0 ) { $workingDir = $outDir }
        else {
        	$workingDir = "$outDir/DRiteration$i";
                print "\nRunning double recombinant masker iteration $i\n\n";
        }

        system("mkdir $workingDir") if !(-d $workingDir);
	system("mkdir $workingDir/logs") if !(-d "$workingDir/logs");

	#### Filter high rate missing data markers

	print STDERR "Filtering markers with high rate of missing datapoints and coverting missing characters\n\n";

	open(LOCI, $workingLociFile) or die "can't opent $loci_file\n";
	open(CLEANLOCI, ">$workingDir/$loci_file.$iterationNum.missings$missingThreshold.loc") or die "can't opent $loci_file\n";

	while(my $loci = <LOCI>) {
	        if( $loci =~ /^;/ || $loci =~ /^#/) {
	                print CLEANLOCI $loci;
	                next;
	        }

	        chomp($loci);

	        my @genotypes = split("\t", $loci);
	        my $markerName = shift(@genotypes);

	        my $genotype = join("\t", @genotypes);

	        my $missingCount = ($genotype =~ s/$missingCharacter/$missingCharacter/g);

	        $genotype =~ s/$missingCharacter/\-/g;

	        $genotype =~ tr/abh/ABH/;

	        $nInd = scalar(@genotypes);

		if(scalar(@genotypes) < 1) {next}

	        if( ( $missingCount / scalar(@genotypes) ) <= $missingThreshold) {

	                print CLEANLOCI $loci, "\n";
		        ++$nLoci;

	        } #end if less missing

	} #end getting numbers
	close(LOCI);
	close(CLEANLOCI);

	#### Reformat Loci file

	print STDERR "Generating MSTMap input file\n\n";

	#open new file for MSTInput
	open(MSTINPUT, ">$workingDir/$loci_file.$iterationNum.MSTInput.loc") or die "Can't open $workingDir/$loci_file.$iterationNum.MSTInput.loc\n";

	# Print the header sections to the file
	print MSTINPUT "population_type $population_type\n";
	print MSTINPUT " population_name $population_name\n";
	print MSTINPUT " distance_function $distance_function\n";
	print MSTINPUT " cut_off_p_value $cut_off_p_value\n";
	print MSTINPUT " no_map_dist $no_map_dist\n";
	print MSTINPUT " no_map_size $no_map_size\n";
	print MSTINPUT " missing_threshold $missingThreshold\n";
	print MSTINPUT " estimation_before_clustering $estimation_before_clustering\n";
	print MSTINPUT " detect_bad_data $detect_bad_data\n";
	print MSTINPUT " objective_function $objective_function\n";
	print MSTINPUT " number_of_loci $nLoci\n";
	print MSTINPUT " number_of_individual $nInd\n\n";

	print MSTINPUT " locus_name";
	for (my $i = 1; $i <= $nInd; ++$i) { print MSTINPUT "\ti", $i }
	print MSTINPUT "\n";

	# print loci information
	open(LOCI,  "$workingDir/$loci_file.$iterationNum.missings$missingThreshold.loc") or die "can't opent $workingDir/$loci_file\n";
	while(my $loci = <LOCI>) {
	        if( $loci =~ /^;/ || $loci =~ /^#/) {next}

	        chomp($loci);

	        my @genotypes = split("\t", $loci);
	        my $markerName = shift(@genotypes);

                my $genotype = join("\t", @genotypes);

	        $genotype =~ s/$missingCharacter/\-/g;

	        $genotype =~ s/(H|h)/U/g;

	        print MSTINPUT $markerName, "\t", $genotype, "\n";
	} #end getting numbers
	close(LOCI);



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



	#### Remove tails if asked

#        if($tail_threshold > 0) {
	        print STDERR "Removing tails\n\n";
	        if(system("Map_tails_remover.pl $workingDir/$loci_file.$iterationNum.mstmap.map $tail_threshold > $workingDir/logs/$loci_file.De-tailer.log") ) {
	                die "Couldn't run de-tailer in \n";
	        } #end running de-tailer
#	} #end if remove tails


	#### Select markers and sort them (ready)

	print STDERR "Selecting and sorting markers from loci file \n\n";

	if(system("Mapped_Markers_selecter.pl $workingDir/$loci_file.$iterationNum.missings$missingThreshold.loc $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map > $workingDir/logs/$loci_file.MapMarker_selecter.log") ) {
	        die "Couldn't run marker selecter and sorter\n";
	} #end running de-tailer



	### Get stats

	print STDERR "Calculating general stats of the map\n\n";

	# Get number of double recombinants
	if(system("Get_mapStats.pl  $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.loc $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map > $workingDir/logs/$loci_file.GetStats.log") ) {
	        die "Couldn't run get map stats\n";
	} #end of double recombinant counter


	### Run CheckMatrix

	print STDERR "Running MadMapper\n\n";

	#Running MadMapper

	# Preparing MadMapper input
	open(SORTEDLOCIFILE, "$workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.loc");
	open(MADMAPPERINPUT, ">$workingDir/$loci_file.$iterationNum.MadMapperInput.loc");

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

	if(system("python $madMapperPath/Python_MadMapper_V248_RECBIT_016NR.py $workingDir/$loci_file.$iterationNum.MadMapperInput.loc $workingDir/MadMapper.temp 0.2 100 25 X 0.33 50 NOTRIO 3 > $workingDir/logs/$loci_file.MadMapper.log 2> $workingDir/logs/$loci_file.MadMapper.error.log")) {
	        die "Couldn't run MadMapper\n";
	} #end MadMapper

	system("mv $workingDir/MadMapper.temp.pairs_all $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.MadMapper.pairs_all");
	system("mv $workingDir/MadMapper.temp.x_tree_clust  $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.MadMapper.tree_clust");

	system("rm $workingDir/MadMapper.temp*");

	#Running ChkMtx

	print STDERR "Generating heatplots\n\n";

	my $popType;

	if( $population_type =~ /(R|r)(I|i)(L|l)/ ) {$popType = "RIL"} else { $popType = "F2"}

	if(system("python $ChkMtxPath/py_matrix_2D_V254_RECBIT_V090710.py $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.MadMapper.pairs_all $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.ChkMtx X Y $workingDir/$loci_file.$iterationNum.MadMapperInput.loc REC NOGRAPH 0.9 SMALL $popType > $workingDir/logs/$loci_file.ChkMtx.log 2> $workingDir/logs/$loci_file.ChkMtx.error.log")) {
	        die "Couldn't run ChkMtx\n";
	} #end running ChkMtx


	### Loop through masking stages

	### Remove double recombinants (ready)

	if(system("Double_recombinants_masker.pl $workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.loc") ) {
		die "Couldn't run double recombinants masker\n";
	} #end running de-tailer

        $workingLociFile = "$workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map.loc.clean.loc";
	$lastMap = "$workingDir/$loci_file.$iterationNum.mstmap.map.withoutTails.map";
	$lastLoc = "$workingDir/$loci_file.$iterationNum.MadMapperInput.loc";

} #end loop iterations

system("ln -s $lastLoc $outDir/$loci_file.final.loc");

system("ln -s $lastMap $outDir/$loci_file.final.map");

if(system("BinMapper_filler.pl $outDir/$loci_file.MadMapperBins.bins.txt $outDir/$loci_file.final.map > $outDir/logs/$loci_file.BinMapper_filler.log 2> $outDir/logs/$loci_file.BinMapper_filler.error.log") ) {
	die "Couldn't run BinMapper\n";
} #end

exit;



#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_l $opt_x $opt_o $opt_a $opt_D $opt_f $opt_r);
getopts('hl:x:o:a:D:f:r');

########################################################
########################################################
#                                                      #
#                                                      #
#           Exon Evidence Counter v1.3                 #
#                   08-09-2012                         #
#                                                      #
#             Sebastian Reyes-Chin-Wo                  #
#              sreyesch@ucdavis.edu                    #
#                                                      #
#                 Genome Center                        #
#         University of California at Davis            #
#                                                      #
#                                                      #
########################################################
########################################################


### Difference for version 1.0 ###
#
# Modified adjustment value
# Add section for handle annotated files
#
### Difference for version 1.2 ###
#
# Section to handle full evidence data (not only exon by exon)
# Split for annotated gene models
#
### Difference for version 1.21 ###
#
# Add recursive mode for looking gff files
#
### Difference from version 1.22 ###
#
# Eliminated strand check for repeats (repeats from both strand are considered for matching)
# Added source into type record

#### User ask for extended help documentation ####

if (defined($opt_h)) {

	#option h was provide, print full description

        usage();

        exit;

} elsif ( !(defined($opt_l)) || !(defined($opt_x)) || ( !(defined($opt_D)) && !(defined($opt_f))) || !(defined($opt_o)) ) {

	# At least one of the requried options is missing

	print "
	usage: .pl -l source-types.txt -x 5 -a string -o output -f/-D
	options:


          Required Options
	  -l <file>     Source/Types file
	  -x <int>	Bases allowance
	  -o <string>   Output directory to be created
          -f <file>     GFF file to be analyze (only one)
	  -D <string>   Directory with the files to analyze (bacth mode) (/home/user/directory/)

          Optional options
          -h            Print extended help
          -a            Flag for annotated gene models (has is writen in gff attributes line of the mRNA features)
          -w            Weights file for evidence scoring (this option will trigger the scoring section of the script)
	  -r		Activate recursive mode

        For a more detailed description and quick tips for tuning the script use option -h.
	";
	die"\n	Error: Missing Arguments\n\n";

} #end if options


####################################
#                                  #
#       Read in arguments          #
#                                  #
####################################

#Get inputed arguments
my $sourcesFile = $opt_l;
my $adj = $opt_x;
my $output = $opt_o;


#Verify if the output directory exist

if(-d "$output") {

	print STDERR "\n$output directory already exist, all the information will be overwritten.\nDo you still want to continue with the analysis?? {Y/N} : ";

	my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
	    print STDERR "Analysis will continue and information will be overwritten.\n";
	} else {
	    die "Analysis terminated by user, please input a different output directory.\n";
	}

} else {

	system("mkdir $output");

} #end else

#Print require arguments to screen
print "\nAnalysis options:\n";
print "Sources file = $sourcesFile\n";
print "Adjustment value = $adj\n";
print "Output directory = $output\n";

#Open and print arguments to log file
open (LOG, ">$output/log");
print LOG "Analysis options:\n";
print LOG "Sources file = $sourcesFile\n";
print LOG "Adjustment value = $adj\n";
print LOG "Output directory = $output\n";

#Print begining of command line to log file
print LOG "Command line : Exon_evidence_counter_v1.0.pl -l $sourcesFile -x $adj -o $output ";

#Define if it will be working for one file or a directory
my @files;

#will analyze a directory of files
if (defined($opt_D)) {

        print "Analyzing directory = $opt_D";

        print LOG "-D $opt_D";

	#Open working directory and get files to analyze
        open(LOGFILES, ">$output/Analyzed_Files.log");
        opendir(DIR, $opt_D) || die "Can't open directory $opt_D\n";

	if(defined($opt_r)) {
		print ": it will be analyzed recursively";
	        print LOG " -r";
	} #end if r on

	print "\n";

        while (my $input = readdir(DIR)) {

		# If recursive mode if active use check_file subroutine
		if(defined($opt_r)) {

			check_file($input, $opt_D);

		} else {

			#Verify extension of the file
			my $extension = substr($input,-3,);
			if( $extension eq "gff") {

				push(@files, "$opt_D/$input");
			        print LOGFILES $input, "\n";

			} #end if for gff files

		} #else recursively or not

   	} #end foreach directory

        if( !(defined($files[0])) ) {

        	print "Directory provided doesn't contain any .gff files.\n";
                die "Couldn't load any files to analyze";

        } #end if for empty files array

#will analyze only one file
} elsif (defined($opt_f)) {

	push(@files, $opt_f);

        print "Analyzing file = $opt_f\n\n";

        print LOG "-f $opt_f";

#None of the flags were present
} else {

       print "Please input a file (-f) or a directory (-D).\n\n";
       die "No input data";

} #end else

#Print to log file and screen if annotation flag was provided
if( defined($opt_a) ) {

	print LOG " -a $opt_a";
	print "Annotation flag: '$opt_a'\n";

} #end if

print LOG "\n\n";

my $numFiles = scalar(@files);


#### Prepare sources information ####

# open source related variables
my @tracks;
my @geneModelsSources;
my $inputValues;
my %datasets;
my %types;

open(SOURCES, $sourcesFile);

print LOG "Inputed evidence sources-types and categories\n";
print LOG "Sources\tTypes\tCategories\n";
print "\nReading sources and types from inputed file.\n";

#Load sources information into hash
while(my $source = <SOURCES>) {

	chomp($source);

	my @valuePairs = split("\t", $source);

	if(!defined($valuePairs[2])) {next}

	push(@tracks, [$valuePairs[0], $valuePairs[1], $valuePairs[2]]);
	++$inputValues;

	$datasets{$valuePairs[0]} = [$valuePairs[1], $valuePairs[2]];

	print LOG "$valuePairs[0]\t$valuePairs[1]\t$valuePairs[2]\n";

	$types{$valuePairs[2]} = 1;

	if($valuePairs[2] eq "GeneModels") { push(@geneModelsSources, $valuePairs[0])}

} #end while

my @sources = keys %datasets;

print LOG "\n\n";

#my @categories = keys(%datasets);

#check existing categories
if(!(exists($types{'GeneModels'}))) {
	die "Please define at least one track on the 'GeneModels' category on your input file.\n";

} #genemodel check

if ( !(exists($types{'Transcriptome'})) ) {

	print STDERR "Transcriptome category wasn't found.\nDo you want to continue with the analysis? {Y/N} : ";

	my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
		print STDERR "Analysis will continue without transcriptome information.\n";
	} else {
		die "Analysis terminated by user, missing 'Transcriptome' tag on sources file.\n";
	}
} #transcriptome date check

if ( !(exists($types{'Transposons'})) ) {
	print STDERR "Transposons category wasn't found, analysis will continue without it.\nDo you want to continue with the analysis? {Y/N} : ";

	my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
		print STDERR "Analysis will continue without transposable elements information.\n";
	} else {
		die "Analysis terminated by user, missing 'Transposons' tag on sources file.\n";
	}
} #transposon data check

if ( !(exists($types{'Repeats'})) ) {
	print STDERR "Repeats category wasn't found, analysis will continue without it.\nDo you want to continue with the analysis? {Y/N} : ";

	my $answer = <STDIN>;
	if ($answer =~ /^y(?:es)?$/i) {
		print STDERR "Analysis will continue without tandem repeats information.\n";
	} else {
		die "Analysis terminated by user, missing 'Repeats' tag on sources file.\n";
	}

} #repeat data check

print "Sources and Types loaded from $sourcesFile, ", scalar(@sources)," lines of evidence will be use for counting.\n\n";


####################################
#                                  #
#      Prepare output files        #
#                                  #
####################################


#### Error log ####

open(ERRORS, ">$output/error.log");


#### Files for annotated gene Models ####

#If annotation flag provided open annotated results files

my $annotationFlag;
my @outputAnnotatedFiles = ("$output/Annotated_Counters_1_ExonModels.txt",
	                    "$output/Annotated_Counters_2_ExonModels.txt",
	                    "$output/Annotated_Counters_3_ExonModels.txt",
	                    "$output/Annotated_Counters_4_ExonModels.txt",
	                    "$output/Annotated_Counters_more4_ExonModels.txt"
	                   );

if( defined($opt_a) ) {

	$annotationFlag = $opt_a;

	print "Annotation flag was provided, will generate indepedent output for annotated gene models with for models with $annotationFlag\n\n";

	foreach my $outputFile (@outputAnnotatedFiles) {

		open(ANNOTATED, ">", $outputFile);

		print ANNOTATED "#Data Definition Line\n";
		print ANNOTATED "#Source_file\tGene_Model\n";
		print ANNOTATED "#Organization Levels\n";
		print ANNOTATED "#First level: exons\n";
		print ANNOTATED "#Exon1\tExon2\tExon3\t..\tExon(n)\n";
		print ANNOTATED "#Second level: categories\n";

		foreach my $category (@sources) {

			print ANNOTATED $category,"_",$datasets{$category}[1], ":";

		} #end foreach

		print ANNOTATED "\n";

		print ANNOTATED "#Third level: status\n";
		print ANNOTATED "#matched_exons/partial_match/contradictory_exons/unsupported_exons/micro_exons/not_aligned_exons\n";

		close (ANNOTATED);

	} #end foreach annotated files

} #end if for annotation flag


#### Counters for all gene models per exon ####

my @outputEvidenceFiles = ("$output/Evidence_Counters_1_ExonModels.txt",
                           "$output/Evidence_Counters_2_ExonModels.txt",
                           "$output/Evidence_Counters_3_ExonModels.txt",
                           "$output/Evidence_Counters_4_ExonModels.txt",
                           "$output/Evidence_Counters_more4_ExonModels.txt"
                           );

foreach my $outputFile (@outputEvidenceFiles) {

	open(EXON, ">", $outputFile);

	print EXON "#Data Definition Line\n";
	print EXON "#Source_file\tGene_Model\n";
	print EXON "#Organization Levels\n";
	print EXON "#First level: exons\n";
	print EXON "#Exon1\tExon2\tExon3\t..\tExon(n)\n";
	print EXON "#Second level: categories\n";

	foreach my $category (@sources) { print EXON $category,"_",$datasets{$category}[1], ":" } #end foreach

	print EXON "\n";
	print EXON "#Third level: status\n";
	print EXON "#matched_exons/partial_match/contradictory_exons/unsupported_exons/micro_exons/not_aligned_exons\n";

	close(EXON);

} #end foreach


#### Gene model counters file ####

open(MODEL,">$output/Evidence_Counters_perGeneModel.txt");

print MODEL "GeneModel\tEvidence_Source\tEvidence_Name\tOverlaping\tExons_evidence\tExons_gene\tNum_Matched_Exons\tNum_Unsupported_Exons\tNum_Contradictory_Exons\n";

print MODEL "#Data Definition Line\n";
print MODEL "#Source_file\tGene_Model\n";
print MODEL "#Organization Levels\n";
print MODEL "#First level: exons\n";
print MODEL "#Exon1\tExon2\tExon3\t..\tExon(n)\n";
print MODEL "#Second level: categories\n";

foreach my $category (@sources) { print MODEL $category,"_",$datasets{$category}[1], ":" } #end foreach

print MODEL "\n";
print MODEL "#Third level: status\n";
print MODEL "#match/partialmatch/fuzzy/contradictory/\n";


#### Gene model against evidence overlapping file ####

open(OVERLAP,">$output/GeneModel_evidence_overlapping.txt");

print OVERLAP "GeneModel\tEvidence_Source\tEvidence_Name\tOverlaping\tExons_evidence\tExons_gene\tNum_Matched_Exons\tNum_Unsupported_Exons\tNum_Contradictory_Exons\n";


#### Directory to store all the comparison information ####

system("mkdir $output/Evidence_status") unless -d "$output/Evidence_status";


#Set global storage variables

my %all_gene_models;

#Only to check input information, enable this early exit function
#exit;

####################################
#                                  #
#     Main analysis pipeline       #
#                                  #
####################################

print STDERR "Starting the analysis.\n\n";

my $analizedfile = 0;

my $m = int $numFiles/100;

#When less than 100 files are analyzed m need to be corrected
if($m == 0) {$m = 1};

print STDERR "Percentage advance\n";

foreach my $file (@files) {

	#print advance status
	print STDERR int(100 * $analizedfile / $numFiles), " " if $analizedfile++ % $m == 0;

	my %evidence;

	my @filePath = split("/", $file);

	my $fileName = pop(@filePath);

	#open input file and new stats file
	if(!(-f $file)) {print "File $file not found\n"; next}
	open(INPUT, $file);
        open(STATS, ">$output/Evidence_status/$fileName.evidence_stats.txt");

	my $annotated = "FALSE";

	while (my $line = <INPUT>) {

		chomp($line);

		my @fields = split("\t", $line);

		my $score = 0;

		#if line contains proper number of fields
		if(defined($fields[8])) {

			#if it is a mRNA line (only use in case that annotation flag as been probided)
			if($fields[2] eq "mRNA" && defined($opt_a) ) {

				$annotated = "FALSE";

				my @attributes = split(";", $fields[8]);

				foreach my $attribute (@attributes) {

					if ($attribute eq $annotationFlag) {

						$annotated = "TRUE";

						last if $annotated eq "TRUE";

					} else {

						$annotated = "FALSE";

					} #end if annotated

				} #end foreach

			} #end if for mRNA features

# Loop to be eliminated and use the datasets hash
			#loop through inputed tracks (sources file)
#			for(my $i = 0; $i < scalar(@tracks); ++$i) {

				#check if source is need it and record score
				if(defined($datasets{$fields[1]}) ) {

					if($fields[5] ne ".") {

						$score = $fields[5];

					} #end if score present

				} #end if

				#check if source/type set is need it
				if( defined($datasets{$fields[1]}) && $datasets{$fields[1]}[0] eq $fields[2] ) {

					if($fields[6] ne ".") {

						$score = $fields[5];

					} #end if score present

					my @attributes = split(";", $fields[8]);

					my $parent;

					#Look for parent or name tags on attribute fields
					foreach my $attribute(@attributes) {

						if ($attribute =~ /^Parent=(\S+)/) {

							my @data = split ("=", $attribute);

							$parent = $data[1];

						} elsif ( ($attribute =~ /^Name=(\S+)/) && !(defined($parent)) ) {

							my @data = split ("=", $attribute);

							$parent = $data[1];

						} #end elsif parent

					} #end foreach

					#Determine stranding (change the order in which data should be recorded) and add data into the hash of hashes of arrays
                                        if ($fields[6] eq "+") {

						push(@{$evidence{$fields[1]}{$parent}}, [$fields[3], $fields[4], $score, $fields[6], $annotated]);

                                        } else {

						unshift(@{$evidence{$fields[1]}{$parent}}, [$fields[3], $fields[4], $score, $fields[6], $annotated]);

                                        } #end else

				} #end if

#			} #end for tracks

		} #end if for gff line exists

	} #end while files

#	print "Analyzing file $fileName.\n";

	#Loop through all gene model sources
	foreach my $source (@geneModelsSources) {

		#Loop through all the gene models in the hash
		foreach my $geneModel (keys (%{$evidence{$source}})) {

			#calculate the number of exons
			my $exons_gene = scalar(@{$evidence{$source}{$geneModel}});
			--$exons_gene;

			#load exons from main hash
			my @gene_exon_structure = @{$evidence{$source}{$geneModel}};

			#start local stat recording variable
		        my %gene_model_total_stats;

#	                print "\nAnalyzing $geneModel\n";

			my $done_flag = "FALSE";

			# Splitting by numbers of exons
			for(my $i = 0; $i <= 4; ++$i) {

				#check the numbers of exons
				if ($exons_gene eq $i || ( ($i eq 4) && ($exons_gene >= $i) ) ) {

					#Loop through the different categories
					foreach my $category (@sources) {

						#determine the adjustment value
						my $adjustment;
						if ($datasets{$category} eq "Predictors") {$adjustment = 0}
						elsif ($datasets{$category} eq "Transcriptome") {$adjustment = $adj}
						else {$adjustment = $adj+$adj;
						} #end if adjustment settings

						#identigy if working with repeats (different check subroutine)
			                	if( $datasets{$category}[1] ne "Repeats" && $datasets{$category}[1] ne "Transposons") {
							$gene_model_total_stats{$category} = { check_evidence($fileName, $category, $exons_gene, $geneModel, $adjustment, \@gene_exon_structure, \%evidence) };
			                        } else {
#							print "Going to check repeats for $geneModel\n";
							$gene_model_total_stats{$category} = { check_repeats($category, $exons_gene, $geneModel, $adjustment, \@gene_exon_structure, \%evidence) };
			                        } #end else

			                	#print "Category $category done \n";

					} #end foreach

#			                print "Summarizing evidence\n";

	#Modify function to use the datasets hash
					#find proper outfile and print data
			                my $outputFile = $outputEvidenceFiles[$i];
 		       	                summarize_evidence($fileName, $outputFile, $geneModel, $exons_gene, \@sources, \%datasets, \%gene_model_total_stats);

					#print to annotated file if annotated flag is TRUE
					if($gene_exon_structure[0][4] eq "TRUE") {

						#print "Output Send to annotated gene models\n";
				                my $outputFile = $outputAnnotatedFiles[$i];
		        	                summarize_evidence($fileName, $outputFile, $geneModel, $exons_gene, \%datasets, \%gene_model_total_stats);

					} #end if for annotated outputing

				} # if check number exons

			} #end for exons number

		} #end foreach genemodel

	} #end foreach source

#        foreach my $evidenceType (keys %evidence) {
#        	print $evidenceType, "\t", scalar (keys %{$evidence{$evidenceType}}), "\n";
#        } #end foreach evidencetype

} #end for (files)

print "\nAnalysis done\n\n";

exit;

###########################################################################################
###########################################################################################
#
#   END OF SCRIPT
#
###########################################################################################
###########################################################################################


####################################
#                                  #
#         Check file sub           #
#                                  #
####################################

#This sub is only use if recursive mode is active

sub check_file {

	my ($file, $path) = @_;

	#filter out . and ..
	if($file ne "." && $file ne "..") {

#		print "$file, $path\n";

		#Check if it is a directory
		if(-d "$path/$file") {

			#open subdir and make new path
			opendir(SUBDIR, "$path/$file");
			my $subpath = $path."/".$file;
			my @dirs;

			while(my $subinput = readdir(SUBDIR) ) {

				#avoid . and ..
				if($subinput ne "." && $subinput ne "..") {

					#if is a directory add to dir list or add to gff list if .gff extension
					if(-d "$subpath/$subinput") { push(@dirs,$subinput) }
					else {

						my $extension = substr($subinput,-3,);
						if( $extension eq "gff") {

							push(@files, "$subpath/$subinput");
							print LOGFILES "$subinput\t$subpath\n";

						} #end if for gff files

					} #end else is dir

				} #end if

			} #end while subdir

			#after analyzed all the files, analyzed the directories
			foreach my $dir (@dirs) {check_file($dir, $subpath)} #end foreach dirs

		#if not a directory evaluate extenssion
		} else {

			# if it is a gff add to file list
			my $extension = substr($file,-3,);
			if( $extension eq "gff") {

				push(@files, "$path/$file");
				print LOGFILES "$file\t$path\n";

			} #end if for gff files

		} #end else is dir

	} #end if is a real directory

} #end check_file sub


####################################
#                                  #
#       Check evidence sub         #
#                                  #
####################################

sub check_evidence {

	my ($fileName, $category, $exons_gene, $geneModel, $adjustment, $geneExons_ref, $evidence) = @_;

	my (@geneExons) = @$geneExons_ref;

        my %geneModel_stats;

	foreach my $transcript ( keys %{$$evidence{$category}}) {

		my $exons_transcript = scalar(@{$$evidence{$category}{$transcript}});
		--$exons_transcript;

		my @eviExons = @{$$evidence{$category}{$transcript}};

                #Check if evidence and gene are in the same strand or evidence don't have strand

                if($eviExons[0][3] eq $geneExons[0][3] || $eviExons[0][3] eq "." || $geneExons[0][3] eq ".") {

	                my $run = "TRUE";

	                #Check if evidence overlaps with gene model
	                if( ($geneExons[$exons_gene][1] < $eviExons[0][0]) || ($geneExons[0][0] > $eviExons[$exons_transcript][1]) ) {

				#If there isn't overlap skip
                        	#problems setting the oposite condition

                        } elsif ( ($geneModel ne $transcript) && ($run eq "TRUE") ) {

	                        #Exon locators
	                        my $exonEvi = 0;
	                        my $exonGene = 0;

	                        #Status definers
	                        my $within = "FALSE";

	                        #Exons stat recorders
	                        my %evi_exon_stats;
	                        my %gene_exon_stats;

	                        #Scan exons in evidence against exons in gene model
	                        while($exonEvi <= $exons_transcript) {

	                                if($exonGene > $exons_gene) {

	                                        #If we are in the end of gene and evidence, add to terminate the loop
	                                        if ($exonEvi >= $exons_transcript) {

	                                                last;

	                                        } #end if

	                                        $exonGene = $exons_gene;

	                                } #end if

	                                #Compare start and end of both exons and determine the status
	                                #Comparison is made using the gene model location has reference
	                                #Nomenclature is has follows: gene_evidence
	                                #To understand the meaning read like this, for example $start_start = "before", means start of the evidence its before the start of the gene model (for a particular exon)

	                                my $start_start = "";
	                                my $start_end = "";
	                                my $end_start = "";
	                                my $end_end = "";

	                                #Start of the gene against start of the evidence
	                                if( ($geneExons[$exonGene][0]-$adjustment) > $eviExons[$exonEvi][0]) { $start_start = "before"}

	                                elsif ( ($geneExons[$exonGene][0]+$adjustment) < $eviExons[$exonEvi][0]) { $start_start = "after"}

	                                else { $start_start = "match"

	                                } #end else for start-start comparison


	                                #Start of the gene against end of the evidence
	                                if( ($geneExons[$exonGene][0]-$adjustment) > $eviExons[$exonEvi][1]) { $start_end = "before"}

	                                elsif ( ($geneExons[$exonGene][0]+$adjustment) < $eviExons[$exonEvi][1]) { $start_end = "after"}

	                                else { $start_end = "match"

	                                } #end else for start-end comparison


	                                #End of the gene against start of the evidence
	                                if( ($geneExons[$exonGene][1]-$adjustment) > $eviExons[$exonEvi][0]) { $end_start = "before"}

	                                elsif ( ($geneExons[$exonGene][1]+$adjustment) < $eviExons[$exonEvi][0]) { $end_start = "after"}

	                                else { $end_start = "match"

	                                } #end else for end-start comparison


	                                #End of the gene against end of the evidence
	                                if( ($geneExons[$exonGene][1]-$adjustment) > $eviExons[$exonEvi][1]) { $end_end = "before"}

	                                elsif ( ($geneExons[$exonGene][1]+$adjustment) < $eviExons[$exonEvi][1]) {    $end_end = "after"}

	                                else { $end_end = "match"

	                                } #end else for end-end comparison


	                                #Check if we are within the gene model
	                                if(($geneExons[0][0]-$adjustment > $eviExons[$exonEvi][0]) && ($geneExons[$exons_gene][1]+$adjustment < $eviExons[$exonEvi][1]) ) {

	                                        $within = "TRUE";

	                                } #end if within gene model


	                                #Set_status for the exons

	                                #Check for micro gene exons
	                                if ( ($geneExons[$exonGene][1] - $geneExons[$exonGene][0]) < 20 ) {

	                                        #print "Micro exon\n";

	                                        $gene_exon_stats{$exonGene} = "Micro_exon";

						#Status
						if($exonGene eq $exons_gene) {++$exonEvi}

						#Walking
	                                        ++$exonGene;

	                                #Exon evidence before exon of gene
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "after") && ($end_end eq "after") ) {

						#Example
						#
						#Gene model        |---|
						#
						#Evidence   |---|

#	                                        print "1 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($within eq "TRUE") {$gene_exon_stats{$exonGene} = "Unsupported_exon"}

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi}
	                                        ++$exonGene;

	                                #Exon gene before exon of evidence
	                                } elsif ( ($start_start eq "before") && ($start_end eq "before") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence          |---|

#	                                        print "2 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($within eq "FALSE") {$evi_exon_stats{$exonEvi} = "Extra_exon_ends"}
	                                        else {$evi_exon_stats{$exonEvi} = "Extra_exon_within"}

						#Walking
	                                        ++$exonEvi;

	                                #Overlaping exons (evidence first)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model   |---|
						#
						#Evidence   |---|

#	                                        print "3 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($category eq "Transcriptome" && $exons_transcript eq 0 && $exonGene eq 0) {
	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match"
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory"
	                                        } #else

						#Walking
	                                        ++$exonEvi;

	                                #Overlaping exons (gene first)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence     |---|

#	                                        print "4 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($category eq "Transcriptome" && $exons_transcript eq 0 && $exonGene eq $exons_gene) {
	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #else

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi}
	                                        ++$exonGene;

	                                #Partial match (end matches) (evidence shorter)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence    |--|

#	                                        print "5 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($exonEvi eq 0) {
	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #else

						#Walking
	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Partial match (end matches) (gene shorter)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#
						#Gene model  |--|
						#
						#Evidence   |---|

#	                                        print "6 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($category eq "Transcriptome" && $exonGene eq 0 && $exonEvi eq 0) {
	                                                $evi_exon_stats{$exonEvi} = "Match";
	                                                $gene_exon_stats{$exonGene} = "Match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #else

						#Walking
	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Partial match (start matches) (evidence shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence   |--|

#	                                        print "7 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($exonEvi eq $exons_transcript) {
	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #else

						#Walking
	                                        ++$exonEvi;

	                                #Partial match (start matches) (gene shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#
						#Gene model |--|
						#
						#Evidence   |---|

#	                                        print "8 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($category eq "Transcriptome" && $exonGene eq $exons_gene && $exonEvi eq $exons_transcript) {
	                                                $evi_exon_stats{$exonEvi} = "Match";
	                                                $gene_exon_stats{$exonGene} = "Match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #else

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi} #end if
	                                        ++$exonGene;

	                                #Evidence exon within gene exon
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence    |-|

#	                                        print "9 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        if($exons_transcript eq 0) {
	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";
	                                        } else {
	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";
	                                        } #end else

						#Walking
	                                        ++$exonEvi;

	                                #Gene exon within evidence exon
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#
						#Gene model  |-|
						#
						#Evidence   |---|

#	                                        print "10 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        $evi_exon_stats{$exonEvi} = "Contradictory";
	                                        $gene_exon_stats{$exonGene} = "Contradictory";

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi} #end if
	                                        ++$exonGene;

	                                #Overlapping end-start
	                                } elsif ( ($start_start eq "before") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model     |---|
						#
						#Evidence   |---|

#	                                        print "11 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Walking
	                                        ++$exonEvi;

	                                #Overlapping start-end
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "after") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence       |---|

#	                                        print "12 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi} #end if
	                                        ++$exonGene;

	                                #
	                                } elsif ( ($start_start eq "match") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence   |

#	                                        print "13 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Walking
	                                        ++$exonEvi;

	                                #
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "match") ) {

						#Example
						#
						#Gene model     |
						#
						#Evidence   |---|

#	                                        print "14 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Walking
						if($exonGene eq $exons_gene) {++$exonEvi} #end if
	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Perfect match
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#
						#Gene model |---|
						#
						#Evidence   |---|

#	                                        print "Match $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

						#Status
	                                        $evi_exon_stats{$exonEvi} = "Match";
	                                        $gene_exon_stats{$exonGene} = "Match";

						#Walking
	                                        ++$exonGene;
	                                        ++$exonEvi;

					#Anything else
	                                } else {

#	                                        print "Whatever $exonEvi $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1] $eviExons[$exonEvi][0] $eviExons[$exonEvi][1]\n";
#	                                        print "$start_start $start_end $end_start $end_end\n";

						chomp($geneModel);

						#Status
	                                        print ERRORS "In $fileName none of the comparisons return TRUE for $transcript, exon $exonEvi, against $geneModel, exon $exonGene.\n";
						print ERRORS "        $start_start $start_end $end_start $end_end\n";

						#Walking
	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                } # end else check status

	                        } #end while evidence exons

				#Set status counters
	                        my $match_exons = 0;
	                        my $partial_match = 0;
	                        my $contradictory_exons = 0;
	                        my $unsupported_exons = 0;
	                        my $micro_exons = 0;

	                        print STATS $geneModel, "\t", $transcript;

				#Count exon stats
	                        for(my $i = 0; $i <= $exons_gene; ++$i) {
                                        if( defined($gene_exon_stats{$i}) ) {

	                                        print STATS "\t", $gene_exon_stats{$i};

                                                if($gene_exon_stats{$i} eq "Match") { ++$geneModel_stats{"exons"}{$i}{'matched_exons'}; ++$match_exons }
                                                elsif($gene_exon_stats{$i} eq "Partial_match") { ++$geneModel_stats{"exons"}{$i}{'partial_match'}; ++$match_exons }
                                                elsif($gene_exon_stats{$i} eq "Contradictory") { ++$geneModel_stats{"exons"}{$i}{'contradictory_exons'}; ++$contradictory_exons }
                                                elsif($gene_exon_stats{$i} eq "Unsupported_exon") { ++$geneModel_stats{"exons"}{$i}{'unsupported_exons'}; ++$unsupported_exons }
                                                elsif($gene_exon_stats{$i} eq "Micro_exon") { ++$geneModel_stats{"exons"}{$i}{'micro_exons'}; ++$micro_exons }

	                                } else {

	                                        print STATS "\t";
                                                ++$geneModel_stats{"exons"}{$i}{'not_aligned_exons'}

	                                } #end else for defined exon stat

	                        } #end foreach

				#Create gene model stats
				if ( ( $unsupported_exons == 0 ) && ( $contradictory_exons == 0 ) && ( $match_exons > 0 ) ) {

					print OVERLAP $geneModel,"\t",$category,"\t",$transcript,"\tmatch\t",($exons_transcript+1), "\t",($exons_gene+1), "\t", $match_exons, "\t", $unsupported_exons, "\t", $contradictory_exons,"\n";

					++$geneModel_stats{"total"}{"match"};

				} elsif ( ( $unsupported_exons > 0 ) && ( $unsupported_exons < $match_exons ) && ( $contradictory_exons == 0 ) ) {

					print OVERLAP $geneModel,"\t",$category,"\t",$transcript,"\tpartialMatch\t",($exons_transcript+1), "\t",($exons_gene+1), "\t",$match_exons,"\t", $unsupported_exons, "\t", $contradictory_exons,"\n";

					++$geneModel_stats{"total"}{"parcialMatch"};

				} elsif ( ( $unsupported_exons > 0 ) && ( $unsupported_exons > $match_exons ) && ( $contradictory_exons == 0 ) ) {

					print OVERLAP $geneModel,"\t",$category,"\t",$transcript,"\tfuzzy\t",($exons_transcript+1), "\t",($exons_gene+1), "\t",$match_exons, "\t", $unsupported_exons, "\t", $contradictory_exons,"\n";

					++$geneModel_stats{"total"}{"fuzzy"};

				} elsif ( ($contradictory_exons > 0) && ( $contradictory_exons < $match_exons ) ) {

					print OVERLAP $geneModel,"\t",$category,"\t",$transcript,"\tfuzzy\t",($exons_transcript+1), "\t",($exons_gene+1), "\t",$match_exons, "\t", $unsupported_exons, "\t", $contradictory_exons,"\n";

					++$geneModel_stats{"total"}{"fuzzy"};

				} elsif ( ($contradictory_exons > 0) && ( $contradictory_exons > $match_exons ) ) {

					print OVERLAP $geneModel,"\t",$category,"\t",$transcript,"\tcontradictory\t",($exons_transcript+1), "\t",($exons_gene+1), "\t",$match_exons, "\t", $unsupported_exons, "\t", $contradictory_exons,"\n";

					++$geneModel_stats{"total"}{"contradictory"};

				} #elsif

	                        print STATS "\n";

	                        #Exon stats counters
	                        my $extra_exons_ends = 0;
	                        my $extra_exons_within = 0;

	                } #end if for overlapping

                } #end if for strand

	} #end foreach transcripts

        return(%geneModel_stats);

} #end sub check_evidence


####################################
#                                  #
#       Check repeats sub          #
#                                  #
####################################

sub check_repeats {

	my ($category, $exons_gene, $geneModel, $adjustment, $geneExons_ref, $evidence) = @_;

	my (@geneExons) = @$geneExons_ref;

        my %geneModel_stats;

#        print "Checking repeats for $geneModel, ", keys %{$$evidence{$category}}, "\n";

	foreach my $repeat ( keys %{$$evidence{$category}}) {

#        	print ".";

		my @repeatExons = @{$$evidence{$category}{$repeat}};

                #Check if evidence and gene are in the same strand or evidence don't have strand
                if($repeatExons[0][3] eq $geneExons[0][3] || $repeatExons[0][3] eq ".") {

	                #Check if evidence overlaps with gene model
	                if( ($geneExons[$exons_gene][1] < $repeatExons[0][0]) || ($geneExons[0][0] > $repeatExons[0][1]) ) {

				#If there isn't overlap skip
                        	#problems setting the oposite condition

                        } else {

#	                        print $geneModel, "\t", $repeat, "\n";

	                        print STATS $geneModel, "\t", $repeat;

	                        #Exon locators
	                        my $exonGene = 0;
	                        my %gene_exon_stats;
				my $repeatOverlap = 0;

                                my $start_end = "";
                                my $end_start = "";

	                        #Scan exons in evidence against exons in gene model
	                        while($exonGene <= $exons_gene) {

	                                #Start of the gene against end of the repeat
	                                if( ($geneExons[$exonGene][0]-$adjustment) > $repeatExons[0][1]) { $start_end = "before"

	                                } elsif ( ($geneExons[$exonGene][0]+$adjustment) < $repeatExons[0][1]) {    $start_end = "after"

	                                } else { $start_end = "match"

	                                } #end else for start-end comparison


	                                #End of the gene against start of the repeat
	                                if( ($geneExons[$exonGene][1]-$adjustment) > $repeatExons[0][0]) { $end_start = "before"

	                                } elsif ( ($geneExons[$exonGene][1]+$adjustment) < $repeatExons[0][0]) {    $end_start = "after"

	                                } else { $end_start = "match"

	                                } #end else for end-start comparison

                                        #if there is any overlap with the repeat
                                        if( ($start_end eq "after" && $end_start eq "before") || ($start_end eq "before" && $end_start eq "after") || $start_end eq "match" || $end_start eq "match" ) {

                                        	print STATS "\t", $category;

                                        	$gene_exon_stats{$exonGene} = $category;

                                                ++$geneModel_stats{$exonGene};

						$repeatOverlap = 1;

                                        } else {

	                                        print STATS "\t";

	                                } #end else for defined exon stat

                                        ++$exonGene;

                                } #end while loop for exons

				#Add counter to total genemodel counters
				if($repeatOverlap == 1) {

					print OVERLAP $geneModel,"\t",$category,"\t",$repeat,"\tOverlap\t1\t",($exons_gene+1), "\t1\n";

					++$geneModel_stats{"total"};

				} #end if repeats were found


	                        print STATS "\n";

                        }#end if overlapping

                } #end if for strand

        } #end foreach repeat

        return(%geneModel_stats);

} #end sub check_repeats


####################################
#                                  #
#     Summarize evidence sub       #
#                                  #
####################################

sub summarize_evidence {

	my ($fileName,$outputFile, $geneModel, $exons_gene, $sources_reference, $datasets_reference, $gene_model_total_stats_reference) = @_;

        my @sources = @$sources_reference,

	my (%datasets) = %$datasets_reference;

        my %geneModel_stats = %$gene_model_total_stats_reference;


	#Set arrays for exons classes
        my @statsExons = (
                     "matched_exons",
                     "partial_match",
                     "contradictory_exons",
                     "unsupported_exons",
                     "micro_exons",
                     "not_aligned_exons"
                    );

	#Set arrays for genes classes
        my @statsModel = (
                     "match",
                     "partialMatch",
                     "fuzzy",
                     "contradictory",
                    );

#Keys checkpoint
#        foreach my $key (keys %geneModel_stats ) {
#        	foreach my $key2 (sort (keys %{$geneModel_stats{$key}} ) ) {
#	        	foreach my $key3 (sort (keys %{$geneModel_stats{$key}{$key2}} ) ) {
#		        	print $key, " ", $key2, " ",  $key3, " ", $geneModel_stats{$key}{$key2}{$key3}, "\n";
#                        } #end foreach
#                } #end foreach
#        } #end foreach


	#Print to the exons file
        open(EXONS, ">>", $outputFile) || die "can't open file $outputFile\n";
	chomp($geneModel);
        print EXONS $fileName, "\t", $geneModel;

	#loop throught the exons
	for(my $i = 0; $i <= $exons_gene; ++$i) {

        	print EXONS "\t";

		#loop through the categories
        	foreach my $category (@sources) {

			# if is not a repeat
                	if( $datasets{$category}[1] ne "Repeats" && $datasets{$category}[1] ne "Transposons" ) {

                        	foreach my $stat (@statsExons) {

	                            	if( defined($geneModel_stats{$category}{"exons"}{$i}{$stat}) ) {
			                        print EXONS $geneModel_stats{$category}{"exons"}{$i}{$stat}, "/";
					} else {
                                		print EXONS "0/";
                        	        } #end else match exons

                                } #end foreach stats loop

			} else {

                            	if( defined($geneModel_stats{$category}{"exons"}{$i}) ) {
	                        	print EXONS $geneModel_stats{$category}{"exons"}{$i};
				} else {
                                	print EXONS "0";
                                } #end else match exons

                        } #end else

	        	print EXONS ":";

                } #end foreach

	} #end for

       	print EXONS "\n";

	close (EXONS);


	#Print to the gene models file
        open(MODEL, ">>$output/Evidence_Counters_perGeneModel.txt") || die "can't open file $output/Evidence_Counters_perGeneModel.txt\n";
	chomp($geneModel);
        print MODEL $fileName, "\t", $geneModel, "\t", (1+$exons_gene), "\t";

	#loop through the categories
       	foreach my $category (@sources) {

		#if it is not a repeat
               	if( $datasets{$category}[1] ne "Repeats" && $datasets{$category}[1] ne "Transposons" ) {

			foreach my $stat (@statsModel) {

		                if( defined($geneModel_stats{$category}{"total"}{$stat}) ) {
					print MODEL $geneModel_stats{$category}{"total"}{$stat}, "/";
				} else {
		                      	print MODEL "0/";
				} #end else match exons

			} #end foreach stats loop

		} else {

	                if( defined($geneModel_stats{$category}{"total"}) ) {
				print MODEL $geneModel_stats{$category}{"total"};
			} else {
	                      	print MODEL "0";
			} #end else match exons


		} #end else for category

		print MODEL "\t";

	} #end foreach

	print MODEL "\n";

} #end sub


####################################
#                                  #
#       Usage printing sub         #
#                                  #
####################################

sub usage {

	print "
	usage: .pl -l source-types.txt -x 5 -o output -f/-D
	options:
	  -1 <file>     Source/Types file
                        Tab Delimited file with three colums, two first colums should be taken from the gff files with this command line
                        'cut -2,3 <gffFIle> | sort | uniq > source-types.txt'
                        When selecting the rows, make sure to extract the ones that contain the exons, not the main feature (Usually denoted by having a 'Parent' attribute.
                        Third column should be manually add with the categories for the evidence, there are three predefine categories that
                        have a differente behavior, need to be type has it follow
                        	GeneModels: Automated gene models that will be analyze
				Predictors: Annotation generated from the Ab Initio prediction programs.
                        	Transcriptome: If among the evidence are EST, RNAseq data or transcriptome assemblies they should be input has transcriptome
                                       On the validation section the last exons are allow to be longer than the CDS.
                        	Transposons: Not much description need it, are transposons.
                        	Repeats: This category is specially define for tandem repeats, the weighting is different than for transposons.

                        Exampl of the Sources/Types file

                        AUGUSTUS	match_part	Predictors
                        Cuff	mRNA	GeneModels
                        est2genome	match_part	Transcriptome
                        ESTs	match_part	Transcriptome
                        GLEAN	CDS	GeneModels
                        GlimmerHMM	match_part	Predictors
                        maker	CDS	GeneModels
                        RepeatMasker	match	Transposons
                        RepeatProteinMask	match	Transposons
                        TRF	match	Repeats

	  -x <int>      Adjustment value for the exon/intron boundaries. Value depend of category of evidence as it follow
				Predictors: no adjustment is permitted, boundaries of the gene model need to do a perfect match with the prediction.
				Transcriptome: no modification to the inputted adjustment value.
				Everything else: two times the adjustment is allow between the boundaries of the gene model and the boundaries of the evidence.

	  -o <string>   Output directory to be created
                        A set of files will be created on this directory

                        blablabla.txt:
                        blablabla.txt:
                        blablabla.txt:


          -f <file>     GFF file to be analyze (only one)

	  -D <string>   Directory with the files to analyze (bacth mode) (/home/user/directory/)
                        In case of batch mode, the file extension should be gff, only this ones will be considere for analysys.
                        If yoy have other file extension, this can be change in the line 98 of this script.

          -a            Flag for annotated gene models (has is writen in gff attributes line of the mRNA features).

	  -r		Activate recursive mode (if -D is provided, it will look within folders in the input folder for .gff files

	";

} #end of sub usage
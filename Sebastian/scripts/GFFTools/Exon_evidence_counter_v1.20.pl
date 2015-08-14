#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_l $opt_x $opt_o $opt_a $opt_D $opt_f);
getopts('hl:x:o:a:D:f:');

########################################################
########################################################

#                                                      #
#                                                      #
#           Exon Evidence Counter v1.20                #
#                   05-23-2012                         #
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

##################################
### Difference for version 1.0 ###
##################################
#
# Modified adjustment value
# Add section for handle annotated files
#


#### User ask for extended help documentation ####

if (defined($opt_h)) {

	#option h was provide, print full description

        usage();

        die "";

} elsif ( !(defined($opt_l)) || !(defined($opt_x)) || ( !(defined($opt_D)) && !(defined($opt_f))) || !(defined($opt_o)) ) {

	# At least one of the requried options is missing

	die "
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

        For a more detailed description and quick tips for tuning the script use option -h.
	";
} #end if options


#### Check arguments and prepare log file ####

#Get inputed arguments
my $sourcesFile = $opt_l;
my $adj = $opt_x;
my $output = $opt_o;


#Create output directory unless it already exists
system("mkdir $output") unless -d "$output";

#Print require arguments to screen
print "Analyze options:\n";
print "Sources file = $sourcesFile\n";
print "Adjustment value = $adj\n";
print "Output directory = $output\n\n";

#Open and print arguments to log file
open (LOG, ">$output/log");
print LOG "Analyze options:\n";
print LOG "Sources file = $sourcesFile\n";
print LOG "Adjustment value = $adj\n";
print LOG "Output directory = $output\n";

#Print begining of command line to log file
print LOG "Command line : Exon_evidence_counter_v1.0.pl -l $sourcesFile -x $adj -o $output ";

#Define if it will be working for one file or a directory
my @files;

#will analyze a directory of files
if (defined($opt_D)) {

        print "Analyzing directory = $opt_D\n";

        print LOG "-D $opt_D\n\n";

	#Open working directory and get files to analyze
        open(LOGFILES, ">$output/Analyzed_Files.log");
        opendir(DIR, $opt_D) || die "Can't open directory $opt_D\n";
        while (my $input = readdir(DIR)) {

		#Verify extension of the file
		my $extension = substr($input,-3,);
                if( $extension eq "gff") {

                	push(@files, "$opt_D/$input");
                        print LOGFILES $input, "\n";

                } #end if for gff files

   	} #end foreach directory

        if( !(defined($files[0])) ) {

        	print "Directory provided don't contain any .gff files.\n";
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

print "Reading sources and types from inputed file.\n";

# open source related variables
my @tracks;
my $numTracks;
my $inputValues;
my %datasets;

open(SOURCES, $sourcesFile);

print LOG "Inputed evidence sources-types and categories\n";
print LOG "Sources\tTypes\tCategories\n";

#Load sources information into hash
while(my $source = <SOURCES>) {

	chomp($source);

	my @valuePairs = split("\t", $source);

	push(@tracks, [$valuePairs[0], $valuePairs[1], $valuePairs[2]]);

	++$inputValues;

	$datasets{$valuePairs[2]} = 1;

	print LOG "$valuePairs[0]\t$valuePairs[1]\t$valuePairs[2]\n";

} #end while

print LOG "\n\n";

my @categories = keys(%datasets);

#check existing categories
if(!(exists($datasets{'GeneModels'}))) {

	print "Please define at least one track on the 'GeneModels' category on your input file.\n";
	die;

} elsif ( !(exists($datasets{'Transcriptome'})) ) {

	print "Transcriptome category wasn't found, analysis will continue without it.\n";
	print "If you do have transcriptome-like date, please check the sources-types file.\n\n";

}

if ( !(exists($datasets{'Transposons'})) ) {

	print "Transposons category wasn't found, analysis will continue without it.\n";
	print "If you do have transposon data, please check the sources-types file.\n\n";

}

if ( !(exists($datasets{'Repeats'})) ) {

	print "Repeats category wasn't found, analysis will continue without it.\n";
	print "If you do have repeats data, please check the sources-types file.\n\n";

} #end if for categories

print "Sources and Types loaded from $sourcesFile, $inputValues lines of evidence will be use for counting.\n\n";



open(ERRORS, ">$output/error.log");


#### Prepare output files ####

#If annotation flag provided open annotated results files

my $annotationFlag = "";
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
		print ANNOTATED "#First level\n";
		print ANNOTATED "#Exon1\tExon2\tExon3\t..\tExon(n)\n";
		print ANNOTATED "#Second level\n";

		foreach my $category (@categories) {

			print ANNOTATED $category, ":";

		} #end foreach

		print ANNOTATED "\n";

		print ANNOTATED "#Third level\n";
		print ANNOTATED "#matched_exons/partial_match/contradictory_exons/unsupported_exons/micro_exons/not_aligned_exons\n";

		close (ANNOTATED);

	} #end foreach annotated files

} #end if for annotation flag

#Open Counters files

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
	print EXON "#First level\n";
	print EXON "#Exon1\tExon2\tExon3\t..\tExon(n)\n";
	print EXON "#Second level\n";

	foreach my $category (@categories) {
		print EXON $category, ":";
	} #end foreach

	print EXON "\n";
	print EXON "#Third level\n";
	print EXON "#matched_exons/partial_match/contradictory_exons/unsupported_exons/micro_exons/not_aligned_exons\n";

	close(EXON);

} #end foreach

#Create directory for evidence status files
system("mkdir $output/Evidence_status") unless -d "$output/Evidence_status";



####################################
#                                  #
#     Main analysis pipeline       #
#                                  #
####################################

print "Starting the analysis.\n\n";

my $fileProgress = 1;
my $analizedfile = 0;

foreach my $file (@files) {

	if( (int($analizedfile/$numFiles) ) > $fileProgress) {

		print int($analizedfile/$numFiles), " ";

		$fileProgress = int($analizedfile/$numFiles);

	} #end if

	my %evidence;

	my @filePath = split("/", $file);

	my $fileName = pop(@filePath);

	open(INPUT, $file);
        open(STATS, ">$output/Evidence_status/$fileName.evidence_stats.txt");

	my $annotated = "FALSE";

	while (my $line = <INPUT>) {

		chomp($line);

		my @fields = split("\t", $line);

		my $score = 0;

		if(scalar(@fields) > 1) {

				if($fields[2] eq "mRNA") {

					$annotated = "FALSE";

					my @attributes = split(";", $fields[8]);

					foreach my $attribute (@attributes) {

						if ($attribute eq $annotationFlag) {

							$annotated = "TRUE";

							last if $annotated == "TRUE";

						} else {

							$annotated = "FALSE";

						} #end if annotated

					} #end foreach

				} #end if for mRNA features

			for(my $i = 0; $i < scalar(@tracks); ++$i) {

				if($fields[1] eq $tracks[$i][0]) {

					if($fields[5] ne ".") {

						$score = $fields[5];

					} #end if score present

				} #end if

				if($fields[1] eq $tracks[$i][0] && $fields[2] eq $tracks[$i][1]) {

					if($fields[6] ne ".") {

						$score = $fields[5];

					} #end if score present

					my @attributes = split(";", $fields[8]);

					my $parent;

					foreach my $attribute(@attributes) {

						if ($attribute =~ /^Parent=(\S+)/) {

							my @data = split ("=", $attribute);

							$parent = $data[1];

						} elsif ( ($attribute =~ /^Name=(\S+)/) && !(defined($parent)) ) {

							my @data = split ("=", $attribute);

							$parent = $data[1];

						} #end elsif parent

					} #end foreach

                                        if ($fields[6] eq "+") {

						push(@{$evidence{$tracks[$i][2]}{$parent}}, [$fields[3], $fields[4], $score, $fields[6], $annotated]);

                                        } else {

						unshift(@{$evidence{$tracks[$i][2]}{$parent}}, [$fields[3], $fields[4], $score, $fields[6], $annotated]);

                                        } #end else

				} #end if

			} #end for tracks

		} #end if for gff line

	} #end while

	#print "Analyzing file $fileName.\n";

	foreach my $geneModel (keys (%{$evidence{'GeneModels'}})) {

		my $exons_gene = scalar(@{$evidence{'GeneModels'}{$geneModel}});

		--$exons_gene;

		my @gene_exon_structure = @{$evidence{'GeneModels'}{$geneModel}};

                my %gene_model_total_stats;

                #print "\nAnalyzing $geneModel\n";

		my $done_flag = "FALSE";

		for(my $i = 0; $i <= 4; ++$i) {

			if ($exons_gene eq $i || ( ($i eq 4) && ($exons_gene >= $i) ) ) {

				foreach my $category (@categories) {

					my $adjustment;

					if ($category eq "Predictors") {

						$adjustment = 0;

					} elsif ($category eq "Transcriptome") {

						$adjustment = $adj;

					} else {

						$adjustment = $adj+$adj;

					} #end if adjustment settings

        	                	if( $category ne "Repeats" && $category ne "Transposons") {

						$gene_model_total_stats{$category} = { check_evidence($fileName, $category, $exons_gene, $geneModel, $adjustment, \@gene_exon_structure, \%evidence) };

        	                        } else {

						$gene_model_total_stats{$category} = { check_repeats($category, $exons_gene, $geneModel, $adjustment, \@gene_exon_structure, \%evidence) };

        	                        } #end else

        	                	#print "Category $category done \n";

				} #end foreach

        	                my $outputFile = $outputEvidenceFiles[$i];

        	                #print "Summarizing evidence\n";

        	                summarize_evidence($fileName, $outputFile, $geneModel, $exons_gene, \@categories, \%gene_model_total_stats);

#				print $gene_exon_structure[0][4], "\n";

				if($gene_exon_structure[0][4] eq "TRUE") {

					#print "Output Send to annotated gene models\n";

	        	                my $outputFile = $outputAnnotatedFiles[$i];

	        	                summarize_evidence($fileName, $outputFile, $geneModel, $exons_gene, \@categories, \%gene_model_total_stats);

				} #end if for annotated outputing


			} # if check number exons

		} #end for exons number

	} #end foreach

	++$analizedfile;

} #end for (files)


exit;

###########################################################################################
###########################################################################################


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

		        	#print "Transcript $transcript within the gene model, been analyzed\n";

	                        #print $geneModel, "\t", $transcript, "\n";

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
	                                if( ($geneExons[$exonGene][0]-$adjustment) > $eviExons[$exonEvi][0]) { $start_start = "before"

	                                } elsif ( ($geneExons[$exonGene][0]+$adjustment) < $eviExons[$exonEvi][0]) {    $start_start = "after"

	                                } else { $start_start = "match"

	                                } #end else for start-start comparison


	                                #Start of the gene against end of the evidence
	                                if( ($geneExons[$exonGene][0]-$adjustment) > $eviExons[$exonEvi][1]) { $start_end = "before"

	                                } elsif ( ($geneExons[$exonGene][0]+$adjustment) < $eviExons[$exonEvi][1]) {    $start_end = "after"

	                                } else { $start_end = "match"

	                                } #end else for start-end comparison


	                                #End of the gene against start of the evidence
	                                if( ($geneExons[$exonGene][1]-$adjustment) > $eviExons[$exonEvi][0]) { $end_start = "before"

	                                } elsif ( ($geneExons[$exonGene][1]+$adjustment) < $eviExons[$exonEvi][0]) {    $end_start = "after"

	                                } else { $end_start = "match"

	                                } #end else for end-start comparison


	                                #End of the gene against end of the evidence
	                                if( ($geneExons[$exonGene][1]-$adjustment) > $eviExons[$exonEvi][1]) { $end_end = "before"

	                                } elsif ( ($geneExons[$exonGene][1]+$adjustment) < $eviExons[$exonEvi][1]) {    $end_end = "after"

	                                } else { $end_end = "match"

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

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #Exon evidence before exon of gene
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "after") && ($end_end eq "after") ) {

	                                        #print "1 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

	                                        if($within eq "TRUE") {

	                                                $gene_exon_stats{$exonGene} = "Unsupported_exon";

	                                        } #end if

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #Exon gene before exon of evidence
	                                } elsif ( ($start_start eq "before") && ($start_end eq "before") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "2 $exonEvi $exonGene\n";

	                                        if($within eq "FALSE") {

	                                                $evi_exon_stats{$exonEvi} = "Extra_exon_ends";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Extra_exon_within";

	                                        } #end else

	                                        ++$exonEvi;

	                                #Overlaping exons (evidence first)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "3 $exonEvi $exonGene\n";

	                                        if($category eq "Transcriptome" && $exons_transcript eq 0 && $exonGene eq 0) {

	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

	                                        ++$exonEvi;

	                                #Overlaping exons (gene first)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

	                                        #print "4 $exonEvi $exonGene\n";

	                                        if($category eq "Transcriptome" && $exons_transcript eq 0 && $exonGene eq $exons_gene) {

	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #Partial match (end matches) (evidence shorter)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

	                                        #print "5 $exonEvi $exonGene\n";

	                                        if($exonEvi eq 0) {

	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Partial match (end matches) (gene shorter)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

	                                        #print "6 $exonEvi $exonGene\n";

	                                        if($category eq "Transcriptome" && $exonGene eq 0 && $exonEvi eq 0) {

	                                                $evi_exon_stats{$exonEvi} = "Match";
	                                                $gene_exon_stats{$exonGene} = "Match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Partial match (start matches) (evidence shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "7 $exonEvi $exonGene\n";

	                                        if($exonEvi eq $exons_transcript) {

	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

	                                        ++$exonEvi;

	                                #Partial match (start matches) (gene shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

	                                        #print "8 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

	                                        if($category eq "Transcriptome" && $exonGene eq $exons_gene && $exonEvi eq $exons_transcript) {

	                                                $evi_exon_stats{$exonEvi} = "Match";
	                                                $gene_exon_stats{$exonGene} = "Match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #else

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #Evidence exon within gene exon
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "9 $exonEvi $exonGene\n";

	                                        if($exons_transcript eq 0) {

	                                                $evi_exon_stats{$exonEvi} = "Partial_match";
	                                                $gene_exon_stats{$exonGene} = "Partial_match";

	                                        } else {

	                                                $evi_exon_stats{$exonEvi} = "Contradictory";
	                                                $gene_exon_stats{$exonGene} = "Contradictory";

	                                        } #end else

	                                        ++$exonEvi;

	                                #Gene exon within evidence exon
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

	                                        #print "10 $exonEvi $eviExons[$exonEvi][0] $eviExons[$exonEvi][1] $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1]\n";

	                                        $evi_exon_stats{$exonEvi} = "Contradictory";
	                                        $gene_exon_stats{$exonGene} = "Contradictory";

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #Overlapping end-start
	                                } elsif ( ($start_start eq "before") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "11 $exonEvi $exonGene\n";

	                                        ++$exonEvi;

	                                #Overlapping start-end
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "after") ) {

	                                        #print "12 $exonEvi $exonGene\n";

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;

	                                #
	                                } elsif ( ($start_start eq "match") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

	                                        #print "13 $exonEvi $exonGene\n";

	                                        ++$exonEvi;

	                                #
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "match") ) {

	                                        #print "14 $exonEvi $exonGene\n";

						if($exonGene eq $exons_gene) {

							++$exonEvi;

						} #end if

	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                #Perfect match
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

	                                        #print "Match $exonEvi $exonGene\n";

	                                        $evi_exon_stats{$exonEvi} = "Match";
	                                        $gene_exon_stats{$exonGene} = "Match";

	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                } else {

#	                                        print "Whatever $exonEvi $exonGene $geneExons[$exonGene][0] $geneExons[$exonGene][1] $eviExons[$exonEvi][0] $eviExons[$exonEvi][1]\n";
#	                                        print "$start_start $start_end $end_start $end_end\n";

						chomp($geneModel);

	                                        print ERRORS "In $fileName none of the comparisons return TRUE for $transcript, exon $exonEvi, against $geneModel, exon $exonGene.\n";
						print ERRORS "        $start_start $start_end $end_start $end_end\n";
	                                        ++$exonGene;
	                                        ++$exonEvi;

	                                } # end else check status

	                        } #end while evidence exons

	                        my $match_exons = 0;
	                        my $partial_match = 0;
	                        my $contradictory_exons = 0;
	                        my $unsupported_exons = 0;
	                        my $micro_exons = 0;

	                        print STATS $geneModel, "\t", $transcript;

	                        for(my $i = 0; $i <= $exons_gene; ++$i) {
                                        if( defined($gene_exon_stats{$i}) ) {

	                                        print STATS "\t", $gene_exon_stats{$i};

                                                if($gene_exon_stats{$i} eq "Match") { ++$geneModel_stats{$i}{'matched_exons'} }

                                                elsif($gene_exon_stats{$i} eq "Partial_match") { ++$geneModel_stats{$i}{'partial_match'} }

                                                elsif($gene_exon_stats{$i} eq "Contradictory") { ++$geneModel_stats{$i}{'contradictory_exons'} }

                                                elsif($gene_exon_stats{$i} eq "Unsupported_exon") { ++$geneModel_stats{$i}{'unsupported_exons'} }

                                                elsif($gene_exon_stats{$i} eq "Micro_exon") { ++$geneModel_stats{$i}{'micro_exons'}}


	                                } else {

	                                        print STATS "\t";
                                                ++$geneModel_stats{$i}{'not_aligned_exons'}

	                                } #end else for defined exon stat

#	                                print $category, " ", $geneModel_stats{$i}{'match_exons'}, "\n";

	                        } #end foreach

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

	foreach my $repeat ( keys %{$$evidence{$category}}) {

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

                                        } else {

	                                        print STATS "\t";

	                                } #end else for defined exon stat

                                        ++$exonGene;

                                } #end while loop for exons

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

	my ($fileName,$outputFile, $geneModel, $exons_gene, $categories_reference, $gene_model_total_stats_reference) = @_;

	my (@categories) = @$categories_reference;

        my %geneModel_stats = %$gene_model_total_stats_reference;

        my @stats = (
                     "matched_exons",
                     "partial_match",
                     "contradictory_exons",
                     "unsupported_exons",
                     "micro_exons",
                     "not_aligned_exons"
                    );

#        foreach my $key (keys %geneModel_stats ) {
#        	foreach my $key2 (sort (keys %{$geneModel_stats{$key}} ) ) {
#	        	foreach my $key3 (sort (keys %{$geneModel_stats{$key}{$key2}} ) ) {
#		        	print $key, " ", $key2, " ",  $key3, " ", $geneModel_stats{$key}{$key2}{$key3}, "\n";
#                        } #end foreach
#                } #end foreach
#        } #end foreach

        open(OUTPUT, ">>", $outputFile) || "can't open file $outputFile\n";

	chomp($geneModel);

        print OUTPUT $fileName, "\t", $geneModel;

	for(my $i = 0; $i <= $exons_gene; ++$i) {

        	print OUTPUT "\t";

        	foreach my $category (@categories) {

                	if( $category ne "Repeats" && $category ne "Transposons" ) {

                        	foreach my $stat (@stats) {

	                            	if( defined($geneModel_stats{$category}{$i}{$stat}) ) {
			                        print OUTPUT $geneModel_stats{$category}{$i}{$stat}, "/";
					} else {
                                		print OUTPUT "0/";
                        	        } #end else match exons

                                } #end foreach stats loop

			} else {

#                        	print "Problems with repeats\n";

                            	if( defined($geneModel_stats{$category}{$i}) ) {
	                        	print OUTPUT $geneModel_stats{$category}{$i};
				} else {
                                	print OUTPUT "0";
                                } #end else match exons

                        } #end else

	        	print OUTPUT ":";

                } #end foreach

	} #end for

       	print OUTPUT "\n";

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

	  -x <int>      Adjustment value for the exon/intron boundaries

	  -o <string>   Output directory to be created
                        A set of files will be created on this directory

                        blablabla.txt:
                        blablabla.txt:
                        blablabla.txt:


          -f <file>     GFF file to be analyze (only one)

	  -D <string>   Directory with the files to analyze (bacth mode) (/home/user/directory/)
                        In case of batch mode, the file extension should be gff, only this ones will be considere for analysys.
                        If yoy have other file extension, this can be change in the line 98 of this script.
	";

} #end of sub usage
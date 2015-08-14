#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_o);
getopts('h:i:o:');

########################################################
########################################################
#                                                      #
#                                                      #
#             GFF geneStats calculator                 #
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


#### User ask for extended help documentation ####

if (defined($opt_h)) {

	#option h was provide, print full description

        usage();

        exit;

} elsif ( !(defined($opt_i)) || !(defined($opt_o)) ) {

	# At least one of the requried options is missing

	print "
	usage: .pl -i GFF_geneFile
	options:


          Required Options
          -i <file>     GFF file to be analyze (only one)
	  -o <prefof>	Output prefix to be used

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
my $input = $opt_i;
my $prefix = $opt_o;



#Open output files
open (LOG, ">$prefix.log");
open (GENEMODELSTATS, ">$prefix.genemodelstats.txt");
open (EXONSIZES, ">$prefix.ExonSizes.txt");
open (INTRONSIZES, ">$prefix.IntronSizes.txt");
open (GENERALSTATS, ">$prefix.generalStats.txt");

print GENEMODELSTATS join("\t", ("gene_lodel", "num_exons", "transcript_length", "intron_length", "CDS_length", "mean_exon_length", "mean_intron_length") ), "\n";

print GENERALSTATS join("\t", ("Num_gene_models", "Num_monoexonic_gene_models", "mean_transcript_length", "mean_CDS_length", "mean_num_exons", "mean_exon_length", "max_exon_length", "mean_intron_length", "max_intron_length") ), "\n";

#Print require arguments to screen
print "\nAnalysis options:\n";
print "Input GFF file = $input\n";
print "Output prefix = $prefix\n";

#Open and print arguments to log file
print LOG "Analysis options:\n";
print LOG "Input GFF file = $input\n";
print LOG "Output prefix = $prefix\n";

#Print begining of command line to log file
print LOG "Command line : Exon_evidence_counter_v2.1.pl -i $input -o $prefix ";



####################################
#                                  #
#     Main analysis pipeline       #
#                                  #
####################################

print STDERR "Starting the analysis.\n\n";

my %evidence;
my @filePath = split("/", $input);
my $fileName = pop(@filePath);

#open input file
open(INPUT, $input) or die "File $input not found\n";

while (my $line = <INPUT>) {

	chomp($line);

	my @fields = split("\t", $line);

	#if line contains proper number of fields
	if(defined($fields[8])) {

		#check if source/type set is need it
		if( ($fields[2] eq "CDS") || ($fields[2] eq "exon") ) {

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

				push(@{$evidence{$parent}{$fields[2]}}, [$fields[3], $fields[4], $fields[6]]);

                        } else {

				unshift(@{$evidence{$parent}{$fields[2]}}, [$fields[3], $fields[4], $fields[6]]);

                        } #end else

		} #end if

	} #end if for gff line exists

} #end reading file

print STDERR "Finish reading file, while generate exon statistics\n\n";

# Create data storage variables
my @exonLengths;
my @intronLengths;
my @CDSLengths;
my @transcriptLengths;

my $monoExonicGeneModels;

my $exonTotalSum;
my $intronTotalSum;
my $CDSTotalSum;

my $numExonTotalSum;

#Loop through all the gene models in the hash
foreach my $geneModel (keys (%evidence)) {

	my $exonsSum = 0;
	my $intronsSum = 0;
	my $CDSsSum = 0;

	### Calculating values for the exons and introns

	#load exons from main hash
	if(!defined($evidence{$geneModel}{'exon'})) {
		print STDERR "Error: 'exon' feature missing for $geneModel\n";
		delete ($evidence{$geneModel});
		next;
	} #end if exon exists
	my @gene_exon_structure = @{$evidence{$geneModel}{'exon'}};

	#calculate the number of exons
	my $exons_gene = scalar(@gene_exon_structure);

	$numExonTotalSum += $exons_gene;

	if($exons_gene == 1) {++$monoExonicGeneModels}

	my $firstExonSize = $gene_exon_structure[0][1] - $gene_exon_structure[0][0];
	push(@exonLengths, $firstExonSize);
	$exonsSum += $firstExonSize;

	my $prevEnd = $gene_exon_structure[0][1];

	for( my $exon = 1; $exon < $exons_gene; ++$exon) {

		my $exonSize = $gene_exon_structure[$exon][1] - $gene_exon_structure[$exon][0];
		push(@exonLengths, $exonSize);
		$exonsSum += $exonSize;

		my $intronSize = abs($gene_exon_structure[$exon][0] - $prevEnd);
		push(@intronLengths, $intronSize);
		$intronsSum += $intronSize;

		$prevEnd = $gene_exon_structure[$exon][1];
			
	} #end for each exon

	my @gene_CDS_structure = @{$evidence{$geneModel}{'CDS'}};
	my $CDSs_gene = scalar(@gene_CDS_structure);

	for( my $CDS = 0; $CDS < $CDSs_gene; ++$CDS) {

		my $CDSSize = $gene_CDS_structure[$CDS][1] - $gene_CDS_structure[$CDS][0];
		push(@CDSLengths, $CDSSize);
		$CDSsSum += $CDSSize;

	} #end for each exon

	my $meanExonLength = $exonsSum/$exons_gene;

	my $meanIntronLength = 0;
	if( ($exons_gene-1) > 0) {$meanIntronLength = $intronsSum/$exons_gene-1}

	print GENEMODELSTATS join( "\t", ($geneModel, $exons_gene, $exonsSum, $intronsSum, $CDSsSum, sprintf("%.2f", $meanExonLength), sprintf("%.2f", $meanIntronLength)) ), "\n";

	$exonTotalSum += $exonsSum;
	$intronTotalSum += $intronsSum;
	$CDSTotalSum += $CDSsSum;

} #end foreach genemodel

my $numGeneModels = scalar(keys %evidence);

foreach my $exonSize (@exonLengths) { print EXONSIZES $exonSize, "\n" }
foreach my $intronSize (@intronLengths) { print INTRONSIZES $intronSize, "\n" }

my @sortedExonLengths = sort { $a <=> $b} @exonLengths;
my @sortedIntronLengths = sort { $a <=> $b} @intronLengths;

my $maxExon = pop(@sortedExonLengths);
my $maxIntron = pop(@sortedIntronLengths);

my $meanCDSlength = $CDSTotalSum/$numGeneModels;
my $meanTranscriptlength = $exonTotalSum/$numGeneModels;

my $meanNumberExons = $numExonTotalSum/$numGeneModels;

my $meanExonSize = $exonTotalSum/$numExonTotalSum;
my $meanIntronSize = $intronTotalSum/$numExonTotalSum;

print GENERALSTATS join("\t", ($numGeneModels, $monoExonicGeneModels, sprintf("%.2f", $meanTranscriptlength), sprintf("%.2f", $meanCDSlength), sprintf("%.2f", $meanNumberExons), sprintf("%.2f", $meanExonSize), $maxExon, sprintf("%.2f", $meanIntronSize), $maxIntron) ), "\n";

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
#       Usage printing sub         #
#                                  #
####################################

sub usage {

	print "
	usage: .pl -l source-types.txt -x 5 -o output -f/-D
	";

} #end of sub usage

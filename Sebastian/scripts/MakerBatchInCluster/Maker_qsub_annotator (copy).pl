#!/usr/bin/perl

######################################################
#                                                    #
#               Maker qsub Annotator                 #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to run the maker annotation pipeline on a inputed fasta file
# splitting the file in a set of subfasta and running then on parallel in a cluster system under qsub

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use Getopt::Long;
use POSIX;

#########################################
#                                       #
#          Reading in variables         #
#                                       #
#########################################

### Working variables
my $pwd = getcwd;
my $time = getTime();

my $inputFile;
my $gffprefix;
my $threads = 1;
my $numseq;
my $help;
my $retry = "no";
my $qsubsFile;
my $queue;

### Maker behavior variables
my $again = " ";
my $force = " ";


### Maker opts variables

# EST evidence
my $est;
my $est_reads;
my $altest;
my $est_gff;
my $altest_gff;

# Protein evidence
my $protein;
my $protein_gff;

# Repeat masking
my $model_org;
my $rmlib;
my $repeat_protein;
my $rm_gff;

# Ab initio prediction
my $snaphmm;
my $gmhmm;
my $augustus_species;
my $pred_gff;
my $model_gff;
my $est2genome;

# Other information
my $other_gff;

# Behavior options
my $max_dna_len;
my $min_contig;

my $pred_flank;
my $AED_threshold;
my $min_protein;
my $alt_splice;
my $always_complete;
my $map_forward;
my $keep_preds;

my $split_hit;
my $softmask;
my $single_exon;
my $single_length;

my $tmp = "maker_tmp";

GetOptions (
		'help' => \$help,
		'input=s' => \$inputFile,
		'gffprefix=s' => \$gffprefix,
                'threads=i' => \$threads,
                'numseq=i' => \$numseq,
                'retry=s' => \$retry,
		'qsubsFile=s' => \$qsubsFile,
		'queue=s' => \$queue,

		'again' => \$again,
		'force' => \$force,

	 	'est=s' => \$est,
	 	'est_reads=s' => \$est_reads,
	 	'altest=s' => \$altest,
	 	'est_gff=s' => \$est_gff,
	 	'altest_gff=s' => \$altest_gff,

                'protein=s' => \$protein,
		'protein_gff=s' => \$protein_gff,

                'model_org=s' => \$model_org,
                'rmlib=s' => \$rmlib,
                'repeat_protein=s' => \$repeat_protein,
		'rm_gff=s' => \$rm_gff,

                'snaphmm=s' => \$snaphmm,
                'gmhmm=s' => \$gmhmm,
                'augustus_species=s' => \$augustus_species,
                'pred_gff=s' => \$pred_gff,
                'model_gff=s' => \$model_gff,
                'est2genome=i' => \$est2genome,

                'other_gff=s' => \$other_gff,

                'max_dna_len=i' => \$max_dna_len,
                'min_contig=i' => \$min_contig,

                'pred_flank=i' => \$pred_flank,
                'AED_threshold=i' => \$AED_threshold,
                'min_protein=i' => \$min_protein,
                'alt_splice=i' => \$alt_splice,
                'map_forward=i' => \$map_forward,
                'keep_preds=i' => \$keep_preds,

                'split_hit=i' => \$split_hit,
                'softmask=i' => \$softmask,
                'single_exon=i' => \$single_exon,
                'single_length=i' => \$single_length,

		'TMP=s' => \$tmp

           );

if (defined($help)) {

	#option help was provide, print full description

        usage();

        exit;

} elsif ( (!defined($inputFile) && !defined($qsubsFile)) || !(defined($gffprefix)) || ( ($inputFile =~ /.*\\.((f|F)(a|A))|((f|F)(a|A)(s|S)(t|T)(a|A))$/) && !defined($numseq) ) ) {

	print "Missing require arguments\n";
	print "Usage:\n";
	print "maker_automatic_annotator_v1.pl -input <FASTA/TXT> -numseq <INT> -gffprefix <prefix for outputgff> other_options\n";
	print "Use -help for detail help\n";
        die "Please try again.\n"

} #end if usage

print STDERR "\nInputted variables\n\n";

if(defined($inputFile)) {print STDERR "input		= $inputFile\n"}
if(defined($qsubsFile)) {print STDERR "input		= $qsubsFile\n"}

print STDERR "prefix for gff's	= $gffprefix\n";

if(defined($numseq)) {print STDERR "Num Sequences	= $numseq\n"}

print STDERR "Num threads	= $threads\n";

print STDERR "Retry failed sequences	= $retry\n\n";


print STDERR "Maker Behavior Variables\n\n";

if(defined($again)) {print STDERR "Annotation files will be regenerate it\n"; $again = '-again'} else {$again = ''}
if(defined($force)) {print STDERR "All results will be deleted it before reRun (Blas analysis will need to be redone)\n"; $force = '-force'} else {$force = ''}

#### Hashing variables ####
print STDERR "\nMaker Variables for opts file\n\n";

my %maker_variables;

if (defined($est)) {$maker_variables{'est'} = check_path($est); print STDERR "est		= $est\n"}
if (defined($est_reads)) {$maker_variables{'est_reads'} = check_path($est_reads); print STDERR "est_reads		= $est_reads\n"}
if (defined($altest)) {$maker_variables{'altest'} = check_path($altest); print STDERR "altest		= $altest\n"}
if (defined($est_gff)) {$maker_variables{'est_gff'} = check_path($est_gff); print STDERR "est_gff		= $est_gff\n"}
if (defined($protein)) {$maker_variables{'protein'} = check_path($protein); print STDERR "protein		= $protein\n"}
if (defined($protein_gff)) {$maker_variables{'protein_gff'} = check_path($protein_gff); print STDERR "protein_gff		= $protein_gff\n"}
if (defined($model_org)) {$maker_variables{'model_org'} = $model_org; print STDERR "model_org		= $model_org\n"}
if (defined($rmlib)) {$maker_variables{'rmlib'} = check_path($rmlib); print STDERR "rmlib		= $rmlib\n"}
if (defined($repeat_protein)) {$maker_variables{'repeat_protein'} = check_path($repeat_protein); print STDERR "repeat_protein		= $repeat_protein\n"}
if (defined($rm_gff)) {$maker_variables{'rm_gff'} = check_path($rm_gff); print STDERR "rm_gff		= $rm_gff\n"}
if (defined($snaphmm)) {$maker_variables{'snaphmm'} = check_path($snaphmm); print STDERR "snaphmm		= $snaphmm\n"}
if (defined($gmhmm)) {$maker_variables{'gmhmm'} = check_path($gmhmm); print STDERR "gmhmm		= $gmhmm\n"}
if (defined($augustus_species)) {$maker_variables{'augustus_species'} = $augustus_species; print STDERR "augustus_species = $augustus_species\n"}
if (defined($pred_gff)) {$maker_variables{'pred_gff'} = check_path($pred_gff); print STDERR "pred_gff		= $pred_gff\n"}
if (defined($model_gff)) {$maker_variables{'model_gff'} = check_path($model_gff); print STDERR "model_gff		= $model_gff\n"}
if (defined($est2genome)) {$maker_variables{'est2genome'} = $est2genome; print STDERR "est2genome		= $est2genome\n"}
if (defined($other_gff)) {$maker_variables{'other_gff'} = check_path($other_gff); print STDERR "other_gff		= $other_gff\n"}
if (defined($pred_flank)) {$maker_variables{'pred_flank'} = $pred_flank; print STDERR "pred_flank		= $pred_flank\n"}
if (defined($AED_threshold)) {$maker_variables{'AED_threshold'} = $AED_threshold; print STDERR "AED_threshold		= $AED_threshold\n"}
if (defined($min_protein,)) {$maker_variables{'min_protein'} = $min_protein; print STDERR "min_protein		= $min_protein\n"}
if (defined($alt_splice)) {$maker_variables{'alt_splice'} = $alt_splice; print STDERR "alt_splice		= $alt_splice\n"}
if (defined($map_forward,)) {$maker_variables{'map_forward,'} = $map_forward,; print STDERR "map_forward,		= $map_forward\n"}
if (defined($keep_preds)) {$maker_variables{'keep_preds'} = $keep_preds; print STDERR "keep_preds		= $keep_preds\n"}
if (defined($split_hit)) {$maker_variables{'split_hit'} = $split_hit; print STDERR "split_hit		= $split_hit\n"}
if (defined($softmask)) {$maker_variables{'softmask'} = $softmask; print STDERR "softmask		= $softmask\n"}
if (defined($single_exon)) {$maker_variables{'single_exon'} = $single_exon; print STDERR "single_exon		= $single_exon\n"}
if (defined($single_length)) {$maker_variables{'single_length'} = $single_length; print STDERR "single_length		= $single_length\n"}
if (defined($max_dna_len)) {$maker_variables{'max_dna_len'} = $max_dna_len; print STDERR "max_dna_len		= $max_dna_len\n"}
if (defined($min_contig)) {$maker_variables{'min_contig'} = $min_contig; print STDERR "min_contig		= $min_contig\n"}
if (defined($tmp)) {$maker_variables{'TMP'} = $tmp; print STDERR "TMP		= $tmp\n"}

print STDERR "\n";

if( !defined($est) && !defined($est_reads) && !defined($altest) && !defined($est_gff) ) { print STDERR "NON-FATAL ERROR: EST data is missing, for best results provide a file for at least one\n" }
if( !defined($protein) && !defined($protein_gff) ) { print STDERR "NON-FATAL ERROR: Protein data is missing, for best results provide a file for at least one\n" }
if( !defined($snaphmm) && !defined($gmhmm) && !defined($augustus_species) && !defined($est2genome) ) { print STDERR "NON-FATAL ERROR: Undefinied gene model prediction method, need to define at least one ab initio predictor or set est2genome to 1\n" }

print STDERR "\nPlease verify all the inputted variables were properly read in, only values list it above would be use into the analysis\n";
print STDERR "Analysis will be perform in the current folder, all output files and directories will be written here, do you wish to continue?? {Y/N} : ";
my $answer = <STDIN>;
if ($answer =~ /^n(?:o)?$/i) {
	die "Analysis terminated by user.\n\n";
} else {
	print STDERR "\n"
} #end else for output warning


#########################################
#                                       #
#        Preparing Maker CTL File       #
#                                       #
#########################################

system("maker -CTL");

open(MAKEROPTS, "maker_opts.ctl");
open(TEMPOPTS, ">maker_opts.temp");

print STDERR "Modifying the maker_opts.ctl file\n\n";

while(my $optsLine = <MAKEROPTS>) {

       	chomp($optsLine);

        if( ($optsLine =~ /^[[:alpha:]]/) ) {

		my @line = split ("#", $optsLine);

		foreach my $variable (keys(%maker_variables)) {

			if($optsLine =~ /^$variable=/) {
                         	print TEMPOPTS $variable, "= ", $maker_variables{$variable}," #",$line[1]; @line = (); last
                        } #end if match
		} #end foreach

		print TEMPOPTS join("#", @line);

	} else {

                print TEMPOPTS $optsLine;

        } #end else for modification of line

        print TEMPOPTS "\n";

} #end of while maker_opts.ctl

system("mv maker_opts.temp maker_opts.ctl");


#### Create output folders ####

system("mkdir maker_run_logs") unless -d "maker_run_logs";

system("mkdir qsubLogs") unless -d "qsubLogs";
system("mkdir qsubs") unless -d "qsubs";

system("mkdir GFF_files") unless -d "GFF_files";
system("mkdir maker_datastore_logs") unless -d "maker_datastore_logs";

system("mkdir maker_runs") unless -d 'maker_runs';

system("mkdir maker_runs/$maker_variables{'TMP'}") unless -d "maker_runs/".$maker_variables{'TMP'};

if(-d "flags") {

	system("rm -r flags");
	system("mkdir flags");

} else {

	system("mkdir flags");

} #end else for flags folder



#########################################
#                                       #
#           Run setup script            #
#                                       #
#########################################

if(!defined($qsubsFile)) {

	open(SETUP, ">qsubs/setup.sh");

##	print SETUP <<setuptext

#!/bin/bash

#\$ -S /bin/bash
#\$ -o $pwd/qsubLogs/setup.stdout
#\$ -e $pwd/qsubLogs/setup.stderr
#\$ -N maker_cluster_annotator_setup
#\$ -pe threaded 1

##cd $pwd
##source ~/.bash_profile
##source ~/.bashrc
##hostname
##Maker_qsub_annotator_setup.pl $inputFile $numseq $threads $again $force $tmp
##touch flags/setupflag

##setuptext
##;
##	system("qsub qsubs/setup.sh");

	system("Maker_qsub_annotator_setup.pl $inputFile $numseq $threads $again $force $tmp");


#	while (!(-e "flags/setupflag")) {

#		sleep 10;

#		print STDERR "waiting for setup to end\n";

#	} #end while

#	system("rm flags/setupflag");

	print STDERR "Setup of files completed\n\n";

	$qsubsFile = "qsub_list.txt";

} else {

	print STDERR "Qsubs list from $qsubsFile will be use\n\n";

} #end else if qsub file provided




#########################################
#                                       #
#       Running maker from qsubs        #
#                                       #
#########################################

open(QSUBS, $qsubsFile);

while(my $qsub = <QSUBS>) {

	system("qsub $qsub");

} #end

my $allDone = 'FALSE';

while($allDone eq 'FALSE') {

	$allDone = 'TRUE';

	system("qstat > qstat.log");

	open(QSTAT, "qstat.log");

	while(my $qstatLine = <QSTAT>) {

		if($qstatLine =~ m/.*maker.*/) {

			$allDone = 'FALSE';
			print STDERR "Waiting for maker\n";
			sleep 120;
			last;

		} #end if

	} #end while qstat reading

} #end while waiting system

system("rm qstat.log");

print STDERR "Maker Jobs completed\n\n";


#########################################
#                                       #
#        Verifying output stats         #
#                                       #
#########################################

opendir(DATASTORELOGS, "maker_datastore_logs");

while(my $datastore = readdir(DATASTORELOGS) ) {

	if ($datastore =~ /.*master_datastore_index.*/) {

		system("cat maker_datastore_logs/$datastore >> master_datastore_index.log");

	} #end if

} #end while datastore files


##open(CHECKER, ">qsubs/checker.sh");

##print CHECKER <<checkertext

#!/bin/bash

#\$ -S /bin/bash
#\$ -o $pwd/qsubLogs/checker.stdout
#\$ -e $pwd/qsubLogs/checker.stderr
#\$ -N maker_exit_status_checker
#\$ -pe threaded 1

##cd $pwd
##source ~/.bash_profile
##source ~/.bashrc
##hostname
##Runs_checker.pl
##touch flags/checkflag

##checkertext
##;

##system("qsub qsubs/checker.sh");

system("Runs_checker.pl");

##while (!(-e "flags/checkflag")) {

##	sleep 20;

##	print STDERR "waiting for checker to end\n";

##} #end while

print STDERR "Validation of the exit status completed\n\n";

##system("rm flags/checkflag");


#########################################
#                                       #
#         ReRun maker if asked          #
#                                       #
#########################################

if($retry eq "yes") {

	open(FAILEDRUNS, "failed_runs.log");

	while( my $failed = <FAILEDRUNS>) {

		chomp($failed);

		system("qsub qsubs/$failed.sh");

	} #end while rerun failed ones

	my $allFailedDone = 'FALSE';

	while($allFailedDone eq 'FALSE') {

		$allFailedDone = 'TRUE';

		system("qstat > qstat.log");

		open(QSTAT, "qstat.log");

		while(my $qstatLine = <QSTAT>) {

			if($qstatLine =~ m/.*maker.*/) {

				$allFailedDone = 'FALSE';
				print STDERR "Waiting for maker again\n";
				sleep 120;
				last;

			} #end if

		} #end while qstat reading

	} #end while waiting system

	system("rm qstat.log");

	print STDERR "reRun of failed files completed\n\n";

} #end if re run is need it


#########################################
#                                       #
#       Consolidate resulting gff       #
#                                       #
#########################################

##open(CONSOLIDATOR, ">qsubs/consolidator.sh");

##print CONSOLIDATOR <<consolidatortext

#!/bin/bash

#\$ -S /bin/bash
#\$ -o $pwd/qsubLogs/consolidator.stdout
#\$ -e $pwd/qsubLogs/consolidator.stderr
#\$ -N maker_annotation_consolidator
#\$ -pe threaded 1

##cd $pwd
##source ~/.bash_profile
##source ~/.bashrc
##hostname
##Annotation_consolidator.pl $gffprefix
##touch flags/consolidatorflag

##consolidatortext
##;

##system("qsub qsubs/consolidator.sh");

system("Annotation_consolidator.pl $gffprefix");

##while (!(-e "flags/consolidatorflag")) {

##	sleep 5;

##	print STDERR "waiting for consolidator to end\n";

##} #end while

##system("rm flags/consolidatorflag");

print STDERR "End with the consolidation step\n\n";

system("rm -r flags");

print STDERR "MAKER ANNOTATION COMPLETED\n\n";

print STDERR "For details in the output files please use -h option for extended help\n\n";

system("mv maker_opts.ctl $time.maker_opts.ctl");
system("mv master_datastore_index.log $time.master_datastore_index.log");

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
#         Check path sub           #
#                                  #
####################################

sub check_path {

	my ($path) = @_;

	if ($path eq "") { return $path;

	} elsif ($path =~ /^\//) { return $path;

        } else { return "../".$path } #end if absolute path

} #end of check_path sub


####################################
#                                  #
#          Get time sub            #
#                                  #
####################################

sub getTime {

	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	return "$hour-$minute-$second.$months[$month]-$dayOfMonth-$year";

} #end of getTime sub


####################################
#                                  #
#         Print Usage sub          #
#                                  #
####################################

sub usage {

print <<usagetext

Usage:
maker_automatic_annotator_v1.pl (-input <FASTA/TXT>|-qsubs <qsubs.txt>) -threads <INT> -numseq <INT>  other_options

Require inputs:

input|qsubs	

numseq		if input its .fasta

################### Input #########################

input		input file, can be a fasta file or a text file
		for fasta file, this would be splitted using numseq as value
		for text file, list of fasta files that would be use for analysis (file "FastaFiles.txt" generate from this same script
qsubs		qsubs list file (one qsub configuration file). If the script was run before inputting the qsub file will skip the setup 
		step and improve the system. Will also ensure that the same maker commands are run.
threads		number of slots to be requested per job
numseq		in case that input is fasta, this specify how many sequences should be use per subfasta file
retry		Should the program retry fail sequences (yes/no), defaul = no
queue		queue were the jobs will be launch, default it's to don't specify queue and therefore would be launch in the defaul for the scheduler

###### MAKER VARIABLES ######

again		Caculate all annotations and output files again even if no settings have changed. Does not delete old analyses.
force		Forces maker to delete old files before running again. This will require all blast analyses to be re-run.

As in the maker_opts file from maker 2.10

#-----EST Evidence (for best results provide a file for at least one)
est		non-redundant set of assembled ESTs in fasta format (classic EST analysis)
est_reads	unassembled nextgen mRNASeq in fasta format (not fully implemented)
altest		EST/cDNA sequence file in fasta format from an alternate organism
est_gff		EST evidence from an external gff3 file
altest_gff	Alternate organism EST evidence from a separate gff3 file

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein		protein sequence file in fasta format
protein_gff	protein homology evidence from an external gff3 file

#-----Repeat Masking (set values to blank to skip repeat masking)
model_org	select a model organism for RepBase masking in RepeatMasker
rmlib		provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein	provide a fasta file of transposable element proteins for RepeatRunner
rm_gff		repeat elements from an external GFF3 file

#-----Gene Prediction
snaphmm		SNAP HMM file
gmhmm		GeneMark HMM file
augustus_species Augustus gene prediction species model
pred_gff	ab-initio predictions from an external GFF3 file
model_gff	annotated gene models from an external GFF3 file (annotation pass-through)
est2genome	infer gene predictions directly from ESTs, 1 = yes, 0 = no, default = 0

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff	features to pass-through to final output from an extenal GFF3 file

#-----MAKER Behavior Options
max_dna_len	length for dividing up contigs into chunks (increases/decreases  memory usage), default = 10000
min_contig	skip genome contigs below this length (under 10kb are often useless), default = 1

pred_flank	flank for extending evidence clusters sent to gene predictors, default = 200
AED_threshold	Maximum Annotation Edit Distance allowed (bound by 0 and 1), default = 1
min_protein	require at least this many amino acids in predicted proteins, default = 0
alt_splice	Take extra steps to try and find alternative splicing, 1 = yes, 0 = no, default = 0
always_complete	force start and stop codon into every gene, 1 = yes, 0 = no, default = 0
map_forward	map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no, default = 0
keep_preds	Add unsupported gene prediction to final annotation set, 1 = yes, 0 = no, default = 0

split_hit	length for the splitting of hits (expected max intron size for evidence alignments), default = 1000
softmask	use soft-masked rather than hard-masked (seg filtering for wublast), default = 1
single_exon	consider single exon EST evidence when generating annotations, 1 = yes, 0 = no, default = 0
single_length	min length required for single exon ESTs if 'single_exon is enabled', default = 250

################### Output ########################

After running the program it will generate a set of files and folders to store all the results

FILES

If fasta file is inputed the system will generate the next files
	Sequence_index.txt	contains the location of each individual fasta sequence on the subFasta's
	FastaFiles.txt		contains a list of all the fasta's create it from the main fasta file. This file can be re-inputed in further runs to avoid re-splitting the main fasta.

Maker logs
	maker_threads.log	log file of when the threads were started and when they finish
	maker_commands.log	log file of the commands that were use to run maker (if you want to reRun a specific one, the command need to be lunch in maker_runs directory)

Maker control files
	maker_opts.ctl		opts file generate it by maker, this file its modified base in user inputted variables
	maker_bopts.ctl		bopts file, is not modified
	maker_exe.ctl		exe file, is not modified

Results logs
	master_datastore_index.log	concatenation of all the datastore_index from each individual run
	Failed_Sequences.txt		list of sequences that report not finished in datastore_index log

Results GFF
	*.maker.output.*.gff			resulting gff, create by source


DIRECTORIES

Fastas			fasta files generate it from the splitting of the main fasta (if a fasta file is inputed)
maker_run_logs		log of each run (piping of standard error for each maker iteration)
maker_runs		folder with the maker results
GFF_files		all the gff files generate it from maker

usagetext

} #end of usage sub


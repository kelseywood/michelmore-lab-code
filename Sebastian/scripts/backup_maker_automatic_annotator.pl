#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use Parallel::ForkManager;

### Working variables
my $inputFile;
my $threads;
my $numseq;
my $help;
my $retry = 'no';

### Maker variables
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
                'threads=i' => \$threads,
                'numseq=i' => \$numseq,
                'retry=s' => \$retry,

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

} elsif ( !defined($inputFile) || !defined($threads) || ( ($inputFile =~ /.*\\.(fa)|(fasta)$/) && !defined($numseq) ) ) {

	print "Missing require arguments\n";
        print "Usage:\n";
        print "wakawaka.pl -input <input FILE> -threads <Number of threads> -numseq <Number of Sequences per file> ...\n";
	print "Use --help for detail help\n";
        die "Please try again.\n"

} #end if usage

print STDERR "\nInputted variables\n\n";

print STDERR "input		= $inputFile\n";
if(defined($numseq)) {print STDERR "Num Sequences	= $numseq\n"}
print STDERR "Num threads	= $threads\n";

#### Hashing variables ####

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

print STDERR "\nPlease verify all the variables input were properly read, only values list it above would be change into the maker_opts file\n";
print STDERR "Analysis will be perform in the current folder, all output files and directories will be written here, do you wish to continue?? {Y/N} : ";
my $answer = <STDIN>;
if ($answer =~ /^n(?:o)?$/i) {
	die "Analysis terminated by user.\n\n";
} else {
	print STDERR "\n"
} #end else for output warning


#########################################
#                                       #
#        Splitting the Fasta File       #
#                                       #
#########################################

open(FASTALOG, ">FastaFiles.txt");
open(FASTALOCATION, ">Sequence_index.txt");
my @fastas;

my $outputFile = 1;
my $sequence = 0;

open(FASTA, $inputFile);

my @inputPath = split("/", $inputFile);
my $inputName = pop(@inputPath);

if ( $inputName =~ /.*\\.(fa)|(fasta)$/) {

	print STDERR "Inputed file it's a fasta file\nSplitting file into subsets\n";

	system("mkdir Fastas") unless -d "Fastas";

	open(SUBFASTA, ">Fastas/$inputName.$outputFile.fasta");
	print FASTALOG "Fastas/$inputName.$outputFile.fasta\n";
	push(@fastas, ["Fastas/$inputName.$outputFile.fasta", "$inputName.$outputFile.fasta"]);

	while (my $fastaLine = <FASTA>) {
	        if ((substr($fastaLine,0,1)) ne ">") { print SUBFASTA $fastaLine }

	        else {
	                if($sequence < $numseq) {
	                        print SUBFASTA $fastaLine;
				print FASTALOCATION  $inputName.$outputFile, ".fasta", "\t",substr($fastaLine,1,100);
	                        ++$sequence;

	                } else {
	                        close(SUBFASTA);
				print FASTALOCATION  $inputName.$outputFile, ".fasta", "\t",substr($fastaLine,1,100);
	                        ++$outputFile;
	                        $sequence = 1;
	                        open(SUBFASTA, ">Fastas/$inputName.$outputFile.fasta");
	                        print FASTALOG "Fastas/$inputName.$outputFile.fasta\n";
	                        push(@fastas, ["Fastas/$inputName.$outputFile.fasta", "$inputName.$outputFile.fasta"]);

	                        print SUBFASTA $fastaLine;
	                        ++$sequence;

	                } #end if new file its open
	        } #end if its header
	} #end while reading fasta file

} elsif ( $inputName =~ /.*\\.(txt)$/) {

	print STDERR "Inputed file it's a text file\nFilenames will be loaded into memory\n";

	while (my $fastaLine = <FASTA>) {

                chomp($fastaLine);

                my @path = split("/", $fastaLine);

                push(@fastas, [$fastaLine, pop(@path)]);

        } #end while input file

} else {
	print "Unrecognized extension in the input sequences file\n";
        die "Please verify your input file\n\n";
} #end else input is fasta or txt

close(FASTA);
close(SUBFASTA);

print STDERR "Done splitting the fasta file\n";
print STDERR "$outputFile files were created\n\n";


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


#########################################
#                                       #
#      Setting the run conditions       #
#                                       #
#########################################

system ("mkdir maker_run_logs") unless -d "maker_run_logs";
open(MAKERJOBSLOG, ">>maker_threads.log");
open(MAKERCOMMANDSLOG, ">maker_commands.log");

my $pm = new Parallel::ForkManager($threads);
##### Parallelizing Conditions #####
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print MAKERJOBSLOG "$ident finish with PID $pid, at ", getTime(), "\n"; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; print MAKERJOBSLOG "$ident started, pid: $pid, at ", getTime(),"\n"; } );
$pm->run_on_wait( sub {print STDERR "."}, 10 );


#########################################
#                                       #
#        Running Maker in Fastas        #
#                                       #
#########################################

my $numFiles = scalar(@fastas);
my $m = int $numFiles/100;
#When less than 100 files are analyzed m need to be corrected
if($m == 0) {$m = 1};
my $analizedfile = 0;

print STDERR "Percentage advance\n";

system("mkdir maker_runs") unless -d 'maker_runs';

$tmp = "maker_runs/".$maker_variables{'TMP'};

system("mkdir $tmp") unless -d $tmp;

foreach my $fasta (@fastas) {

	print STDERR int(100 * $analizedfile / $numFiles), " " if $analizedfile++ % $m == 0;

	print MAKERCOMMANDSLOG "maker -nodatastore -cpus 1 -genome ../@$fasta[0] ../maker_opts.ctl ../maker_exe.ctl ../maker_bopts.ctl 2>> ../maker_run_logs/@$fasta[1].log\n";

#        my $pid = $pm->start(@$fasta[1]) and next;

#	chdir("maker_runs");

#	system("maker -nodatastore -cpus 1 -genome ../@$fasta[0] ../maker_opts.ctl ../maker_exe.ctl ../maker_bopts.ctl 2>> ../maker_run_logs/@$fasta[1].log");

#        $pm->finish;

} #end of running maker

$pm -> wait_all_children;

print STDERR "\n\n";

#########################################
#                                       #
#      Consolidating the Results        #
#                                       #
#########################################

##### Find GFF's #####

print STDERR "Extracting gff out of datastore\n";

system("mkdir GFF_files") unless -d "GFF_files";

my @resultGFFs;
my $makerOutputDir = "maker_runs";

opendir(DIR, $makerOutputDir);
while (my $makerRun = readdir(DIR)) { find_gff($makerRun, $makerOutputDir) } #end foreach directory
closedir(DIR);

print STDERR scalar(@resultGFFs), " number of gff's were generated\n\n";

open(FAILED, ">>Failed_Sequences.txt");

open(DATASTORELOG, "master_datastore_index.log");

my %filesToRun;

while( my $datastoreLogLine = <DATASTORELOG>) {

       	my @logLine = split("\t", $datastoreLogLine);

	if( !( ($datastoreLogLine =~ /.*FINISHED$/ ) || ($datastoreLogLine =~ /.*STARTED$/ ) ) && defined($logLine[1]) ) {

		print FAILED $datastoreLogLine;

                my @failedFile = split("/", $logLine[1]);

                $failedFile[0] =~ s/_datastore//;

		$filesToRun{$failedFile[0]} = "";

        } #end if is not finished or started

} #end while for datastorelog

if(%filesToRun) {

	print STDERR "Some sequences failed\n\n";

        if($retry eq 'yes') {
        	print STDERR "Sequences will be rerun\n";

		foreach my $fileToRun (keys (%filesToRun) ) {

			print STDERR "Re-analyzing $fileToRun\n";

#                        my $pid = $pm->start($fileToRun) and next;

#			chdir("maker_runs");

#			system("maker -nodatastore -cpus 1 -genome ../@$fasta[0] ../maker_opts.ctl ../maker_exe.ctl ../maker_bopts.ctl 2>> ../maker_run_logs/@$fasta[1].log");

#		        $pm->finish;

                } #end foreach file to rerun

#		$pm -> wait_all_children;

		print STDERR "\n";

        } else {

        	print STDERR "Sequences will be skipped\n\n";

        } #end else retry yes

} else {

	print STDERR "Congratulations all sequences were completed\n\n";

} # end of checking if failed sequences

##### Generate complete GFF #####

print STDERR "Consolidating annotation into by source GFF\n";

foreach my $gff (@resultGFFs) {

        my %annotation;

	open(GFF, $gff);

        while(my $gffLine = <GFF>) {

        	if($gffLine =~ /^[[:alnum:]]/) {

                	my @fields = split("\t", $gffLine);

                        if($fields[1] =~ /^[[:alnum:]]/) {

                        	push(@{$annotation{$fields[1]}}, $gffLine);

                        } #end if source starts with a letter

                } elsif ($gffLine =~ /^>/) {

                	last;

                } #end elsif line

        } #end while

        close(GFF);

        foreach my $source (keys(%annotation)) {

        	open(SOURCE, ">>$inputName.maker.output.$source.gff");

                foreach my $line (@{$annotation{$source}}) { print SOURCE $line }

                close(SOURCE);

        } #end foreach

} #end foreach gff files

### Add separator lines to growing files

print MAKERJOBSLOG "\n\nEND OF RUN\nOn ", getTime(),"\n\n";
print FAILED "\n\nEND OF RUN\nOn ", getTime(),"\n\n";
open(DATASTORE, '>>master_datastore_index.log');
print DATASTORE "\n\nEND OF RUN\nOn ", getTime(),"\n\n";

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

        if ($path =~ /^\//) { return $path;

        } else { return "../".$path } #end if absolute path

} #end of check_path sub

####################################
#                                  #
#         Find GFF file sub        #
#                                  #
####################################

sub find_gff {

	my ($makerRun, $makerOutputDir) = @_;

	#Find output directories
	if($makerRun =~ /output$/ && -d "$makerOutputDir/$makerRun") {

		#open subdir and make new path
		opendir(RUNSRESULTS, "$makerOutputDir/$makerRun");

		while(my $datastore = readdir(RUNSRESULTS) ) {

			#Find the datastore
			if($datastore =~ /datastore$/) {

				opendir(DATASTORE, "$makerOutputDir/$makerRun/$datastore");

				#if is a directory add to dir list or add to gff list if .gff extension

				while(my $sequence = readdir(DATASTORE) ) {

					if( -d "$makerOutputDir/$makerRun/$datastore/$sequence") {

						opendir(SEQUENCERESULT, "$makerOutputDir/$makerRun/$datastore/$sequence");

						while( my $result = readdir(SEQUENCERESULT) ) {

							my $extension = substr($result,-3);

							if( $extension eq "gff") {

								system("cp $makerOutputDir/$makerRun/$datastore/$sequence/$result GFF_files/$result");
			                                	push(@resultGFFs, "GFF_files/$result");

							} #end if for gff files
						} #end while reading result dir
					} #end if it is a sequence directory
				} #end reading datastore directory
			} elsif ($datastore =~ /master_datastore_index.log$/) {

				system("cat $makerOutputDir/$makerRun/$datastore >> master_datastore_index.log");

			} #end if is datastore or log
		} #end while maker run directory
	} #end if is an output directory
} #end find_gff sub


sub getTime {

	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	return "$hour:$minute:$second, $months[$month] $dayOfMonth, $year";

} #end of getTime sub


####################################
#                                  #
#         Print Usage sub          #
#                                  #
####################################

sub usage {

print <<usagetext

retry		Should the program retry fail sequences (yes/no), defaul = no


###### MAKER VARIABLES ######

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

usagetext

} #end of usage sub

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_f $opt_d $opt_r $opt_t $opt_k $opt_M $opt_G $opt_w $opt_v);
getopts('hf:d:r:t:k:G:w:v');

my $mismatch = 0.05;
my $seed = 19;
my $band = 100;
my $verbosity = 3;
my $cores = 1;

my $files = $opt_f if $opt_f;
my $outDir = $opt_d if $opt_d;
my $reference = $opt_r if $opt_r;

$cores = $opt_t if $opt_t;
$seed = $opt_k if $opt_k;
$band = $opt_w if $opt_w;
$verbosity = $opt_v if $opt_v;

if( (defined($opt_h)) ) {

	print "BWA-MEM_batch_paired_aligner.pl HELP\n\n";
	print "BWA-MEM_batch_paired_aligner.pl -h -f FastqFiles.txt -r reference.fasta -d OutputDir <other options>\n";
	print "-h output help\n\n";

	print "INPUT\n\n";

	print "REQUIRED VARIABLES\n";
	print "-f tab-delimited file with the list of fastq files to align (each line with two columns, forward and reverse sense)\n";
	print "-d output directory for alignment files\n";
	print "-r Name of the indexed reference (prefix)\n\n";

	print "OPTIONAL VARIABLES\n";
	print "-t INT number of cores to use (Default = 1)\n";
	print "-k INT minimum seed length (Default = 19)\n";
	print "-w INT band width for banded alignment, maximum gap size (Default = 100)\n";
	print "-v INT verbose level: 1=error, 2=warning, 3=message, 4+=debugging (Default = 3)\n";

	print "Example for FastqFiles.txt\n";
	print "Read_1_forward.fq	Read_2_reverse.fq\n";
	print "Read_2_forward.fq	Read_3_reverse.fq\n";
	print "Read_2_forward.fq	Read_3_reverse.fq\n\n";

	print "OUTPUT\n";
	print "All results are stored in OutputDir folder\n";
	print "Three files are created for each pair of reads, the results file get their name from the forward file\n";
	print "		Forward.sam	SAM alignment output from bwa\n";
	print "		Forward.sam.log	Log of the bwa alignment (useful in case of erros)\n";
	print "		Forward.bam	BAM format file of the alignment\n\n";

	print "Help menu trigger\n\n";
	exit;

} #end if


if( !(defined($files)) || !(defined($reference)) || !(defined($outDir)) ) {

	print "BWA_bath_alginer.pl -h -f fileList -f FastqFiles.txt -r reference.fasta -d OutputDir <other options>\n";
	print "-h output help\n";

	print "REQUIRED VARIABLES\n";
	print "-f tab-delimited file with the list of fastq files to align (each line with two columns, forward and reverse sense)\n";
	print "-d output directory for alignment files\n";
	print "-r Name of the indexed reference (prefix)\n";

	print "Other BWA variables\n";
	print "-t INT number of cores to use (Default = 1)\n";
	print "-k INT minimum seed length (Default = 19)\n";
	print "-w INT band width for banded alignment, maximum gap size (Default = 100)\n";
	print "-v INT verbose level: 1=error, 2=warning, 3=message, 4+=debugging (Default = 3)\n";

	print "Missing Arguments.\n";
	print "Please input all the arguments need it\n";
	exit;

} #end if

open(INPUT, $files) or die "Can't open $files\n";

system("mkdir $outDir") unless -d $outDir;

open(LOG, ">$outDir/$files.aligner.log.txt");

print LOG "Input fastq files: ", $files, "\n";
print LOG "Reference index: ", $reference, "\n";
print LOG "Output dir: ", $outDir, "\n";
print LOG "Number of threads to use: ", $cores, "\n";
print LOG "Seed length: ", $seed, "\n";
print LOG "Bad with (maximun gap size_: ", $band, "\n";
print LOG "Verbosity level: ", $verbosity, "\n";

open(SAMS, ">$outDir/$files.SAMlist.txt");
open(TABS, ">$outDir/$files.TABlist.txt");

while(my $line = <INPUT>) {

	chomp($line);

	my @files = split ("\t", $line);

	print STDERR "\nProcessing $line files.\n";
	print  LOG "\nProcessing $files[0] files.\n";

	if(!defined($files[1])) { print "$line doesn't contain two files separated by tab, please check\n"; next}

	my @arrayName1 = split("/", $files[0]);
	my $filename1 = pop(@arrayName1);

	my @arrayName2 = split("/", $files[1]);
	my $filename2 = pop(@arrayName2);

	#Align the reads
	print LOG "/home/sreyesch/MichelmoreBin/bwa-0.7.4/bwa mem -M -t $cores -w $band -v $verbosity $reference $files[0] $files[1] > $outDir/$filename1.paired.sam 2> $outDir/$filename1.paired.sam.log\n";

	if( system("/home/sreyesch/MichelmoreBin/bwa-0.7.4/bwa mem -M -t $cores -w $band -v $verbosity $reference $files[0] $files[1] > $outDir/$filename1.paired.sam 2> $outDir/$filename1.paired.sam.log") != 0 ) {

		print STDERR "ERROR: Couldn't align reads for $files[0] and $files[1].\n";
		print LOG "ERROR: Couldn't align reads for $files[0] and $files[1].\n";

	} else {

		print SAMS "$files[1].paired.sam\n";

	} #end running bwa


	#Filter SAM file
	print LOG "/home/sreyesch/scripts/BWA_Analysis/SAM_filterer_v2.pl $outDir/$filename1.paired.sam\n";

	if( system("/home/sreyesch/scripts/BWA_Analysis/SAM_filterer_v2.pl $outDir/$filename1.paired.sam") != 0 ) {

		print STDERR "ERROR: Couldn't filter SAM for $files[0] and $files[1].\n";
		print LOG "ERROR: Couldn't filter SAM for $files[0] and $files[1].\n";

	} 

	#Generate TAB file
	print LOG "perl /home/sreyesch/MichelmoreBin/SSPACE-BASIC-2.0_linux-x86_64/tools/sam_bam2tab.pl $outDir/$filename1.paired.sam.filterMappedPaired.sam /1 /2 $outDir/$filename1.paired.sam.filterMappedPaired.tab\n";

	if( system("perl /home/sreyesch/MichelmoreBin/SSPACE-BASIC-2.0_linux-x86_64/tools/sam_bam2tab.pl $outDir/$filename1.paired.sam.filterMappedPaired.sam /1 /2 $outDir/$filename1.paired.sam.filterMappedPaired.tab") != 0 ) {

		print STDERR "ERROR: Couldn't generate TAB for $files[0] and $files[1].\n";
		print LOG "ERROR: Couldn't generate TAB for $files[0] and $files[1].\n";

	} else {

		print TABS "$filename1.paired.sam.filterMappedPaired.tab\n";

	} #end generating TAB files

} #end while for input files

exit;


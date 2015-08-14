#!/usr/bin/perl

######################################################
#                                                    #
#             Maker qsub Annotator setup             #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to prepare the qsub files to run maker
# Part of the main pipeline "Maker_qsub_annotator.pl"

### Difference from v1
#
# Script designed to generate qsub files to run on the Lonestar HPC
#

use strict;
use warnings;
use Cwd;

if(@ARGV < 4) {

	print "Missing arguments\n";
	print "Usage: Maker_qsub_annotator_setup.pl <inputfile> <numseq> <amountseq> <threads> <again> <force> <queue> <account> <email> <tmp>\n";
	print "For detailed explanation on required arguments please direct to Maker_qsub_prepareRuns.pl or Maker_qsub_annotator_v2.pl\n";
	die "Please input all arguments\n";

} #end if missing arguments


## Retrieving input arguments

my $inputFile = $ARGV[0];
my $numseq = $ARGV[1];
my $amountseq = $ARGV[2];
my $threads = $ARGV[3];
my $again = $ARGV[4];
my $force = $ARGV[5];
my $queue = $ARGV[6];
my $account = $ARGV[7];
my $email =  $ARGV[8];
my $gffprefix = $ARGV[9];
my $tmp = $ARGV[10];


my $megabasesAllowed = $amountseq*1000000;

my $pwd = getcwd;

my $username = getpwuid( $< );

#########################################
#                                       #
#         Detecting Input File          #
#                                       #
#########################################

my %fastas;

my @inputPath = split("/", $inputFile);
my $inputName = pop(@inputPath);

if ( $inputName =~ /.*\\.((f|F)(a|A))|((f|F)(a|A)(s|S)(t|T)(a|A))$/) {

	open(FASTA, $inputFile) or die "FATAL ERROR: Can't open file $inputFile\n\n";

	open(FASTALOCATION, ">Sequence_index.txt");

	open(FASTALOG, ">FastaFiles.txt");

	my $outputFile = 1;
	my $sequenceCount = 0;
	my $sequenceSum = 0;

	print STDERR "Inputed file it's a fasta file\nSplitting file into subsets\n";

	system("mkdir Fastas") unless -d "Fastas";

	open(SUBFASTA, ">Fastas/$inputName.$outputFile.fasta");
        my $subFastaName = $inputName.".".$outputFile;
	print FASTALOG "$subFastaName\t$subFastaName.fasta\tFastas/$subFastaName.fasta\n";
     	@{$fastas{$subFastaName}} = ("Fastas/$subFastaName.fasta", "$subFastaName.fasta");

	while (my $fastaLine = <FASTA>) {
	        if ( !($fastaLine =~ m/^>/) ) { print SUBFASTA $fastaLine; $sequenceSum += length($fastaLine)-1 }

	        else {

			print $sequenceSum, "\t";

	                if( ($sequenceCount < $numseq) && ($sequenceSum < $megabasesAllowed) ) {

	                        print SUBFASTA $fastaLine;
				print FASTALOCATION  $inputName.$outputFile, ".fasta", "\t",substr($fastaLine,1,100);
	                        ++$sequenceCount;

	                } else {
	                        close(SUBFASTA);
				print FASTALOCATION  $inputName.$outputFile, ".fasta", "\t",substr($fastaLine,1,100);
	                        ++$outputFile;
	                        $sequenceCount = 0;
				$sequenceSum = 0;
	                        open(SUBFASTA, ">Fastas/$inputName.$outputFile.fasta");
			        $subFastaName = $inputName.".".$outputFile;
				print FASTALOG "$subFastaName\t$subFastaName.fasta\tFastas/$subFastaName.fasta\n";
			     	@{$fastas{$subFastaName}} = ("Fastas/$subFastaName.fasta", "$subFastaName.fasta");

	                        print SUBFASTA $fastaLine;
	                        ++$sequenceCount;

	                } #end if new file its open
	        } #end if its header
	} #end while reading fasta file


	close(SUBFASTA);

	print STDERR "Done splitting the fasta file\n";
	print STDERR "$outputFile files were created\n\n";

} elsif ( $inputName =~ /.txt$/) {

	open(TXT, $inputFile) or die "FATAL ERROR: Can't open file $inputFile\n\n";

	print STDERR "Inputed file it's a text file\nFilenames will be loaded into memory\n";

	my $numFiles = 0;

	while (my $fastaLine = <TXT>) {

                chomp($fastaLine);

                my @fileData = split("\t", $fastaLine);

	     	@{$fastas{$fileData[0]}} = ($fileData[2],$fileData[1]);

		++$numFiles;

        } #end while input file

	print STDERR "$numFiles filenames were loaded\n\n";

} else {
	print "Unrecognized extension in the input sequences file\n";
        die "Please verify your input file\n\n";
} #end else input is fasta or txt

#########################################
#                                       #
#         Preparing qsub sh Files       #
#                                       #
#########################################

if ($again eq "YES") {$again = "-again"} else {$again = ""}
if ($force eq "YES") {$force = "-force"} else {$force = ""}

system("mkdir qsubs") unless -d "qsubs";

open(QSUBS, ">qsub_list.txt");

foreach my $fasta (keys (sort {$a <=> $b} %fastas) ) {

	print QSUBS "qsubs/$fasta.sh\n";

	open(QSUB, ">qsubs/$fasta.sh");

	print QSUB qq{
#!/bin/bash
#submit_maker.sh
                           # The following are options passed to qsub
#\$ -V                      # Inherit the submission environment
#\$ -cwd                    # Start job in submission directory
#\$ -e $fasta.stdout.e\$JOB_ID   # Combine stderr and stdout
#\$ -o $fasta.stderr.o\$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#\$ -q $queue               # Queue name normal
#\$ -m bes                  # Email at Begin and End of job or if suspended
#\$ -A $account       # Name of account to charge this job to
};

	if ($email ne "") {print QSUB "#\$ -M $email      # E-mail address (change to your e-mail)"}

	print QSUB <<qsubtext2

#\$ -N makerjob         # The name of your job
#\$ -pe $threads            # Specifies the number of CPU to use
 
# The maker command is called using ibrun, which implements MPI
# Add additional maker options to this line if needed

source \$HOME/.bash_profile

mkdir \$SCRATCH/$gffprefix
mkdir \$SCRATCH/$gffprefix/maker_runs
mkdir \$SCRATCH/$gffprefix/maker_runs/$tmp

cd \$SCRATCH/$gffprefix/maker_runs

hostname > \$SCRATCH/$tmp/$fasta.hostname.txt

date > \$SCRATCH/$tmp/$fasta.date.txt
echo "$fasta" "started" > \$SCRATCH/$tmp/$fasta.name.txt
paste \$SCRATCH/$tmp/$fasta.name.txt \$SCRATCH/$tmp/$fasta.hostname.txt \$SCRATCH/$tmp/$fasta.date.txt >> \$SCRATCH/$gffprefix/maker_threads.log

ibrun maker -fix_nucleotides -nodatastore $again $force -genome  \$SCRATCH/$gffprefix/@{$fastas{$fasta}}[0] \$SCRATCH/$gffprefix/maker_opts.ctl \$SCRATCH/$gffprefix/maker_exe.ctl \$SCRATCH/$gffprefix/maker_bopts.ctl 2>> \$SCRATCH/$gffprefix/maker_run_logs/@{$fastas{$fasta}}[1].log
rm -r $fasta.maker.output/mpi_blastdb/
rm $fasta.maker.output/$fasta.db
find_gff_v2.pl $fasta.maker.output \$SCRATCH/$gffprefix YES

tar -zcf $fasta.maker.output.tar.gz $fasta.maker.output

date > \$SCRATCH/$tmp/$fasta.date.txt
echo "$fasta" "finished" > \$SCRATCH/$tmp/$fasta.name.txt
paste \$SCRATCH/$tmp/$fasta.name.txt \$SCRATCH/$tmp/$fasta.hostname.txt \$SCRATCH/$tmp/$fasta.date.txt >> \$SCRATCH/$gffprefix/maker_threads.log
rm \$SCRATCH/$tmp/$fasta.hostname.txt \$SCRATCH/$tmp/$fasta.date.txt \$SCRATCH/$tmp/$fasta.name.txt

qsubtext2

} #end of preparing sh files

system("split -d 50 qsub_list.txt qsub_list.batch");

exit;

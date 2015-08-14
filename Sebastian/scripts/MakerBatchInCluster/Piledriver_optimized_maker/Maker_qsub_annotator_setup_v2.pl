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
my $gffprefix = $ARGV[6];
my $workingDir = $ARGV[7];


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

#			print $sequenceSum, "\t";

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

foreach my $fasta (sort keys (%fastas) ) {

	print QSUBS "qsubs/$fasta.sh\n";

	open(QSUB, ">qsubs/$fasta.sh");

	print QSUB <<qsubtext2
#!/bin/bash

#\$ -S /bin/bash
#\$ -o $workingDir/qsubLogs/$fasta.stdout
#\$ -e $workingDir/qsubLogs/$fasta.stderr
#\$ -N makerjob
#\$ -pe threaded $threads
#\$ -l piledriver
 
source \$HOME/.bash_profile

cd $workingDir/maker_runs

hostname > $workingDir/maker_runs/tmp/$fasta.hostname.txt

date > $workingDir/maker_runs/tmp/$fasta.date.txt
echo "$fasta" "started" > $workingDir/maker_runs/tmp/$fasta.name.txt
paste $workingDir/maker_runs/tmp/$fasta.name.txt $workingDir/maker_runs/tmp/$fasta.hostname.txt $workingDir/maker_runs/tmp/$fasta.date.txt >> $workingDir/maker_threads.log

maker -fix_nucleotides -nodatastore $again $force -genome  $workingDir/@{$fastas{$fasta}}[0] $workingDir/maker_opts.ctl $workingDir/maker_exe.ctl $workingDir/maker_bopts.ctl 2>> $workingDir/maker_run_logs/@{$fastas{$fasta}}[1].log
rm -r $fasta.maker.output/mpi_blastdb/
rm $fasta.maker.output/$fasta.db
find_gff_v2.pl $fasta.maker.output $workingDir YES

tar -zcf $fasta.maker.output.tar.gz $fasta.maker.output

date > $workingDir/maker_runs/tmp/$fasta.date.txt
echo "$fasta" "finished" > $workingDir/maker_runs/tmp/$fasta.name.txt
paste $workingDir/maker_runs/tmp/$fasta.name.txt $workingDir/maker_runs/tmp/$fasta.hostname.txt $workingDir/maker_runs/tmp/$fasta.date.txt >> $workingDir/maker_threads.log
rm $workingDir/maker_runs/tmp/$fasta.hostname.txt $workingDir/maker_runs/tmp/$fasta.date.txt $workingDir/maker_runs/tmp/$fasta.name.txt

qsubtext2

} #end of preparing sh files

system("split -d -l 50 qsub_list.txt qsub_list.batch");

exit;

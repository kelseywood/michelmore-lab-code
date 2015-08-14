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

use strict;
use warnings;
use Cwd;

if(@ARGV < 6) {

	print "Missing arguments\n";
	print "Usage: Maker_qsub_annotator_setup.pl <inputfile> <numseq> <amountsqe> <threads> <again> <force> <tmp>\n";
	die "Please input all arguments\n";

} #end if missing arguments


## Retrieving input arguments

my $inputFile = $ARGV[0];
my $numseq = $ARGV[1];
my $amountseq = $ARGV[2];
my $threads = $ARGV[3];
my $again = $ARGV[4];
my $force = $ARGV[5];
my $tmp = $ARGV[6];


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
	        if ((substr($fastaLine,0,1)) ne ">") { print SUBFASTA $fastaLine; $sequenceSum += length($fastaLine)-1 }

	        else {
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

foreach my $fasta (sort {$a <=> $b} keys %fastas) {

	print QSUBS "qsubs/$fasta.sh\n";

	open(QSUB, ">qsubs/$fasta.sh");

	print QSUB qq{
#!/bin/bash

#\$ -S /bin/bash
#\$ -o $pwd/qsubLogs/$fasta.stdout
#\$ -e $pwd/qsubLogs/$fasta.stderr
#\$ -N makerjob
#\$ -pe threaded $threads

mkdir /state/partition1/$username
mkdir /state/partition1/$username/maker_runs

cd /state/partition1/$username/maker_runs
source ~/.bash_profile
source ~/.profile
source ~/.bashrc

mkdir $tmp

hostname > $tmp/$fasta.hostname.txt

date > $tmp/$fasta.date.txt
echo "$fasta" "started" > $tmp/$fasta.name.txt
paste $tmp/$fasta.name.txt $tmp/$fasta.hostname.txt $tmp/$fasta.date.txt >> $pwd/maker_threads.log

};

	if(-f "$pwd/maker_runs/$fasta.maker.output.tar.gz") {
		print QSUB "cp $pwd/maker_runs/$fasta.maker.output.tar.gz .\n";
		print QSUB "tar -xvzf $fasta.maker.output.tar.gz\n";
	} #end if maker has been run

	print QSUB <<qsubtext2

cp $pwd/@{$fastas{$fasta}}[0] $fasta.fasta

maker -nodatastore -cpus $threads $again $force -genome $pwd/@{$fastas{$fasta}}[0] $pwd/maker_opts.ctl $pwd/maker_exe.ctl $pwd/maker_bopts.ctl 2>> $pwd/maker_run_logs/@{$fastas{$fasta}}[1].log
rm -r $fasta.maker.output/mpi_blastdb/
rm $fasta.maker.output/$fasta.db
find_gff.pl $fasta.maker.output $pwd YES

tar -zcf $fasta.maker.output.tar.gz $fasta.maker.output
mv $fasta.maker.output.tar.gz $pwd/maker_runs/.
rm -r $fasta.maker.output
rm $fasta.fasta

date > $tmp/$fasta.date.txt
echo "$fasta" "finished" > $tmp/$fasta.name.txt
paste $tmp/$fasta.name.txt $tmp/$fasta.hostname.txt $tmp/$fasta.date.txt >> $pwd/maker_threads.log
rm $tmp/$fasta.hostname.txt $tmp/$fasta.date.txt $tmp/$fasta.name.txt

qsubtext2

} #end of preparing sh files

exit;

#!/usr/bin/perl

######################################################
#                                                    #
#         Exonerate_batchCluster_submitter           #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to run the exonerate aligner through qsubs jobs in a cluster

use strict;
use warnings;
use Cwd;
use File::Basename;

## Retrieving input arguments

if(@ARGV < 4) {
	print "Missing arguments\n";
	die "Missing arguments\n\n";

} #end if


my $inputDir = $ARGV[0];
my $target = $ARGV[1];
my $exonerateOptions = $ARGV[2];
my $prefix = $ARGV[3];

my $foldername = basename($inputDir);

my $jobID = "$target-$foldername-$prefix";

my $pwd = getcwd;

#########################################
#                                       #
#         Preparing qsub sh Files       #
#                                       #
#########################################

system("mkdir qsubsFolder-$jobID") unless -d "qsubsFolder-$jobID";
system("mkdir tmp") unless -d "tmp";
if(-f "$jobID.exonerate.log") {system("rm $jobID.exonerate.log")}
if(-f "$jobID.exonerate.output") {system("rm $jobID.exonerate.output")}

if(-d "alignments-$jobID") {system("rm alignments-$jobID/*")} else {system("mkdir alignments-$jobID")}

my @qsubs;
my @outputs;

opendir(DIR, $inputDir);

while (my $fasta = readdir(DIR)) {

	if( $fasta eq ".." || $fasta eq ".") { next }

	push(@qsubs, "qsubsFolder-$jobID/$fasta.sh");
	push(@outputs, "alignments-$jobID/$fasta.exonerate.output");

	open(QSUB, ">qsubsFolder-$jobID/$fasta.sh") or die "Couldn't open file qsubsFolder-$jobID/$jobID.sh\n\n";

	print QSUB <<qsubtext

#!/bin/bash

#\$ -S /bin/bash
#\$ -e $pwd/qsubsFolder-$jobID/$jobID-$fasta.stderr
#\$ -N exonerate

cd $pwd
source ~/.bash_profile
source ~/.bashrc

hostname > tmp/$jobID-$fasta.hostname.txt

date > tmp/$jobID-$fasta.date.txt
echo "$fasta" "started" > tmp/$jobID-$fasta.name.txt
paste tmp/$jobID-$fasta.name.txt tmp/$jobID-$fasta.hostname.txt tmp/$jobID-$fasta.date.txt >> $jobID.exonerate.log

exonerate $exonerateOptions $inputDir/$fasta $target > alignments-$jobID/$fasta.exonerate.output

date > tmp/$jobID-$fasta.date.txt
echo "$fasta" "finished" > tmp/$jobID-$fasta.name.txt
paste tmp/$jobID-$fasta.name.txt tmp/$jobID-$fasta.hostname.txt tmp/$jobID-$fasta.date.txt >> $jobID.exonerate.log
rm tmp/$jobID-$fasta.hostname.txt tmp/$jobID-$fasta.date.txt tmp/$jobID-$fasta.name.txt 

qsubtext

} #end of preparing sh files

print "Qsubs files are ready\n\n";

foreach my $qsub (@qsubs) {

	system("qsub -q all.q $qsub >> $jobID.exonerate.log");

} #end foreach qsub


my $allDone = 'FALSE';

while($allDone eq 'FALSE') {

	$allDone = 'TRUE';

	system("qstat > qstat.log");

	open(QSTAT, "qstat.log");

	while(my $qstatLine = <QSTAT>) {

		if($qstatLine =~ m/.*exonerate.*/) {

			$allDone = 'FALSE';
			print STDERR "Waiting for exonerate jobs\n";
			sleep 60;
			last;

		} #end if

	} #end while qstat reading

} #end while waiting system

system("rm qstat.log");


print "Finish the alignment process\n\n";

foreach my $result (@outputs) {

	system("cat $result >> $jobID.exonerate.output");

} #end catting results

print "Finish merging the results\n\n";

system("rm -r tmp");
system("rm -r alignments-$jobID");

exit;

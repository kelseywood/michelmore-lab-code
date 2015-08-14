#!/usr/bin/perl

######################################################
#                                                    #
#         Interproscan inCluster qsub writer         #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was written to generate qsub files to run interproscan in a cluster using qsub
# given a directory of fasta files
# using the following options '-f GFF3,XML,HTML,TSV -goterms -iprlookup -ms 50 -t p'

use strict;
use warnings;
use Cwd;

if(@ARGV < 3) {

	print "Missing arguments\n";
	print "Usage: Interproscan_inCluster_qsubWriter.pl <listoffastafiles> <threads> <outputPrefix> <temp>\n";
	print "Require arguments:\n";
	print "listoffastafiles\n";
	print "threads\n";
	print "outputPrefix\n";
	print "\nOptional arguments:\n";
	print "temp\n";
	die "Please input all arguments\n";

} #end if missing arguments


## Retrieving input arguments

my $listoffastafiles = $ARGV[0];
my $threads = $ARGV[1];
my $prefix = $ARGV[2];

my $temp = "/state/partition1/";

if(defined($ARGV[3]) ) { $temp = $ARGV[3] }

my $pwd = getcwd;

my $username = getpwuid( $< );

#########################################
#                                       #
#         Preparing qsub sh Files       #
#                                       #
#########################################

system("mkdir qsubs") unless -d "qsubs";
system("mkdir qsubLogs") unless -d "qsubLogs";
system("mkdir interproscanLogs") unless -d "interproscanLogs";
system("mkdir interproscanOutput") unless -d "interproscanOutput";

open(QSUBS, ">qsub_list.txt");

open(FASTAS, $listoffastafiles);

while( my $fasta = <FASTAS>) {

	chomp($fasta);

	my $fastaPath;

	if ($fasta =~ /^\//) {

		$fastaPath = $fasta;

	} else {

		$fastaPath = "$pwd/$fasta";

	} #end checking path

	my @inputPath = split("/", $fasta);
	my $fastaName = pop(@inputPath);

	$fastaName =~ s/\.[^.]+$//;

	print QSUBS "qsubs/$fastaName.sh\n";

	open(QSUB, ">qsubs/$fastaName.sh");

# for new nodes "#\$ -l piledriver"

	print QSUB qq{
#!/bin/bash

#\$ -S /bin/bash
#\$ -o $pwd/qsubLogs/$fastaName.qsub.stdout
#\$ -e $pwd/qsubLogs/$fastaName.qsub.stderr
#\$ -N interproscan
#\$ -pe threaded $threads

mkdir $temp/$username
mkdir $temp/$username/interproscan

cd $temp/$username/interproscan
source ~/.bash_profile
source ~/.profile
source ~/.bashrc

mkdir $fastaName
cd $fastaName

hostname > tmp.$fastaName.hostname.txt

date > tmp.$fastaName.date.txt
echo "$fastaName" "started" > tmp.$fastaName.name.txt
paste tmp.$fastaName.name.txt tmp.$fastaName.hostname.txt tmp.$fastaName.date.txt >> $pwd/interproscan_threads.log

interproscan.sh -b $temp/$username/interproscan/$fastaName/$fastaName -f GFF3,XML,HTML,TSV -goterms -iprlookup -ms 50 -t p -i $fastaPath -T $temp/$username/interproscan/$fastaName/ > $pwd/interproscanLogs/$fastaName.interproscan.log 2> $pwd/interproscanLogs/$fastaName.interproscan.error.log

date > tmp.$fastaName.date.txt
echo "$fastaName" "finished" > tmp.$fastaName.name.txt
paste tmp.$fastaName.name.txt tmp.$fastaName.hostname.txt tmp.$fastaName.date.txt >> $pwd/interproscan_threads.log

rm tmp.$fastaName*

mv $fastaName* $pwd/interproscanOutput/

cd ../

rm -r $fastaName

};

} #end of preparing sh files

exit;

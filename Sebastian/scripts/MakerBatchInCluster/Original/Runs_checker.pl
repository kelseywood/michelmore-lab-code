#!/usr/bin/perl

######################################################
#                                                    #
#                Maker runs checker                  #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to check the exit status of all the maker runs and determine which contigs had a fail annotation
# Part of the main pipeline "Maker_qsub_annotator.pl"

use strict;
use warnings;

my $datastoreLog = $ARGV[0];

open(DATASTORELOG, $datastoreLog);

open(FAILEDRUNS, ">failed_runs.log");
open(FAILEDSEQ, ">failed_Sequences.log");

while( my $datastoreLogLine = <DATASTORELOG>) {

       	my @logLine = split("\t", $datastoreLogLine);

	if( !( ($datastoreLogLine =~ /.*FINISHED$/ ) || ($datastoreLogLine =~ /.*STARTED$/ ) ) && defined($logLine[1]) ) {

		print FAILEDSEQ $datastoreLogLine;

                my @failedFile = split("/", $logLine[1]);

                $failedFile[0] =~ s/_datastore//;

		print FAILEDRUNS $failedFile[0], "\n";

        } #end if is not finished or started

} #end while for datastorelog

exit;

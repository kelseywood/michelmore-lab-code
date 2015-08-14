#!/usr/bin/perl

######################################################
#                                                    #
#                     Find GFF                       #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to find the results gff's and datastore logs in the maker output (from a "-nodatastore" maker command)
# Part of the main pipeline "Maker_qsub_annotator.pl"

use strict;
use warnings;

my $makerRun = $ARGV[0];
my $outDir = $ARGV[1];
my $nodatastore = $ARGV[2];

system("mkdir $outDir") unless -d $outDir;

#open subdir and make new path
opendir(RUNSRESULTS, "$makerRun");

while(my $datastore = readdir(RUNSRESULTS) ) {

	#Find the datastore
	if($datastore =~ /datastore$/) {

		opendir(DATASTORE, "$makerRun/$datastore");

		#if is a directory add to dir list or add to gff list if .gff extension

		while(my $subdatastore = readdir(DATASTORE) ) {

			if( -d "$makerRun/$datastore/$subdatastore") {

				if($nodatastore eq "YES") {

					opendir(SEQUENCERESULT, "$makerRun/$datastore/$subdatastore");

					while( my $result = readdir(SEQUENCERESULT) ) {

						if( $result =~ /.gff$/) {

							system("cp $makerRun/$datastore/$subdatastore/$result ../GFF_files/$result");

						} #end if for gff files
					} #end while reading result dir

				} else {

					opendir(SUBDIRONE, "$makerRun/$datastore/$subdatastore");

					while( my $subdirtwo = readdir(SUBDIRONE) ) {

						opendir(SUBDIRTWO, "$makerRun/$datastore/$subdatastore/$subdirtwo");

						while( my $sequence = readdir(SUBDIRTWO) ) {

							opendir(SEQUENCERESULT, "$makerRun/$datastore/$subdatastore/$subdirtwo/$sequence");

							while( my $result = readdir(SEQUENCERESULT) ) {

								if( $result =~ /.gff$/) {

									system("cp $makerRun/$datastore/$subdatastore/$subdirtwo/$sequence/$result ../GFF_files/$result");

								} #end if for gff files

							} #end while reading result dir

						} #end while reading subdatastore 2 dir

					} #end while reading subdatastore 1 dir

				} #end if type nodatastore especify

			} #end if it is a sequence directory

		} #end reading datastore directory

	} elsif ($datastore =~ /master_datastore_index.log$/) {
		system("cp $makerRun/$datastore ../maker_datastore_logs/$makerRun.master_datastore_index.log");
	} #end if is datastore or log

} #end while maker run directory

exit;

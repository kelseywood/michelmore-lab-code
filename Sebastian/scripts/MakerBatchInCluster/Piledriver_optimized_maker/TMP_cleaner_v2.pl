#!/usr/bin/perl

######################################################
#                                                    #
#                    TMP cleaner                     #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to send ssh commands to all nodes use over the annotation process and delete the tmp folders created
# Part of the main pipeline "Maker_qsub_annotator.pl"

use strict;
use warnings;

my $threadsLog = $ARGV[0];
my $username = getpwuid( $< );

if(!defined($threadsLog)) { die "missing arguments: Please provide a maker_threads.log file from the main Maker_qsub_annotator pipeline\n"}


open(THREADSLOG, $threadsLog);

my %nodes;

while( my $threadLogLine = <THREADSLOG>) {

	if($threadLogLine =~ /.*started.*/) {

	      	my @logLine = split("\t", $threadLogLine);

		$nodes{$logLine[1]} = 1;

        } #end if it's a starting line

} #end while for threadslog

foreach my $node (keys %nodes) {

	system ("ssh $node 'rm -r /state/partition1/$username'");

	print "ssh $node 'rm -r /state/partition1/$username'\n";

} #end foreach $node

exit;

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

my $nodeList = $ARGV[0];
my $command = $ARGV[1];

if(!defined($nodeList)) { die "missing arguments: Please provide a maker_threads.log file from the main Maker_qsub_annotator pipeline\n"}


open(NODES, $nodeList) or die "Can't open $nodeList\n";

while( my $node = <NODES>) {

	chomp($node);

	system ("ssh $node '$command'");

	print "ssh $node '$command'\n";

} #end foreach $node

exit;

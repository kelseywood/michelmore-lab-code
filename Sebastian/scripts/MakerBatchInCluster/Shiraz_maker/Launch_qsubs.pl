#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use Getopt::Long;
use POSIX;

#########################################
#                                       #
#       Running maker from qsubs        #
#                                       #
#########################################

if (!defined($ARGV[0])) { print "Please provide a list of qsub files\n"; die "Missing arguments\n"; }

my $qsubsFile = $ARGV[0];

my $queue;
if (defined($ARGV[1]) ) {
	$queue = $ARGV[1];
}

open(QSUBS, $qsubsFile);

while(my $qsub = <QSUBS>) {

	chomp($qsub);

	open(MAKERTHREADS, ">>maker_threads.log");

	if(defined($queue) ) {

		print MAKERTHREADS $qsub, " submitted  on ", getTime(), " to queue ", $queue, "\n";

		system("qsub -q $queue $qsub");

	} else {

		print MAKERTHREADS $qsub, " submitted  on ", getTime(), " to default queue\n";

		system("qsub $qsub");

	} #end if specific queue requested

	close(MAKERTHREADS);

} #end


###########################################################################################
###########################################################################################
#
#   END OF SCRIPT
#
###########################################################################################
###########################################################################################


####################################
#                                  #
#          Get time sub            #
#                                  #
####################################

sub getTime {

	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	return "$hour-$minute-$second.$months[$month]-$dayOfMonth-$year";

} #end of getTime sub




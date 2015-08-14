#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;

#########################################
#                                       #
#       Running maker from qsubs        #
#                                       #
#########################################

#### Create output folders ####

system("mkdir maker_run_logs") unless -d "maker_run_logs";

system("mkdir qsubLogs") unless -d "qsubLogs";
system("mkdir qsubs") unless -d "qsubs";

system("mkdir GFF_files") unless -d "GFF_files";
system("mkdir maker_datastore_logs") unless -d "maker_datastore_logs";

system("mkdir maker_runs") unless -d 'maker_runs';


my $qsubsFile;
my $requiredTime;

GetOptions (
		'qsubsFile' => \$qsubsFile,
		'requiredTime=s' => \$requiredTime,
	);

if (!defined($qsubsFile) || !defined($requiredTime)) { print "Please provide a list of qsub files and the requiredtime\n"; die "Missing arguments\n"; }

open(QSUBS, $qsubsFile);

while(my $qsub = <QSUBS>) {

	chomp($qsub);

	open(MAKERTHREADS, ">>maker_threads.log");

	print MAKERTHREADS $qsub, " submitted  on ", getTime(), " to default queue\n";

	system("qsub $qsub -l h_rt=$requiredTime");

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




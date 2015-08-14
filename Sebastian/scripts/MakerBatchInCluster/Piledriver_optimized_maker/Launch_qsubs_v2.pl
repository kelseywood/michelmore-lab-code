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

GetOptions (
		'qsubsFile=s' => \$qsubsFile,
	);

if (!defined($qsubsFile)) {
	print "Required arguments are missing\n";
	print "Usage Launch_qsubs_v2.pl -qsubsFile=<qsub file list>\n";
	die "Missing arguments\n";
} #end if missing arguments

open(QSUBS, $qsubsFile) or die "Can't open $qsubsFile";

while(my $qsub = <QSUBS>) {

	chomp($qsub);

	if ($qsub =~ m/^#/) {next}

	open(MAKERTHREADS, ">>maker_threads.log");

	system("echo 'Launching $qsub' >> maker_launchJobs.log");

	if (system("qsub $qsub >> maker_launchJobs.log")  != 0 ) {

		print MAKERTHREADS $qsub, " couldn't be submitted to the queue  on ", getTime(), ", please check the maker_launchJobs.log for the error state\n";

	} else {

		print MAKERTHREADS $qsub, " submitted  on ", getTime(), "\n";

	} #end of launching

	close(MAKERTHREADS);

	system("echo '' >> maker_launchJobs.log");

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




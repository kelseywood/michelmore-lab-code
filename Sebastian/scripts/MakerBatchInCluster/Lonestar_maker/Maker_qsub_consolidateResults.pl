#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Std;
use Getopt::Long;
use POSIX;

######################################################
#                                                    #
#           Maker qsub consolidate results           #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to consolidate the results from 
# multiple maker runs of Maker_qsub_prepareRuns.pl

my $gffprefix;
my $time = getTime();

GetOptions (
		'gffprefix=s' => \$gffprefix,
	);

if(!defined($gffprefix)) {die "please input a gffprefix\n"}

open(RUNLOG, ">$time.makerRun.log");

#########################################
#                                       #
#        Verifying output stats         #
#                                       #
#########################################

opendir(DATASTORELOGS, "maker_datastore_logs") or die "Can't open directory maker_datastore_logs";

while(my $datastore = readdir(DATASTORELOGS) ) {

	if ($datastore =~ /.*master_datastore_index.*/) {

		system("cat maker_datastore_logs/$datastore >> master_datastore_index.log");

	} #end if

} #end while datastore files

if (-f "master_datastore_index.log") {

	print RUNLOG "Runs_checker_v2.pl master_datastore_index.log\n\n";

	if (system("Runs_checker_v2.pl master_datastore_index.log") != 0) { die "FATAL ERROR: Runs_checker.pl failed\n" }

	print STDERR "Validation of the exit status completed\n\n";

} #end if master_datastore_index.log


#########################################
#                                       #
#       Consolidate resulting gff       #
#                                       #
#########################################

print RUNLOG "Annotation_consolidator_v2.pl $gffprefix GFF_files\n\n";

if (system("Annotation_consolidator_v2.pl $gffprefix GFF_files") != 0 ) { die "FATAL ERROR: Annotation_consolidator.pl failed\n" }

print STDERR "End with the consolidation step\n\n";

print STDERR "MAKER ANNOTATION COMPLETED\n\n";

print STDERR "For details in the output files please use -h option for extended help\n\n";

#########################################
#                                       #
#      Modify and rename log files      #
#                                       #
#########################################

#system("TMP_cleaner.pl maker_threads.log");

system("mv maker_opts.ctl $time.maker_opts.ctl");
if (-f "master_datastore_index.log") { system("mv master_datastore_index.log $time.master_datastore_index.log") }
if (-f "failed_Sequences.log") { system("mv failed_Sequences.log $time.failed_Sequences.log") }
if (-f "failed_runs.log") {system("mv failed_runs.log $time.failed_runs.log") }

open(MAKERTHREADS, ">>maker_threads.log");
print MAKERTHREADS "\n\nEND OF RUN\n\n";
close(MAKERTHREADS);

exit;


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




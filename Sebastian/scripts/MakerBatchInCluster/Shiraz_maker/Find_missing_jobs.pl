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

if(!defined($ARGV[0])) { die "missing arguments\n"}

my $makerThreadsLog = $ARGV[0];

open(MAKERTHREADSLOG, $makerThreadsLog) or die "can't open $makerThreadsLog\n";
open(MISSINGJOBS, ">Missing_jobs.txt");

my %submittedJobs;

while( my $makerThreadsLogLine = <MAKERTHREADSLOG>) {

       	if($makerThreadsLogLine =~ m/submitted/) {

my @logLine = split(" ", $makerThreadsLogLine);

my $qsub = $logLine[0];
my $job = $logLine[0];

$job =~ s/.*\///;
$job =~ s/\.[^.]+$//;

$submittedJobs{$job} = $qsub;

       	} elsif ($makerThreadsLogLine =~ m/finished/) {

my @logLine = split(" ", $makerThreadsLogLine);

if(defined($submittedJobs{$logLine[0]})) { delete($submittedJobs{$logLine[0]}) }

       	} #end elsif

} #end while for datastorelog

foreach my $missedJob (keys %submittedJobs) { print MISSINGJOBS $submittedJobs{$missedJob},"\n" }



exit;

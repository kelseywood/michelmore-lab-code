#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) {
	print "Please provide tab delimited file with files to run VMatch on\n";
	print "Vmatch will be run with the following conditions\n";
	print "-dbcluster 80 10 -d -s -showdesc 0 -seedlength 10 -exdrop 3 -identity 80\n";
	print "<file>	<FolderName>	<Output prefix>\n";
	die "No input provided, please provided something\n";
} #end else

open(INPUT, $ARGV[0]) or die "Can't open $ARGV[0]\n";

open(LOG, ">>$ARGV[0].log") or die "Can't open $ARGV[0].log";

while (my $line = <INPUT>) {

	if ( ($line eq "") || ($line =~ m/^#/) ) {next}

	chomp($line);

	my @data = split("\t", $line);

	if(!(defined$data[2])) {next}

	my $file = $data[0];

	my $folder = $data[1];

	my $prefix = $data[2];

	if( !(-d $folder)) { system("mkdir $folder") }

	if( system("vmatch -dbcluster 80 10 $folder/$prefix -d -s -showdesc 0 -seedlength 10 -exdrop 3 -identity 80 $file > $folder/$prefix.All.txt 2> $folder/$prefix.log") ) {

		print STDERR "Couldn't run vmatch for $file\n";

	} else {

		print LOG "vmatch -dbcluster 80 10 $folder/$prefix -d -s -showdesc 0 -seedlength 10 -exdrop 3 -identity 80 $file > $folder/$prefix.All.txt 2> $folder/$prefix.log\n";

		print "Done with $file\n";

	} #end runnning vmatch

} #end of running to files

print "Done running vmatch\n\n";

close(INPUT);

exit;

	

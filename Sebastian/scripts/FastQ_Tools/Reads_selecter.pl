#!/usr/bin/perl
use strict;
use warnings;

###################################################
#
#
#	Reads selecter
#
#
#
###################################################

#
# This script was prepared to split fasta files into
# smaller files containing a specific number of
# sequences
#

if (@ARGV < 2) {

	print "Three are arguments are needed, please inpu them.\n";
	print "Fastq file where to retrieve the sequence\n";
        print "Text file with the reads to retrieve.\n";
        print "Outfile prefix.\n";

        exit 0;

} #end if

my %reads;

open(READS, $ARGV[1]);

while(my $read = <READS>) {

	chomp($read);

        $reads{$read} = "1";

} #end while

print scalar(keys %reads), " will be look for\n";

open(FASTQ, $ARGV[0]);
open(RESULTS, ">$ARGV[2].fastq");

while (my $header = <FASTQ>) {

	chomp($header);

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	$header =~ s/@//;

#	print $header, "\n";

	if(defined($reads{$header})) {

		print RESULTS "@", $header, "\n";
		print RESULTS $sequence;
		print RESULTS $header2;
		print RESULTS $qual;

	} #end if read to be retrieve

} #end while

close(FASTQ);
close(RESULTS);

exit;

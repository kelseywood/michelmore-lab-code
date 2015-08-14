#!/usr/bin/perl

###################################################
#
#
#	Fasta Sequence Splitter
#
#
#
###################################################

use strict;
use warnings;
use File::Basename;

#
# This script was prepared to split fasta files into
# smaller files containing a specific number of
# sequences
#

if (@ARGV < 2) {

	print "Two are arguments are needed, please inpu them.\n";
        print "Fasta file, number of sequences per file and output directory.\n";

        exit 0;

} #end if

my $fastaFile = $ARGV[0];
my $numSeqs = $ARGV[1];
my $dir = $ARGV[2];
my $outputFile = 1;
my $sequence = 1;

my $basename = basename($fastaFile);

system("mkdir $dir") unless -d $dir;

open(FASTA, $fastaFile) or die "Can't open file $fastaFile\n";
open(RESULTS, ">$dir/$basename.$outputFile.fasta");

while (my $line = <FASTA>) {

	if ((substr($line,0,1)) ne ">") {

	        print RESULTS $line;

        } else {

		if($sequence <= $numSeqs) {

                	print RESULTS $line;
                        ++$sequence;

                } else {

                	close(RESULTS);
                        ++$outputFile;
                        $sequence = 1;
                        open(RESULTS, ">$dir/$basename.$outputFile.fasta");
                	print RESULTS $line;
                        ++$sequence;

                } #end else

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;







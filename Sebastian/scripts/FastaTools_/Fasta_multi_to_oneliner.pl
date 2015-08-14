#!/usr/bin/perl

use strict; use warnings;

###################################################
#
#
#	Fasta multi to one liner
#
#
#
###################################################

#
# This script was prepared to convert multiple line fasta to one liner fastas
#

if (@ARGV < 0) {

	print "One argument is needed, please input them.\n";
        print "Fasta file to convert.\n";

        exit 0;

} #end if

open(FASTA, $ARGV[0]) or die "Can't open file $ARGV[0]\n";
open(RESULTS, ">$ARGV[0].oneLiner.fasta");

while (my $line = <FASTA>) {

	if ((substr($line,0,1)) ne ">") {

		chomp($line);
	        print RESULTS $line;

        } else {

		print RESULTS "\n";
		print RESULTS $line;

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;







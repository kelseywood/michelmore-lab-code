#!/usr/bin/perl

use strict;
use warnings;

if (!(defined($ARGV[1])) ) {

	print "Please input Fasta file where to filter and size threeshold\n";
	die "Missing arguments\n";

} #end if missing arguments

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_biggerthan_$ARGV[1].fasta");

my $header;
my $sequence = "";

while(my $line = <FASTAINPUT>) {

	if ($line =~ /^>/) {

		if( (length($sequence) > $ARGV[1]) && defined($header)) {

			print OUTPUT $header;
                	print OUTPUT $sequence, "\n";

                } #end else

	        $header = $line;
		$sequence = "";

        } else {

		chomp($line);
		$sequence .= $line;

        } #end else

} #end while

if( (length($sequence) > $ARGV[1]) && defined($header)) {

	print OUTPUT $header;
        print OUTPUT $sequence, "\n";

} #end else

close(FASTAINPUT);
close(OUTPUT);

exit;


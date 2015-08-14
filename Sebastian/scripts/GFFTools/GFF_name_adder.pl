#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Usage: print input GFF file from where add the name.\n";
	die;

} #end if


open(GFFINPUT, $ARGV[0]);
open(GFFOUTPUT, ">$ARGV[0].plusName.gff");

while (my $GFFline = <GFFINPUT>) {

	chomp($GFFline);

        my @data = split ("\t", $GFFline);

        my @attributes = split(";", $data[8]);

        if(substr($attributes[0], 0,2) eq "ID") {

	        my $ID = substr($attributes[0], 3,50);

        	print GFFOUTPUT $data[0], "\t", $data[1], "\t", $data[2], "\t", $data[3], "\t", $data[4], "\t", $data[5], "\t", $data[6], "\t", $data[7], "\tName=",$ID, ";", $data[8], "\n";

        } else {

                print GFFOUTPUT $GFFline, "\n";

        } #end else

} #end while

close(GFFINPUT);
close(GFFOUTPUT);

exit;











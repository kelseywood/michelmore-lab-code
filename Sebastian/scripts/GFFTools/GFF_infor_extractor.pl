#!/usr/bin/perl
use strict;
use warnings;

###################################################
#
#
#	GFF3 ID/Name/Parent Changer
#
#
#
###################################################

#

if (@ARGV < 2) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to extract data.\n";
	print "perl Fasta_renamer.pl <GFFFile> <type>\n";
        exit 0;

} #end if

my %loci = ();


open(GFF, $ARGV[0]);
open(COORDS, ">$ARGV[0].$ARGV[1].coordinates.txt");
open(INFO, ">$ARGV[0].$ARGV[1].information.txt");

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	if ($GFF_line[2] eq $ARGV[1]) {

		my @features = split (";", $GFF_line[8]);

		my @features_toPrint;

		my $name;

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "ID") {

			} elsif ($value[0] eq "Name") {
				$name = $value[1];

			} else {
				push(@features_toPrint, $feature);

			} #end else

		} #end foreach

		print COORDS $name, "\t", $GFF_line[1], "\t", $GFF_line[0], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[6], "\n";
		print INFO $name, "\t", $GFF_line[1], "\t", $GFF_line[5], "\t", join(";", @features_toPrint), "\n";

	} #end if

} #end while

close(GFF);
close(COORDS);
close(INFO);
exit;

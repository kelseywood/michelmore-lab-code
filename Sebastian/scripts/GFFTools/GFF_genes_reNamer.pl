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

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to change the ID (ID is the first of the features).\n";
        print "Translation Table with sequence names to change.\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <translationTable>\n";
        exit 0;

} #end if

my %loci = ();

open(TRANSLATIONTABLE, $ARGV[1]);

while(my $codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        my @data = split("\t", $codeTranslation);

        $loci{$data[0]} = $data[1];

} #end while

open(GFF, $ARGV[0]);
open(RESULTS, ">$ARGV[0].rename.gff3");

print RESULTS "##gff-version 3\n";

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	my $ID = "";
	my $name = "";
	my $alias = "";
	my $parent;

	my @otherFeatures;

	foreach my $feature (@features) {

		my @value = split("=", $feature);

	        if($value[0] eq "ID") {

			if(defined($loci{$value[1]})) { $ID = "ID=".$loci{$value[1]}.";" }
			else { $ID = "ID=".$value[1].";" }

	        } elsif ($value[0] eq "Name") {

			if(defined($loci{$value[1]})) { $name = "Name=".$loci{$value[1]}.";"; $alias = $value[1] }
			else { $name = "Name=".$value[1].";" }

	        } elsif ($value[0] eq "Parent") {

			$parent = $value[1]

	        } else {
			push(@otherFeatures, $feature);

	        } #end else

	} #end foreach

	if($name eq "") { $ID = ""}

	if(defined($parent)) {

		my @parents = split(",", $parent);

		foreach my $par (@parents) {

			if(!defined($loci{$par})) {print "Error: Parent not found for $name, parent $par\n"; next}

			print RESULTS $GFF_line[0], "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $ID, $name, "Parent=", $loci{$par}, ";" , join(";", @otherFeatures), "\n";
		} #end printing for each parent

	} else {

		print RESULTS $GFF_line[0], "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $ID, $name, join(";", @otherFeatures), "\n";

	} #end if parent

} #end while

close(GFF);
close(RESULTS);

exit;

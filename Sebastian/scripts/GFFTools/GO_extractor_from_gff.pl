#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Usage: print input GFF file from where to extract data\n";

	die;

} #end if

my $GFFfile = $ARGV[0];

open(GFF, $GFFfile);
open(RESULTS, ">$GFFfile.GOs.txt");

while(my $gffLine = <GFF>) {

	my @fields = split("\t", $gffLine);

	if($fields[2] eq "mRNA") {

		my @attributes = split(";", $fields[8]);

		my $name = "";

		foreach my $attribute (@attributes) {

			my @valuePair = split("=", $attribute);

			if($valuePair[0] eq "ID") {

				$name = $valuePair[1];

			} elsif ($valuePair[0] eq "GO_annotations") {

				my @GOs = split("\\|", $valuePair[1]);

				foreach my $GO (@GOs) {

					my @data = split(":", $GO);

					print RESULTS $data[1], "\t", $data[3], "\t", $data[2], "\t", $name, "\t", $fields[0], "\t", $fields[3], "\t", $fields[4], "\n";

				} #end foreach go

			} #end if go found

		} #end foreach attribute

	} #end if mrna

} #end while

close(GFF);
close(RESULTS);


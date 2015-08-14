#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Usage: two files need it, please input them.\n";
	print "GFF file with the original annotations.\n";
	print "GFF file with the new annotations.\n";
	die;

} #end if

my %originals;

open(GFF1, $ARGV[0]);

while(my $gffLine1 = <GFF1>) {

	chomp($gffLine1);

	my @fields = split("\t", $gffLine1);

	my @attributes = split(";", $fields[8]);

	foreach my $attribute (@attributes) {

		my @data = split("=", $attribute);

		if ($data[0] eq "ID" || $data[0] eq "Name") {

			$originals{$data[1]} = 1;

		} #end if

	} #end foreach

} #end while

close(GFF1);

open(GFF2, $ARGV[1]);

while(my $gffLine2 = <GFF2>) {

	chomp($gffLine2);

	my @fields = split("\t", $gffLine2);

	my @attributes = split(";", $fields[8]);

	my $toOutput = "NO";

	foreach my $attribute (@attributes) {

		my @data = split("=", $attribute);

		if ($data[0] eq "ID") {

			if(exists($originals{$data[1]})) { } else {print $gffLine2; print "\n"; $toOutput = "YES" }

		} elsif ($data[0] eq "Name") {

			if(exists($originals{$data[1]})) { } else {print $gffLine2; print "\n"; $toOutput = "YES" }

		} #end if

	} #end foreach

#	if($toOutput eq "YES") {

#		print $gffLine2;
#		print "\n";

#	} #end if

} #end while

close(GFF2);


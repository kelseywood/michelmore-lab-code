#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Usage: GFF file where to look for the data and output text file name.\n";
	print "GFF_AED_extractr.pl In.gff Out.txt\n";
	die "Missing Arguments";

} #end if

my $in = $ARGV[0];
my $out = $ARGV[1];

open(GFF, $in) or die "Can't open $in\n";
open(OUTPUT, ">$out") or die "Can't open $out\n";

while(my $gffLine = <GFF>) {

	chomp($gffLine);

	my $firstChar = substr($gffLine,0,1);

	my @fields = split("\t", $gffLine);

	my $toprint = "no";

	if($firstChar ne "#" && $firstChar ne "" && defined($fields[2])) {
	
		my @features = split (";", $fields[8]);

		my  $name = "";

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "Name") {

				$name = $value[1];

			} elsif($value[0] eq "_AED") {

				print OUTPUT $name, "\t", $value[1], "\n";

		        } #end if AED

		} #end foreach

	} #end if

} #end while

close(GFF);
close(OUTPUT);


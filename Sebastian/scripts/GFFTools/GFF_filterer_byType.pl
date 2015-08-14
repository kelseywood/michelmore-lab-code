#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Usage: print input GFF file from where to extract data, and keywords to look on type colum.\n";
	print "GFF_filterer_byType.pl file.gff type1 type2 type3\n";

	die;

} #end if

my %types;

foreach my $type (@ARGV) {

	$types{$type} = 1;

} #end foreach


my $GFFfile = shift(@ARGV);

open(GFF, $GFFfile);

while(my $gffLine = <GFF>) {

	my $firstChar = substr($gffLine,0,1);

	my @fields = split("\t", $gffLine);

	if($firstChar ne "#" && $firstChar ne "" && defined($fields[2])) {
	
		if (exists($types{$fields[2]})) {

			print $gffLine;

		} #end if

	} #end if

} #end while

close(GFF);


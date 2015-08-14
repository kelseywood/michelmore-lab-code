#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Usage: print input GFF file from where to extract data, and keywords to look on source colum.\n";
	print "GFF_filterer_bySource.pl file.gff source1 source2 source3\n";

	die;

} #end if

my %sources;

foreach my $source (@ARGV) {

	$sources{$source} = 1;

} #end foreach


my $GFFfile = shift(@ARGV);

open(GFF, $GFFfile);

while(my $gffLine = <GFF>) {

	my $firstChar = substr($gffLine,0,1);

	my @fields = split("\t", $gffLine);

	if($firstChar ne "#" && $firstChar ne "" && defined($fields[1])) {
	
		if (exists($sources{$fields[1]})) {

			print $gffLine;

		} #end if

	} #end if

} #end while

close(GFF);


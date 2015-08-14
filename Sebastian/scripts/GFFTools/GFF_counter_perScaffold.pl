#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Usage: print input GFF file from where were to count.\n";
	print "GFF_counter_perScaffold.pl file.gff\n";

	die;

} #end if

my %counts;
my %sources;

my $GFFfile = ($ARGV[0]);

open(GFF, $GFFfile);

while(my $gffLine = <GFF>) {

	my $firstChar = substr($gffLine,0,1);

	my @fields = split("\t", $gffLine);

	if($firstChar ne "#" && $firstChar ne "" && defined($fields[2])) {
	
		if ($fields[2] eq "mRNA") {

			++$counts{$fields[0]}{$fields[1]};
			$sources{$fields[1]} = 1;

		} #end if

	} #end if

} #end while

close(GFF);


foreach my $source (keys %sources) {

	print "\t", $source;

} #end foreach

print "\n";

foreach my $scaffold(keys %counts) {

	print $scaffold;

	foreach my $source (keys %sources) {

		if(exists($counts{$scaffold}{$source})) {

			print "\t", $counts{$scaffold}{$source};

		} else {

			print "\t0";

		} #end else

	} #end foreach

	print "\n";

} #end foreach

exit;

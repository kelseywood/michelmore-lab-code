#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV < 2) {

	print "Usage: print input GFF file for what to mask, gff input file, gff output file.\n";
	print "GFF_masker.pl <regions.gff> <IN.gff>\n";

	die;

} #end if

my %maskingRegions;

open(GFF, $ARGV[0]);

while(my $gffLine1 = <GFF>) {

	chomp($gffLine1);

	my @fields = split("\t", $gffLine1);

	$maskingRegions{$fields[0]}{$fields[3]} = $fields[4];

} #end while

close(GFF);

open(INGFF, $ARGV[1]) or die "can't open gff input\n";;

open(OUTGFF, ">$ARGV[1].matchedFeatures.gff") or die "can't open gff output\n";

while(my $gffLine2 = <INGFF>) {

	chomp($gffLine2);

        my @data = split("\t", $gffLine2);

        if(defined($maskingRegions{$data[0]})) {

		my $overlap = "FALSE";

		foreach my $region (sort {$a<=>$b} keys %{$maskingRegions{$data[0]}}) {

			if( (($region <= $data[3]) && ($maskingRegions{$data[0]}{$region} >= $data[3]))
                         || (($region <= $data[4]) && ($maskingRegions{$data[0]}{$region} >= $data[4]))
                         || (($region <= $data[3]) && ($maskingRegions{$data[0]}{$region} >= $data[4]))
                         || (($region >= $data[3]) && ($maskingRegions{$data[0]}{$region} <= $data[4])) ) {

				$overlap = "TRUE";

			} #end if

		} #end foreach masking regions

		if($overlap ne "FALSE") {

			print OUTGFF $gffLine2;
		        print OUTGFF "\n";

		} #end else masked

	} #end if


} #end while fasta file

close (INGFF);
close (OUTGFF);

exit;
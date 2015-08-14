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

open(OUTGFF, ">$ARGV[1].masked.gff") or die "can't open gff output\n";
open(REMOVEDGFF, ">$ARGV[1].removedFeatures.gff") or die "can't open gff output\n";

while(my $gffLine2 = <INGFF>) {

	chomp($gffLine2);

        my @data = split("\t", $gffLine2);

        if(defined($maskingRegions{$data[0]})) {

		my $overlap = "FALSE";

		foreach my $region (sort {$a<=>$b} keys %{$maskingRegions{$data[0]}}) {

	                my $start_start = "";
	                my $start_end = "";
	                my $end_start = "";
	                my $end_end = "";

	                #Start of the gene against start of the evidence
	                if($region > $data[3]) { $start_start = "before"}
	                elsif ($region < $data[3]) { $start_start = "after"}
	                else { $start_start = "match"} #end else for start-start comparison

	                #Start of the gene against end of the evidence
	                if($region > $data[4]) { $start_end = "before"}
	                elsif ($region < $data[4]) { $start_end = "after"}
	                else { $start_end = "match"} #end else for start-end comparison

	                #End of the gene against start of the evidence
	                if($maskingRegions{$data[0]}{$region} > $data[3]) { $end_start = "before"}
	                elsif ($maskingRegions{$data[0]}{$region} < $data[3]) { $end_start = "after"}
	                else { $end_start = "match"} #end else for end-start comparison

	                #End of the gene against end of the evidence
	                if($maskingRegions{$data[0]}{$region} > $data[4]) { $end_end = "before"}
	                elsif ($maskingRegions{$data[0]}{$region} < $data[4]) {$end_end = "after"}
	                else { $end_end = "match"} #end else for end-end comparison

			if(
#                        !( ($start_start eq "before") && ($start_end eq "before") && ($end_start eq "before") && ($end_end eq "before") ) &&
#                        !( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") )
			( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) || 
			( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) ||
			( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) ||
			( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) ||
			( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) ||
			( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) ||
			( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) ||
			( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) ||
			( ($start_start eq "before") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) ||
			( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "after") ) ||
			( ($start_start eq "match") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) ||
			( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "match") ) ||
			( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") )
                        ) {
				$overlap = "TRUE";
			} #end if

		} #end foreach masking regions

		if($overlap eq "FALSE") {

			print OUTGFF $gffLine2;
		        print OUTGFF "\n";

		} else {

			print REMOVEDGFF $gffLine2;
		        print REMOVEDGFF "\n";

		} #end else masked

	} else {

		print OUTGFF $gffLine2;
	        print OUTGFF "\n";

	} #end else


} #end while fasta file

close (INGFF);
close (OUTGFF);

exit;

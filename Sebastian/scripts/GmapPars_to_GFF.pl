#!/usr/bin/perl
use strict;
use warnings;

#######################################################
#                                                     #
#                  GMAP Pars to GFF                   #
#                                                     #
#              Sebastian Reyes-Chin-Wo                #
#               sreyesch@ucdavis.edu                  #
#                                                     #
#######################################################

if(!defined($ARGV[0])) {

	print "usage\n";

	die "missing arguments\n";

} #end if arguments

my $type = "Something";

open(GMAPPARS, $ARGV[0]);
open(GFF, ">$ARGV[0].gff");

while (my $gmapLine = <GMAPPARS>) {

	chomp($gmapLine);

	my @data = split("\t", $gmapLine);

	my $sequence = shift(@data);
	my $numPaths = shift(@data);
	my $strand = shift(@data);
	my $reference = shift(@data);
	$reference =~ s/ //g;
	my $coordinates = shift(@data);
	my $exons = shift(@data);
	my $coverage = shift(@data);
	my $identity = shift(@data);

	print $exons, "\n";

	$coordinates =~ s/,//g;

	my @coord = split("--", $coordinates);

	my $exonStructure = shift(@data);
	$exonStructure =~ s/( )+/ /g;

	my @exonSet = split(":", $exonStructure);

	for(my $i=1; $i <= $numPaths;++$i) {

		print GFF "$reference\tgmap\tmRNA\t$coord[0]\t$coord[1]\t$identity\t$strand\t.\tID=$sequence-$i;Name=$sequence-$i;\n";

		my $start = $coord[0];
		my $end = 0;

		#while(@data) {
		for(my $j=1; $j <= $exons; ++$j) {

			my $blank = shift(@exonSet);
			my $info = shift(@exonSet);


			my @exonData = split(" ", $info);

			my $id = $exonData[2];
			$id =~ s/%//g;
			my $coordExon = $exonData[1];
			$coordExon =~ s/\(|\)//g;
			my @coordExonValues = split("-", $coordExon);
			my $intron = 0;
			if(defined($exonData[4])) {$intron = $exonData[4]; $intron =~ s/\.//g;}

			print $intron,"\n";

			$end = $start+($coordExonValues[1]-$coordExonValues[0]);

			print GFF "$reference\tgmap\tCDS\t$start\t$end\t$id\t$strand\t.\tParent=$sequence-$i;\n";

			$start = $end + $intron;

		} #end foreach exon

	} #end foreach path

} #end reading file while



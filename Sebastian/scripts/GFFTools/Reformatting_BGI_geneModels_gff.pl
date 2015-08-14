#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Usage: print input GFF file which to reformat.\n";
	print "Reformatting_BGI_geneModels_gff.pl file.gff\n";

	die;

} #end if

open(GFF, $ARGV[0]) or die "Couldn't open file $ARGV[0]\n";
open(OUTPUTGFF, ">$ARGV[0].makerComplying.gff3");

my %geneModels;

while(my $gffLine = <GFF>) {

	chomp($gffLine);

	my @fields = split("\t", $gffLine);

	if($fields[1] eq "maker") { print OUTPUTGFF $gffLine, "\n"; next }

	my $ID = "";

	my @attributes = split(";", $fields[8]);

	#Look for parent or name tags on attribute fields
	foreach my $attribute(@attributes) {

		if ($attribute =~ /^Parent=(\S+)/) {

			my @data = split ("=", $attribute);

			$ID = $data[1];

		} elsif ( $attribute =~ /^Name=(\S+)/ ) {

			my @data = split ("=", $attribute);

			$ID = $data[1];

		} #end elsif parent

	} #end foreach
	
	push(@{$geneModels{$ID}{$fields[2]}}, [@fields]);

} #end while

close(GFF);


foreach my $geneModel (keys %geneModels) {

#	print $geneModel, "\n";

	# Recoordinate mRNA by UTR_5
	if( defined($geneModels{$geneModel}{"UTR_5"}) ) {
		if($geneModels{$geneModel}{"mRNA"}[0][6] eq "+") { $geneModels{$geneModel}{"mRNA"}[0][3] = $geneModels{$geneModel}{"UTR_5"}[0][3];
		} elsif($geneModels{$geneModel}{"mRNA"}[0][6] eq "-") {	$geneModels{$geneModel}{"mRNA"}[0][4] = $geneModels{$geneModel}{"UTR_5"}[0][4];
		} # if stranding
	} #end if five prime utr present

	# Recoordinate mRNA by UTR_3
	if(defined($geneModels{$geneModel}{"UTR_3"}) ) {
		if($geneModels{$geneModel}{"mRNA"}[0][6] eq "+") { $geneModels{$geneModel}{"mRNA"}[0][4] = $geneModels{$geneModel}{"UTR_3"}[0][4];
		} elsif($geneModels{$geneModel}{"mRNA"}[0][6] eq "-") { $geneModels{$geneModel}{"mRNA"}[0][3] = $geneModels{$geneModel}{"UTR_3"}[0][3];
		} # if stranding
	} #end if five prime utr present

	# Get the mRNA line
	my @geneLine = @{$geneModels{$geneModel}{"mRNA"}[0]};

	# Change the type and set score to 0
	$geneLine[2] = "gene";
	$geneLine[5] = 0;

	# Get attributes
	my $gene;
	my @geneAttributes = split(";", pop(@geneLine));

	print OUTPUTGFF join("\t", @geneLine), "\t";

	# Look for the ID attribute
	foreach my $attribute(@geneAttributes) {

		if ($attribute =~ /^ID=(\S+)/) {

			my @data = split ("=", $attribute);
			$data[1] =~ s/\.[0-9]*$//;
			$gene = $data[1];

			print OUTPUTGFF "Name=", $gene, ";ID=", $gene, ";";

		} #end if ID

	} #end foreach

	print OUTPUTGFF "\n";

	# Get and print the mRNA line
	print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"mRNA"}[0]}), "Parent=", $gene, ";\n";

	my @firstCDS = @{$geneModels{$geneModel}{"CDS"}[0]};

	# Changing the type
	$firstCDS[2] = "exon";

	# Getting the attributes
	my $parent;
	my @CDSAttributes = split(";", pop(@firstCDS) );

	# Finding the parent
	foreach my $attribute (@CDSAttributes) {
		if ($attribute =~ /^Parent=(\S+)/) {
			my @data = split ("=", $attribute);
			$parent = $data[1];
		} #end else attributes
	} #end foreach

	# For single cds features
	if(scalar(@{$geneModels{$geneModel}{"CDS"}}) == 1) {

		# Recoordinating CDS base on UTR_5
		if( defined($geneModels{$geneModel}{"UTR_5"}) ) {
#			print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"UTR_5"}[0]}), "\n";
			if($firstCDS[6] eq "+") { $firstCDS[3] = $geneModels{$geneModel}{"UTR_5"}[0][3];
			} elsif($firstCDS[6] eq "-") { $firstCDS[4] = $geneModels{$geneModel}{"UTR_5"}[0][4];
			} # if stranding
		} #end if five prime utr present

		# Recoordinating CDS base on UTR_3
		if(defined($geneModels{$geneModel}{"UTR_3"}) ) {
#			print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"UTR_3"}[0]}), "\n";
			if($firstCDS[6] eq "+") { $firstCDS[4] = $geneModels{$geneModel}{"UTR_3"}[0][4];
			} elsif($firstCDS[6] eq "-") { $firstCDS[3] = $geneModels{$geneModel}{"UTR_3"}[0][3];
			} # if stranding
		} #end if five prime utr present

		# Remove the frame
		$firstCDS[7] = ".";

		# Printing exon line
		if(defined($CDSAttributes[0]) ) { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $parent, ":exon:1;" , join(";", @CDSAttributes), ";\n";
		} else { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $parent, ":exon:1;\n";
		} #end else extra atributes

		# Print CDS
		my @curCDS = @{$geneModels{$geneModel}{"CDS"}[0]};
		my $curAttributes = pop(@curCDS);

		if(defined($curAttributes) ) { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":cds:1;", $curAttributes, "\n"
		} else { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":cds:1;\n"
		} #end else extra atributes

	} else {

		# For multi cds features

		# Recoordinating CDS base on UTR_5
		if( defined($geneModels{$geneModel}{"UTR_5"}) ) {
#			print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"UTR_5"}[0]}), "\n";
			if($firstCDS[6] eq "+") { $firstCDS[3] = $geneModels{$geneModel}{"UTR_5"}[0][3];
			} elsif($firstCDS[6] eq "-") { $firstCDS[4] = $geneModels{$geneModel}{"UTR_5"}[0][4];
			} # if stranding
		} #end if five prime utr present

		# Remove the frame
		$firstCDS[7] = ".";

		# Printing first exon line
		if(defined($CDSAttributes[0]) ) { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $parent, ":exon:1;" , join(";", @CDSAttributes), ";\n";
		} else { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $parent, ":exon:1;\n";
		} #end else extra atributes

		# Printing the rest of the exons except the last one
		my $exonCount = 2;
		for(my $i = 1; $i < scalar(@{$geneModels{$geneModel}{"CDS"}})-1; ++$i) {

			my @curCDS = @{$geneModels{$geneModel}{"CDS"}[$i]};
			my $curAttributes = pop(@curCDS);
			$curCDS[2] = "exon";
			$curCDS[7] = ".";

			if(defined($curAttributes) ) { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":exon:", $exonCount, ";", $curAttributes, "\n"
			} else { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":exon:", $exonCount, ";\n"
			} #end else extra atributes

			++$exonCount;

		} # end for

		if ($exonCount-1 > 0) {

			my @lastCDS = @{$geneModels{$geneModel}{"CDS"}[$exonCount-2]};

			if(defined($geneModels{$geneModel}{"UTR_3"}) ) {
#				print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"UTR_3"}[0]}), "\n";
				if($lastCDS[6] eq "+") { $lastCDS[4] = $geneModels{$geneModel}{"UTR_3"}[0][4];
				} elsif($lastCDS[6] eq "-") { $lastCDS[3] = $geneModels{$geneModel}{"UTR_3"}[0][3];
				} # if stranding
			} #end if five prime utr present

			$lastCDS[2] = "exon";
			$lastCDS[7] = ".";

			my $lastCDSAttributes = pop(@lastCDS);
	
			if(defined($CDSAttributes[0]) ) { print OUTPUTGFF join("\t", @lastCDS), "\tID=", $parent, ":exon:", $exonCount, ";", $lastCDSAttributes, "\n";
			} else { print OUTPUTGFF join("\t", @lastCDS), "\tID=", $parent, ":exon:", $exonCount, ";\n";
			} #end else extra atributes

		} #avoid duplicating exon in single exon gene models

		my $cdsCount = 1;

		# Print all the CDS's
		for(my $j = 0; $j < scalar(@{$geneModels{$geneModel}{"CDS"}}); ++$j) {
	
			my @curCDS = @{$geneModels{$geneModel}{"CDS"}[$j]};
			my $curAttributes = pop(@curCDS);

			if(defined($curAttributes) ) { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":cds:", $cdsCount, ";", $curAttributes, "\n"
			} else { print OUTPUTGFF join("\t", @curCDS), "\tID=", $parent, ":cds:", $cdsCount, ";\n"
			} #end else extra atributes

			++$cdsCount;

		} #end for printing CDSs

	} #end else 

} #end foreach genemodel

close(OUTPUTGFF);

exit;




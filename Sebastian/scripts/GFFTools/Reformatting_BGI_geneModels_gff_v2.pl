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
#	if($fields[1] eq "maker") { next }

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
	my $mRNA;
	my $gene;
	my @geneAttributes = split(";", pop(@geneLine));

	print OUTPUTGFF join("\t", @geneLine), "\t";

	# Look for the ID attribute
	foreach my $attribute(@geneAttributes) {

		if ($attribute =~ /^ID=(\S+)/) {

			my @data = split ("=", $attribute);
			$mRNA = $data[1];
			$gene = $data[1];
			$gene =~ s/\.[0-9]*$//;

			print OUTPUTGFF "ID=", $gene, ";Name=", $gene, ";";

		} #end if ID

	} #end foreach

	print OUTPUTGFF "\n";

	# Get and print the mRNA line
	print OUTPUTGFF join("\t", @{$geneModels{$geneModel}{"mRNA"}[0]}), "Parent=", $gene, ";\n";

	my @CDSs = @{$geneModels{$geneModel}{"CDS"}};
	my @exons = @{$geneModels{$geneModel}{"CDS"}};

	my @firstExon = @{shift(@exons)};
	my @firstCDS = @firstExon;

	# Changing the type
	$firstExon[2] = "exon";

	# Remove the frame
	$firstExon[7] = ".";

	# Getting the attributes
	my $ExonAttributes = pop(@firstExon);


	# Recoordinating CDS base on UTR_5
	if( defined($geneModels{$geneModel}{"UTR_5"}) ) {
		if($firstExon[6] eq "+") { $firstExon[3] = $geneModels{$geneModel}{"UTR_5"}[0][3];
		} elsif($firstExon[6] eq "-") { $firstExon[4] = $geneModels{$geneModel}{"UTR_5"}[0][4];
		} # if stranding
	} #end if five prime utr present


	# Split analysis between monoexonic gene models and gene models with multiple exons
	if(scalar(@{$geneModels{$geneModel}{"CDS"}}) == 1) {

		# For monoexonic verify the 3 primer UTR
		
		# Recoordinating CDS base on UTR_3
		if(defined($geneModels{$geneModel}{"UTR_3"}) ) {
			if($firstExon[6] eq "+") { $firstExon[4] = $geneModels{$geneModel}{"UTR_3"}[0][4];
			} elsif($firstExon[6] eq "-") { $firstExon[3] = $geneModels{$geneModel}{"UTR_3"}[0][3];
			} # if stranding
		} #end if five prime utr present

		# Printing exon line
		if(defined($ExonAttributes) ) { print OUTPUTGFF join("\t", @firstExon), "\tID=", $mRNA, ":exon:1;" , $ExonAttributes, ";\n";
		} else { print OUTPUTGFF join("\t", @firstExon), "\tID=", $mRNA, ":exon:1;\n";
		} #end else extra atributes

		# Print CDS
		my $firstCDSAttributes = pop(@firstCDS);

		if(defined($firstCDSAttributes) ) { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $mRNA, ":cds:1;", $firstCDSAttributes, "\n"
		} else { print OUTPUTGFF join("\t", @firstCDS), "\tID=", $mRNA, ":cds:1;\n";
		} #end else extra atributes

	} else {

		# For multi cds features

		# Recoordinating exon base on UTR_5
		if( defined($geneModels{$geneModel}{"UTR_5"}) ) {
			if($firstExon[6] eq "+") { $firstExon[3] = $geneModels{$geneModel}{"UTR_5"}[0][3];
			} elsif($firstExon[6] eq "-") { $firstExon[4] = $geneModels{$geneModel}{"UTR_5"}[0][4];
			} # if stranding
		} #end if five prime utr present

		# Printing first exon line
		if(defined($ExonAttributes) ) { print OUTPUTGFF join("\t", @firstExon), "\tID=", $mRNA, ":exon:1;" , $ExonAttributes, ";\n";
		} else { print OUTPUTGFF join("\t", @firstExon), "\tID=", $mRNA, ":exon:1;\n";
		} #end else extra atributes

		my @lastExon = @{pop(@exons)};

		# Printing the rest of the exons except the last one
		my $exonCount = 2;

		foreach my $exon (@exons) {

			my @curExon = @{$exon};
			my $curAttributes = pop(@curExon);
			$curExon[2] = "exon";
			$curExon[7] = ".";

			if(defined($curAttributes) ) { print OUTPUTGFF join("\t", @curExon), "\tID=", $mRNA, ":exon:", $exonCount, ";", $curAttributes, "\n"
			} else { print OUTPUTGFF join("\t", @curExon), "\tID=", $mRNA, ":exon:", $exonCount, ";\n"
			} #end else extra atributes

			++$exonCount;

		} # end for

		# Print the last Exon

		if(defined($geneModels{$geneModel}{"UTR_3"}) ) {
			if($lastExon[6] eq "+") { $lastExon[4] = $geneModels{$geneModel}{"UTR_3"}[0][4];
			} elsif($lastExon[6] eq "-") { $lastExon[3] = $geneModels{$geneModel}{"UTR_3"}[0][3];
			} # if stranding
		} #end if five prime utr present

		$lastExon[2] = "exon";
		$lastExon[7] = ".";

		my $lastExonAttributes = pop(@lastExon);
	
		if(defined($lastExonAttributes) ) { print OUTPUTGFF join("\t", @lastExon), "\tID=", $mRNA, ":exon:", $exonCount, ";", $lastExonAttributes, "\n";
		} else { print OUTPUTGFF join("\t", @lastExon), "\tID=", $mRNA, ":exon:", $exonCount, ";\n";
		} #end else extra atributes

		my $cdsCount = 1;

		# Print all the CDS's
		foreach my $curCDS (@{$geneModels{$geneModel}{"CDS"}}) {
	
			my $curAttributes = pop(@{$curCDS});

			if(defined($curAttributes) ) { print OUTPUTGFF join("\t", @{$curCDS}), "\tID=", $mRNA, ":cds:", $cdsCount, ";", $curAttributes, "\n"
			} else { print OUTPUTGFF join("\t", @{$curCDS}), "\tID=", $mRNA, ":cds:", $cdsCount, ";\n"
			} #end else extra atributes

			++$cdsCount;

		} #end for printing CDSs

	} #end else

} #end foreach genemodel

close(OUTPUTGFF);

exit;

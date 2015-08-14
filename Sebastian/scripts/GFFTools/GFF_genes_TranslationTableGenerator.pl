#!/usr/bin/perl

use strict;
use warnings;

###################################################
#
#
#	GFF3 ID/Name/Parent Changer
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to change the ID (ID is the first of the features).\n";
        print "Translation Table with sequence names to change.\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <Gene name prefix>\n";
        exit 0;

} #end if

my $prefix = $ARGV[1];

open(GFF, $ARGV[0]);

my %genes;

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	if($GFF_line[2] eq "mRNA") {

		my $ID = "";
		my $name = "";
		my $parent;

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "ID") {

				$ID = $value[1]

			} elsif ($value[0] eq "Parent") {

				$parent = $value[1]

			} #end else

		} #end foreach

		push(@{$genes{$parent}}, $ID);

	} #end if it's an mRNA

} #end while

close(GFF);

open(RESULTS, ">$ARGV[0].TranslationTable.txt");
print RESULTS "#OriginalName\tNewName\n";

my $geneID = 1;

foreach my $gene (keys %genes) {

	print RESULTS $gene, "\t", $prefix, sprintf("%06d", $geneID), "\n";

	my $mRNAcount = 1;

	foreach my $mRNA (@{$genes{$gene}}) {

		print RESULTS $mRNA, "\t", $prefix, sprintf("%06d", $geneID), ".", $mRNAcount, "\n";
		++$mRNAcount;

	} #end printing mRNAs

	++$geneID;

} #end printing new ID's

close(RESULTS);

exit










close(RESULTS);

exit;

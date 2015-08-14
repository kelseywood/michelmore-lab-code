#!/usr/bin/perl
use strict; use warnings;
###################################################
#
#
#	GFF3 Orphan feature remover
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "One argument is needed, please input them.\n";
	print "GFF output file from GMAP (using '-f 2' option).\n";
	print "perl GMAP_geneGFF3_to_makerLikeGFF.pl <GFFFile>\n";
        exit 0;

} #end if

my %ids;

open(GFF, $ARGV[0]);

my $gff = $ARGV[0];

$gff =~ s{\.[^.]+$}{};

open(MAKERGFF, ">$gff.makerCompliant.gff");

while (my $line = <GFF>) {

	if ($line =~ /^#/) {next}

	chomp($line);

        my @GFF_line = split ("\t", $line);

	if ( $GFF_line[2] eq "mRNA" ) {

#		my $ID;

#		my @features = split (";", $GFF_line[8]);
#		foreach my $feature (@features) {
#			my @value = split("=", $feature);
#			if($value[0] eq "ID") { $ID = $value[1]}
#		} #end foreach

		$GFF_line[2] = "expressed_sequence_match";

		print MAKERGFF join ("\t", @GFF_line), "\n";

	} elsif ( ($GFF_line[2] eq "exon") ) {

#		my $parent;

#		my @features = split (";", $GFF_line[8]);
#		foreach my $feature (@features) {
#			my @value = split("=", $feature);
#			if($value[0] eq "Parent") { $parent = $value[1]}
#		} #end foreach

		$GFF_line[2] = "match_part";

		print MAKERGFF join ("\t", @GFF_line), "\n";

	} #end if types

} #end while

close(GFF);
close(MAKERGFF);

exit;

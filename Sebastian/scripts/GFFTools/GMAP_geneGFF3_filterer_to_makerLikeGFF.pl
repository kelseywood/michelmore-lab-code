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

if (@ARGV < 3) {

	print "Three arguments are needed, please input them.\n";
	print "GFF output file from GMAP (using '-f 2' option).\n";
	print "Minimun coverage require.\n";
	print "Minimun identity require.\n";
	print "perl GMAP_geneGFF3_to_makerLikeGFF.pl <GFFFile> <minimun coverage> <minimun identity>\n";
        exit 0;

} #end if

my %ids;

open(GFF, $ARGV[0]);

my $gff = $ARGV[0];

my $minCoverage = $ARGV[1];

my $minIdentity = $ARGV[2];

$gff =~ s{\.[^.]+$}{};

open(GOOD, ">$gff.id$minIdentity.cov$minCoverage.makerCompliant.gff");

my %goodAlignmentsParents;
my %goodAlignmentsIDs;

while (my $line = <GFF>) {

	if ($line =~ /^#/) {next}

	chomp($line);

        my @GFF_line = split ("\t", $line);

	if ($GFF_line[2] eq "mRNA") {

		my @features = split (";", $GFF_line[8]);
		
		my $ID;

		my $parent;

		my $coverage;

		my $identity;

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "ID") { $ID = $value[1];

			} elsif ($value[0] eq "Parent") { $parent = $value[1];

			} elsif ($value[0] eq "coverage") { $coverage = $value[1];

			} elsif ($value[0] eq "identity") {$identity = $value[1];

			} #end elsif

		} #end foreach

		if ( ($identity >= $minIdentity) && ($coverage >= $minCoverage) ) {

			$goodAlignmentsIDs{$ID} = 1;
			$goodAlignmentsParents{$parent} = 1;

		} #end if it's a good alignment

	} #end if it's mRNA

} #end while looking mRNA

close(GFF);

print "Done checking alignments\n";
print "Will print the file now\n\n";

open(GFF, $ARGV[0]);

while (my $line = <GFF>) {

	if ($line =~ /^#/) {next}

	chomp($line);

        my @GFF_line = split ("\t", $line);

	if ( $GFF_line[2] eq "mRNA" ) {

		my $ID;

		my @features = split (";", $GFF_line[8]);
		foreach my $feature (@features) {
			my @value = split("=", $feature);
			if($value[0] eq "ID") { $ID = $value[1]}
		} #end foreach

		$GFF_line[2] = "expressed_sequence_match";

		if (defined($goodAlignmentsIDs{$ID}) ) { print GOOD join ("\t", @GFF_line), "\n"}

	} elsif ( ($GFF_line[2] eq "exon") ) {

		my $parent;

		my @features = split (";", $GFF_line[8]);
		foreach my $feature (@features) {
			my @value = split("=", $feature);
			if($value[0] eq "Parent") { $parent = $value[1]}
		} #end foreach

		$GFF_line[2] = "match_part";

		if (defined($goodAlignmentsIDs{$parent}) ) { print GOOD join ("\t", @GFF_line), "\n"}

	} #end if types

} #end while

close(GFF);
close(GOOD);

exit;

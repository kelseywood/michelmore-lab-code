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
	print "GMAP_GFF3_filterer.pl <GFFFile> <Min coverage (int)> <Min identity (int)>\n";
        exit 0;

} #end if

my %ids;

open(GFF, $ARGV[0]) or die "Can't open $ARGV[0]\n";

my $gff = $ARGV[0];

my $minCoverage = $ARGV[1];

my $minIdentity = $ARGV[2];

$gff =~ s{\.[^.]+$}{};

open(GOOD, ">$gff.id$minIdentity.cov$minCoverage.good.gff");
open(BAD, ">$gff.id$minIdentity.cov$minCoverage.bad.gff");
open(STATS, ">$gff.id$minIdentity.cov$minCoverage.stats.txt");

my %goodAlignmentsParents;
my %goodAlignmentsIDs;

my $totalAlignments;

while (my $line = <GFF>) {

	if ($line =~ /^#/) {next}

	chomp($line);

        my @GFF_line = split ("\t", $line);

	if ($GFF_line[2] eq "mRNA") {

		++$totalAlignments;

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

print STATS $gff, "\t", $minIdentity, "\t", $minCoverage, "\t", $totalAlignments, "\t", scalar(keys %goodAlignmentsIDs), "\n";

open(GFF, $ARGV[0]);

while (my $line = <GFF>) {

	if ($line =~ /^#/) {next}

	chomp($line);

        my @GFF_line = split ("\t", $line);

	if( $GFF_line[2] eq "gene") {

		my $ID;

		my @features = split (";", $GFF_line[8]);
		foreach my $feature (@features) {
			my @value = split("=", $feature);
			if($value[0] eq "ID") { $ID = $value[1]}
		} #end foreach

		if (defined($goodAlignmentsParents{$ID}) ) { print GOOD $line, "\n"} else { print BAD $line, "\n"}

	} elsif ( $GFF_line[2] eq "mRNA" ) {

		my $ID;

		my @features = split (";", $GFF_line[8]);
		foreach my $feature (@features) {
			my @value = split("=", $feature);
			if($value[0] eq "ID") { $ID = $value[1]}
		} #end foreach

		if (defined($goodAlignmentsIDs{$ID}) ) { print GOOD $line, "\n"} else { print BAD $line, "\n"}

	} elsif ( ($GFF_line[2] eq "CDS") || ($GFF_line[2] eq "exon") ) {

		my $parent;

		my @features = split (";", $GFF_line[8]);
		foreach my $feature (@features) {
			my @value = split("=", $feature);
			if($value[0] eq "Parent") { $parent = $value[1]}
		} #end foreach

		if (defined($goodAlignmentsIDs{$parent}) ) { print GOOD $line, "\n"} else { print BAD $line, "\n"}

	} #end if types

} #end while

close(GFF);
close(GOOD);
close(BAD);

exit;

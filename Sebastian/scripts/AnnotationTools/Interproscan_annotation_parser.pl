#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[1])) {

	print "Please provide an Interproscan tsv file and an output prefix.\n";
	print "Usage:\n";
	print "Interproscan_annotation_parser.pl <Interproscan tsv> <prefix>\n";
	die "Missing arguments\n\n";

} #end if

my $interproscan = $ARGV[0];
my $prefix = $ARGV[1];

my %annotation;

open(INTERPROSCAN, $interproscan);

while (my $line = <INTERPROSCAN> ) {

	chomp($line);

	my @annotationInfo = split("\t", $line);

	my @proteins = split("\\|", $annotationInfo[0]);

	shift(@annotationInfo);
	shift(@annotationInfo);
	shift(@annotationInfo);

	foreach my $protein (@proteins) {

		push(@{$annotation{$protein}{$annotationInfo[0]}}, join("\t", @annotationInfo) );

	} #end gathering general information

	if($line =~ m/GO:/) {

		my @GOs = split("\\|", $annotationInfo[10]);

		foreach my $GO (@GOs) {	foreach my $protein (@proteins) { $annotation{$protein}{"GOs"}{$GO} = 1 } }

	} # end if have GO annotation

} #end reading annotation info

my @annotationSources = (
	"Coils",
	"Gene3D",
	"Hamap",
	"PANTHER",
	"Pfam",
	"PIRSF",
	"PRINTS",
	"ProSitePatterns",
	"ProSiteProfiles",
	"SMART",
	"SUPERFAMILY",
	"TIGRFAM"
	);

foreach my $database (@annotationSources) {
	open(FILE, ">$prefix.$database.txt");
	close(FILE);
} #end opening files


open(FULLLIST, ">$prefix.completeDataset.txt");
open(GO, ">$prefix.GOs.txt");

foreach my $protein ( keys %annotation ) {

	foreach my $database (@annotationSources) {

		if(defined($annotation{$protein}{$database}) ) {

			open(OUT, ">>$prefix.$database.txt");

			foreach my $annotation (@{$annotation{$protein}{$database}}) {
				print FULLLIST $protein, "\t", $annotation, "\n";
				print OUT $protein, "\t", $annotation, "\n";
			} #end printing

			close(OUT);

		} #end if proteins has database hits

	} #end foreach database

	if(defined($annotation{$protein}{"GOs"}) ) {

		foreach my $go (keys %{$annotation{$protein}{"GOs"}} ) { print GO $protein, "\t", $go, "\n" }

	} #end if has GO's

} # end foreach protein


exit;

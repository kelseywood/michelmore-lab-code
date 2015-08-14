#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[2])) {

	print "Please provide an MCL output file, Interproscan tsv file and an output prefix.\n";
	print "Usage:\n";
	print "OrthoMCL_groups-Interproscan_annotation_merger.pl <MCL groups> <Interproscan tsv> <prefix>\n";
	die "Missing arguments\n\n";

} #end if

my $MCL = $ARGV[0];
my $interproscan = $ARGV[1];
my $prefix = $ARGV[2];


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
open(FAMILIES, ">$prefix.geneFamilies.txt");
open(GO, ">$prefix.GOs.txt");

open(GROUPS, $MCL);


while(my $line = <GROUPS>) {

	chomp($line);

	my @group = split(" ", $line);

	my $groupName = shift(@group);

	$groupName =~ s/://;

	foreach my $member (@group) {

		my @memberInfo = split("\\|", $member);

		print FAMILIES $groupName, "\t", $memberInfo[0], "\t", scalar(@group), "\t", $memberInfo[1], "\n";

		foreach my $database (@annotationSources) {

			if(defined($annotation{$memberInfo[1]}{$database}) ) {

				open(OUT, ">>$prefix.$database.txt");

				foreach my $annotation (@{$annotation{$memberInfo[1]}{$database}}) {
					print FULLLIST $groupName, "\t", $memberInfo[0], "\t", $memberInfo[1], "\t", $annotation, "\n";
					print OUT $groupName, "\t", $memberInfo[0], "\t", $memberInfo[1], "\t", $annotation, "\n";
				} #end printing

			} #end if

		} #end opening files

		if(defined($annotation{$memberInfo[1]}{"GOs"}) ) {

			foreach my $go (keys %{$annotation{$memberInfo[1]}{"GOs"}} ) { print GO $groupName, "\t", $memberInfo[0], "\t", $memberInfo[1], "\t", $go, "\n" }

		} #end if has GO's

	} #end foreach member

} #end while file

close(GROUPS);
















exit;

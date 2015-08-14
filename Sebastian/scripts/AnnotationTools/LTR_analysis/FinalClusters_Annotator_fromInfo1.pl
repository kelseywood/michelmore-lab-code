#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[1])) {
	print "Missing inputs, please provide a FinalClusteredSequence file and an info1 file\n";
	die;
} #end else

my %classifications;

open(INFO1, $ARGV[1]);

while( my $info1Line = <INFO1> ) {

	chomp($info1Line);

	my @info1 = split("\t", $info1Line);

	$classifications{$info1[0]} = $info1[2];

} #end reading info1

close(INFO1);

open(CLUSTERS, $ARGV[0]);

open(ANNOTATED, ">$ARGV[0].annotated.txt");

my %clusters;

my $header = <CLUSTERS>;

while (my $sequence = <CLUSTERS>) {

	chomp($sequence);

	my @data = split("\t", $sequence);

	print 	ANNOTATED $data[0], "\t", $data[1], "\t", $classifications{$data[1]}, "\n";

	++$clusters{$data[0]}{$classifications{$data[1]}};

} #end printing annotations

close(CLUSTERS);
close(ANNOTATED);

open(SUMMARY, ">$ARGV[0].familyCounts.txt");

foreach my $cluster (sort keys %clusters) {
	foreach my $family (sort { $clusters{$cluster}{$b} <=> $clusters{$cluster}{$a} } keys %{$clusters{$cluster}}) { print SUMMARY $cluster, "\t", $family, "\t", $clusters{$cluster}{$family}, "\n" }
} #end printing summary

close(SUMMARY);

exit;

#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) {
	print "Missing inputs, please provide a file with annotated clusters\n";
	print "Cluster\tAnnotation\tMembers\n";
	die;
} #end else

my %clusters;

open(CLUSTERS, $ARGV[0]);

while( my $clusterLine = <CLUSTERS> ) {

	chomp($clusterLine);

	my @data = split("\t", $clusterLine);

#	print STDERR $clusterLine, "\n";

	$clusters{$data[0]}{$data[1]} += $data[2];

} #end reading info1

close(CLUSTERS);

open(CLASSIFICATIONS, ">$ARGV[0].classified.txt");

foreach my $cluster (keys %clusters) {

	my @annotations = (sort { $clusters{$cluster}{$b} <=> $clusters{$cluster}{$a} } keys %{$clusters{$cluster}} );

	if( scalar(@annotations) == 1) {

		print CLASSIFICATIONS $cluster,"\t", $annotations[0], "\t", $clusters{$cluster}{$annotations[0]}, "\n";

	} else {

		if($clusters{$cluster}{$annotations[0]} > $clusters{$cluster}{$annotations[1]}) {

			print CLASSIFICATIONS $cluster,"\t", $annotations[0], "\t", $clusters{$cluster}{$annotations[0]}, "\n";

		} else {

			my %familyCounters;

			foreach my $annotation ( @annotations ) {

				my @family = split("-", $annotation);

				$familyCounters{$family[0]} += $clusters{$cluster}{$annotation}			

			} #end foreach annotation

			my @family = (sort { $familyCounters{$b} <=> $familyCounters{$a} } keys %familyCounters );

			if( scalar(@family) == 1) { print CLASSIFICATIONS $cluster,"\t", $family[0], "-X", "\t", $familyCounters{$family[0]}, "\n" }

			elsif($familyCounters{$family[0]} > $familyCounters{$family[1]}) { print CLASSIFICATIONS $cluster,"\t", $family[0], "-X", "\t", $familyCounters{$family[0]}, "\n" }
			

		} #end if it's actually bigger

	} #end if only have one annotation

} #end foreach cluster

close(CLASSIFICATIONS);


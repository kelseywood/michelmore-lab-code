#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[1])) {
	print "Missing inputs, please provide the output from VMatch_clusters_parser.pl and a output prefix\n";
	die
} #end else

my $input = $ARGV[0];

my $prefix = $ARGV[1];


open(INPUT, $input);
my $header = <INPUT>;
chomp($header);
my @Domains = split("\t", $header);
shift(@Domains);
my $numDomains = scalar(@Domains);


my %sequences;

while(my $sequence = <INPUT>) {

	chomp($sequence);

	my @states = split("\t", $sequence);

	my $seq = shift(@states);

	for(my $i = 0; $i < $numDomains; ++$i) { $sequences{$seq}{$Domains[$i]} = $states[$i] }

	$sequences{$seq}{'sequences'} = $sequence;

} #reading in data

my $numSequences = scalar(keys %sequences);

print "Finish loading data\n";
print $numSequences, " predicted LTR's were loaded\n";

my %clusters;
my $clusterID = 1;

foreach my $seq (keys %sequences) {

	if(!defined($sequences{$seq})) {next}

	#Creating the new cluster
	for(my $i = 0; $i < $numDomains; ++$i) { $clusters{$clusterID}{$Domains[$i]} = $sequences{$seq}{$Domains[$i]}} #end creating new cluster
	push(@{$clusters{$clusterID}{'sequences'}}, $sequences{$seq}{'sequences'});

	delete($sequences{$seq});

	my $completed = "no";

	while ($completed eq "no") {

		$completed = "yes";

		# Going to check all the sequences against the cluster profile
		foreach my $seq2 (keys %sequences) {

			my $isCompatible = "yes";
			my $matchCount = 0;

#			print $seq,"\t",$clusterID, "\t";

			for(my $i = 0; $i < $numDomains; ++$i) {

#				print $i,":", $clusters{$clusterID}{$Domains[$i]}, "|", $sequences{$seq2}{$Domains[$i]}, " ";

				if( ($clusters{$clusterID}{$Domains[$i]} eq "NA") || ($sequences{$seq2}{$Domains[$i]} eq "NA") ) {next}
				if( $clusters{$clusterID}{$Domains[$i]} != $sequences{$seq2}{$Domains[$i]} ) {$isCompatible = "no"}
				if ( $clusters{$clusterID}{$Domains[$i]} == $sequences{$seq2}{$Domains[$i]} ) {++$matchCount}

			} #end checking compatibility

#			print "\t", $isCompatible, "\t", $matchCount, "\n";

			if( ($isCompatible eq "yes") && ($matchCount != 0) ) {

				for(my $j = 0; $j < $numDomains; ++$j) { if($clusters{$clusterID}{$Domains[$j]} eq "NA") { $clusters{$clusterID}{$Domains[$j]} = $sequences{$seq2}{$Domains[$j]} } }

				push(@{$clusters{$clusterID}{'sequences'}}, $sequences{$seq2}{'sequences'});

				delete($sequences{$seq2});

				$completed = "no";

			} #end if it's compatible

		} #end checking clusters;

	} #end while looping around

	++$clusterID;

} #end while reading file

close(INPUT);

print "Completed clustering LTR's\n";
print scalar (keys %clusters), " initial cluster were found\n";

open(PROFILES, ">$prefix.ClusterProfiles.txt");
print PROFILES "ClusterID\t", join("\t", @Domains), "\n";

open(CLUSTERS, ">$prefix.ClusteredSequences.txt");
print CLUSTERS "ClusterID\tSequences\t", join("\t", @Domains), "\n";

open(STATS, ">$prefix.ClusterStats.txt");

foreach my $cluster ( sort {$a <=> $b} keys %clusters) {

	my $numValues = 0;

	print PROFILES "Cluster_", $cluster;
	for(my $i = 0; $i < $numDomains; ++$i) {
		if($clusters{$cluster}{$Domains[$i]} ne "NA") {++$numValues};
		print PROFILES "\t", $clusters{$cluster}{$Domains[$i]};
	} #end printing domain clusters
	print PROFILES "\n";

	foreach my $seq (@{$clusters{$cluster}{'sequences'}}) {	print CLUSTERS "Cluster_", $cluster, "\t", $seq, "\n" }

	print STATS "Cluster_", $cluster, "\t", scalar(@{$clusters{$cluster}{'sequences'}}), "\t", $numValues, "\n";

} #end printing cluster profiles

close(PROFILES);

exit;
















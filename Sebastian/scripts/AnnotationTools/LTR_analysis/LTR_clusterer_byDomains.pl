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

print "\nLoading data\n";

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
print $numSequences, " predicted LTR's were loaded\n\n";


#####################################

######### BUILDING CLUSTERS #########

#####################################

print "Building initial clusters\n";

my %clusters;
my $clusterID = 1;

foreach my $seq (keys %sequences) {

	if(!defined($sequences{$seq})) {next}

	#Creating the new cluster
	for(my $i = 0; $i < $numDomains; ++$i) { $clusters{$clusterID}{$Domains[$i]} = $sequences{$seq}{$Domains[$i]}} #end creating new cluster
	$clusters{$clusterID}{'sequences'}{$seq} = $sequences{$seq}{'sequences'};

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

				$clusters{$clusterID}{'sequences'}{$seq2} = $sequences{$seq2}{'sequences'};

				delete($sequences{$seq2});

				$completed = "no";

			} #end if it's compatible

		} #end checking clusters;

	} #end while looping around

	++$clusterID;

} #end while reading file

close(INPUT);

print "Completed clustering LTR's\n";
print scalar (keys %clusters), " initial cluster were found\n\n";


#####################################

##### PRINTING INITIAL CLUSTERS #####

#####################################

print "Printing initial clusters\n";

open(PROFILES, ">$prefix.InitialClusterProfiles.txt");
print PROFILES "ClusterID\t", join("\t", @Domains), "\n";

open(CLUSTERS, ">$prefix.InitialClusteredSequences.txt");
print CLUSTERS "ClusterID\tSequences\t", join("\t", @Domains), "\n";

open(STATS, ">$prefix.InitialClusterStats.txt");

foreach my $cluster ( sort {$a <=> $b} keys %clusters) {

	if(!defined($clusters{$cluster}{'sequences'})) {next}

	my $numValues = 0;

	print PROFILES "Cluster_", $cluster;
	for(my $i = 0; $i < $numDomains; ++$i) {
		if($clusters{$cluster}{$Domains[$i]} ne "NA") {++$numValues};
		print PROFILES "\t", $clusters{$cluster}{$Domains[$i]};
	} #end printing domain clusters
	print PROFILES "\n";

	foreach my $seq (keys %{$clusters{$cluster}{'sequences'}} ) {	print CLUSTERS "Cluster_", $cluster, "\t", $clusters{$cluster}{'sequences'}{$seq}, "\n" }

	print STATS "Cluster_", $cluster, "\t", scalar(keys %{$clusters{$cluster}{'sequences'}}), "\t", $numValues, "\n";

} #end printing cluster profiles

close(PROFILES);

print "Done printing initial clusters\n\n";

#####################################

####### FINDING AMBIGOUS LTRs #######

#####################################

my %ambigous;

print "Looking for ambigous LTR's\n";

foreach my $cluster ( sort {$a <=> $b} keys %clusters) {

	#Start looping through sequences
	foreach my $seq (keys %{$clusters{$cluster}{'sequences'}} ) {

		#Get states
		my @states = split("\t", $clusters{$cluster}{'sequences'}{$seq});
		my $seqName = shift(@states);

		#Variable to count how many clusters it's the LTR compatible
		my @compatibleClusters;

		#Loop through all the clusters
		foreach my $searchCluster ( sort {$a <=> $b} keys %clusters) {

			my $isCompatible = "yes";
			my $matchCount = 0;

			#Verify compatibility per domain
			for(my $i = 0; $i < $numDomains; ++$i) {

				if( ($clusters{$searchCluster}{$Domains[$i]} eq "NA") || ($states[$i] eq "NA") ) {next}
				if( $clusters{$searchCluster}{$Domains[$i]} != $states[$i] ) {$isCompatible = "no"}
				if ( $clusters{$searchCluster}{$Domains[$i]} == $states[$i] ) {++$matchCount}

			} #end checking compatibility

			#Verify compatibility with the cluster
			if( ($isCompatible eq "yes") && ($matchCount != 0) ) {push(@compatibleClusters,"Cluster_".$searchCluster)}

		} #end foreach cluster

		#Verify it's ambigous
		if(scalar(@compatibleClusters) > 1) {

			$ambigous{$seqName}{'sequence'} = $clusters{$cluster}{'sequences'}{$seq};
			$ambigous{$seqName}{'clusters'} =  join("|", @compatibleClusters);

			delete($clusters{$cluster}{'sequences'}{$seq});

		} #end if it's ambigous
	
	} #finish verifying sequence

} #end printing cluster profiles

print "Done looking for ambigous LTR's\n";
print scalar (keys %ambigous), " ambigous LTR's were found\n\n";


#####################################

###### PRINTING AMBIGOUS LTRs #######

#####################################

print "Time to print those ambigous LTR's\n\n";

open(AMBIGOUSPROFILE, ">$prefix.AmbigousProfiles.txt");
print AMBIGOUSPROFILE "SeqID\t", join("\t", @Domains), "\n";

open(AMBIGOUSSTATS, ">$prefix.AmbigousStats.txt");

foreach my $seq (keys %ambigous) {

	print AMBIGOUSPROFILE $ambigous{$seq}{'sequence'}, "\n";

	print AMBIGOUSSTATS $seq, "\t", $ambigous{$seq}{'clusters'}, "\n";

} #end printing ambigous LTRs

close(AMBIGOUSPROFILE);
close(AMBIGOUSSTATS);


#####################################

######### CLEANING CLUSTERS #########

#####################################

print "Cleaning Clusters\n";

open(PROFILES, ">$prefix.CleanedClusterProfiles.txt");
print PROFILES "ClusterID\t", join("\t", @Domains), "\n";

open(CLUSTERS, ">$prefix.CleanedClusteredSequences.txt");
print CLUSTERS "ClusterID\tSequences\t", join("\t", @Domains), "\n";

open(UNCLUSTERED, ">$prefix.CleanedUnclusteredSequences.txt");
print UNCLUSTERED "Sequences\t", join("\t", @Domains), "\n";

open(STATS, ">$prefix.CleanedClusterStats.txt");

foreach my $cluster ( sort {$a <=> $b} keys %clusters) {

	if(scalar(keys %{$clusters{$cluster}{'sequences'}}) == 0) { delete $clusters{$cluster} ; next

	} elsif (scalar(keys %{$clusters{$cluster}{'sequences'}}) == 1) {

		foreach my $seq (keys %{$clusters{$cluster}{'sequences'}} ) { print UNCLUSTERED "Cluster_", $cluster, "\t", $clusters{$cluster}{'sequences'}{$seq}, "\n"; delete $clusters{$cluster} }

	} else {

		my @profile;
		for(my $i = 0; $i < $numDomains; ++$i) {push(@profile, "NA") }

		foreach my $seq (keys %{$clusters{$cluster}{'sequences'}} ) {

			print CLUSTERS "Cluster_", $cluster, "\t", $clusters{$cluster}{'sequences'}{$seq}, "\n";

			#Get states
			my @states = split("\t", $clusters{$cluster}{'sequences'}{$seq});
			shift(@states);

			for(my $j = 0; $j < $numDomains; ++$j) { if($states[$j] ne "NA") { $profile[$j] = $states[$j] } }		

		} #end printing sequences and rebuilding the profile


		my $numValues = 0;
		print PROFILES "Cluster_", $cluster;
		foreach my $value (@profile) {
			if($value ne "NA") {++$numValues};
			print PROFILES "\t", $value;
		} #end printing domain clusters
		print PROFILES "\n";


		print STATS "Cluster_", $cluster, "\t", scalar(keys %{$clusters{$cluster}{'sequences'}}), "\t", $numValues, "\n";

	} #end printing output

} #end printing cluster profiles

close(PROFILES);
close(CLUSTERS);
close(UNCLUSTERED);
close(STATS);

print "Done cleaning the clusters\n\n";


print "Will split clusters\n";

#####################################

######## Splitting CLUSTERS #########

#####################################

my %splitClusters;

foreach my $cluster ( sort {$a <=> $b} keys %clusters) {

	my @seqs = keys (%{$clusters{$cluster}{'sequences'}});

	my $firstSeq = shift(@seqs);

	#Get states
	my @states = split("\t", $clusters{$cluster}{'sequences'}{$firstSeq});
	shift(@states);

	for(my $i = 0; $i < $numDomains; ++$i) { $clusters{$cluster}{$Domains[$i]} = $states[$i] }

	$splitClusters{$cluster}{$firstSeq} = $clusters{$cluster}{'sequences'}{$firstSeq};

	my @clusterIDs;
	push(@clusterIDs, $cluster);

	foreach my $restSeqs ( @seqs) {

		my $wasClustered = "no";

		#Get states
		my @states = split("\t", $clusters{$cluster}{'sequences'}{$restSeqs});
		my $seqName = shift(@states);

		foreach my $clusterToCheck (@clusterIDs) {

			my $matchCount = 0;

			for(my $j = 0; $j < $numDomains; ++$j) {

				if( ($clusters{$clusterToCheck}{$Domains[$j]} ne "NA") && ($states[$j] ne "NA") ) { ++$matchCount }

			} #end checking compatibility

			if( $matchCount != 0 ) {

				for(my $k = 0; $k < $numDomains; ++$k) { if($clusters{$clusterToCheck}{$Domains[$k]} eq "NA") { $clusters{$clusterToCheck}{$Domains[$k]} = $states[$k] } }

				$splitClusters{$clusterToCheck}{$restSeqs} = $clusters{$cluster}{'sequences'}{$restSeqs};

				$wasClustered = "yes";

			} #end if it has matches

		} #end checking clusters

		if($wasClustered eq "no") {

			++$clusterID;

			for(my $k = 0; $k < $numDomains; ++$k) { $clusters{$clusterID}{$Domains[$k]} = $states[$k] }

			push(@clusterIDs, $clusterID);

			$splitClusters{$clusterID}{$restSeqs} = $clusters{$cluster}{'sequences'}{$restSeqs};

		} #end if wasn't clustered

	} #end checking sequences

} #end splitting cluster

print "Done splitting clusters\n\n";




#####################################

###### PRINTING FINAL CLUSTERS ######

#####################################

print "Will print final clusters\n";

open(PROFILES, ">$prefix.FinalClusterProfiles.txt");
print PROFILES "ClusterID\t", join("\t", @Domains), "\n";

open(CLUSTERS, ">$prefix.FinalClusteredSequences.txt");
print CLUSTERS "ClusterID\tSequences\t", join("\t", @Domains), "\n";

open(UNCLUSTERED, ">$prefix.FinalUnclusteredSequences.txt");
print UNCLUSTERED "Sequences\t", join("\t", @Domains), "\n";

open(STATS, ">$prefix.FinalClusterStats.txt");

if( -d "$prefix-Clusters" ) { system("rm -r $prefix-Clusters") }

system("mkdir $prefix-Clusters");

foreach my $cluster ( sort {$a <=> $b} keys %splitClusters) {

	if (scalar(keys %{$splitClusters{$cluster}}) == 1) {

		foreach my $seq (keys %{$splitClusters{$cluster}} ) { print UNCLUSTERED "Cluster_", $cluster, "\t", $splitClusters{$cluster}{$seq}, "\n"; delete $splitClusters{$cluster} }

	} else {

		open(CLUSTER, ">$prefix-Clusters/Cluster_$cluster.txt");

		foreach my $seq (keys %{$splitClusters{$cluster}} ) { print CLUSTERS "Cluster_", $cluster, "\t", $splitClusters{$cluster}{$seq}, "\n"; print CLUSTER $seq, "\n" }

		my $numValues = 0;

		print PROFILES "Cluster_", $cluster;
		for(my $k = 0; $k < $numDomains; ++$k) { 

			if($clusters{$cluster}{$Domains[$k]} ne "NA") {++$numValues}
			print PROFILES "\t", $clusters{$cluster}{$Domains[$k]};

		} #end foreach domain

		print PROFILES "\n";


		print STATS "Cluster_", $cluster, "\t", scalar(keys %{$splitClusters{$cluster}}), "\t", $numValues, "\n";

	} #end printing output

} #end printing cluster profiles

close(PROFILES);
close(CLUSTERS);
close(UNCLUSTERED);

print "Done printing the final the clusters\n\n";

## DONE ##

print "ALL DONE!!!!\n\n";

exit;
















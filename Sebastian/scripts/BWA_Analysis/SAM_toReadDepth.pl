#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[2])) {

	print "Error, insufficient arguments. Please enter a SAM file and and output prefix\n";
	print "Usage: SAM_toReadDepth.pl <SAM FILE> <OUTPUT PREFIX>\n";
	die "Missing arguments\n";

} #end if

my $sam = $ARGV[0];
my $phred = $ARGV[1];
my $prefix = $ARGV[2];

# Initialize global variables to store information
my %reads;
my $unmapReads = 0;

my %coverages;
my %sequenceStats;

my $totalCount;
my $numPositions;

my $totalSqDev = 0;

my $max = 1;

open(SAM, $sam) or die "Can't open $sam\n";
while (my $samLine = <SAM>) {

	if($samLine =~ /^@/) {next}

	chomp($samLine);

	my @alignment = split("\t", $samLine);

	if ($alignment[2] ne "\*" && defined($alignment[9]) && $alignment[4] >= $phred) {

		++$reads{$alignment[0]}{'reference'}{$alignment[2]};
		$reads{$alignment[0]}{'score'}{$alignment[2]} += $alignment[4];

		my $end = $alignment[3] + length($alignment[9]);

		push(@{$reads{$alignment[0]}{'coordinates'}{$alignment[2]}}, [$alignment[3], $end]);


	} else {

		++$unmapReads;

	} #end of map/unmap

} #end while
close(SAM);

my $mapReads = scalar(keys %reads);

print STDERR "Done loading mapping information, will do coverage calculation\n\n";

foreach my $read (keys %reads) {

	if ( (scalar (keys %{$reads{$read}{'reference'}})) == 1) {

		my @ref = keys %{$reads{$read}{'reference'}};

#		if( $reads{$read}{'reference'}{$ref[0]} > 1) {

			foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

				for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { ++$coverages{$ref[0]}{$i} }
				
				++$sequenceStats{$ref[0]}{'NumMappedReads'};

			} #end foreach

#		} else {  } #split single mappings

	} elsif ( (scalar (keys %{$reads{$read}{'reference'}})) == 2) {

		my @ref = keys %{$reads{$read}{'reference'}};

		if ( $reads{$read}{'reference'}{$ref[0]} > $reads{$read}{'reference'}{$ref[1]} ) {

#			++$references{$ref[0]};
			foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

				for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { ++$coverages{$ref[0]}{$i} }
				
				++$sequenceStats{$ref[0]}{'NumMappedReads'};

			} #end foreach


		} elsif ( $reads{$read}{'reference'}{$ref[0]} < $reads{$read}{'reference'}{$ref[1]} ) {

#			++$references{$ref[1]};
			foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

				for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { ++$coverages{$ref[1]}{$i} }
				
				++$sequenceStats{$ref[1]}{'NumMappedReads'};

			} #end foreach

		} elsif ( $reads{$read}{'reference'}{$ref[0]} == $reads{$read}{'reference'}{$ref[1]} ) {


			if ( $reads{$read}{'score'}{$ref[0]} > $reads{$read}{'score'}{$ref[1]} ) {

#				if( $reads{$read}{'reference'}{$ref[0]} > 1) { 
#					++$references{$ref[0]}
					foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

						for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { ++$coverages{$ref[0]}{$i} }
				
						++$sequenceStats{$ref[0]}{'NumMappedReads'};

					} #end foreach
#				} else { 
#					++$singleEndMappedReferences{$ref[0]}
#				} #split single mappings

			} elsif ( $reads{$read}{'score'}{$ref[0]} < $reads{$read}{'score'}{$ref[1]} ) {

#				if( $reads{$read}{'reference'}{$ref[1]} > 1) {
#					++$references{$ref[1]}
					foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

						for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { ++$coverages{$ref[1]}{$i} }
				
						++$sequenceStats{$ref[1]}{'NumMappedReads'};

					} #end foreach
#				} else { 
#					++$singleEndMappedReferences{$ref[1]}
#				} #split single mappings

			} elsif ( $reads{$read}{'score'}{$ref[0]} == $reads{$read}{'score'}{$ref[1]} ) {

#				print EQUALLY $read, "\t", $ref[0], "\t", $ref[1], "\n";

#				if( $reads{$read}{'reference'}{$ref[0]} > 1) {
#					$references{$ref[0]} += 0.5;
#					$references{$ref[1]} += 0.5;
					foreach my $coordinate (@{$reads{$read}{'coordinates'}{$ref[0]}}) {

						for (my $i = $$coordinate[0]; $i <= $$coordinate[1]; ++$i) { $coverages{$ref[0]}{$i} += 0.5;$coverages{$ref[1]}{$i} += 0.5 }
				
						$sequenceStats{$ref[0]}{'NumMappedReads'} += 0.5;
						$sequenceStats{$ref[1]}{'NumMappedReads'} += 0.5;

					} #end foreach
#				} else {
#					$singleEndMappedReferences{$ref[0]} += 0.5;
#					$singleEndMappedReferences{$ref[1]} += 0.5;
#				} #split single mappings

			} #end elsif


		} #end elsif

	} #end multiple references

	delete($reads{$read});

} #end counting references

open(COVERAGES, ">$prefix.coveragePerPositions.txt") or die "Can't open file $prefix.coveragePerPositions.txt\n";
print COVERAGES "#Sequence\tPosition\tCoverage\n";

print STDERR "Done calculating coverages, will generate coverage file\n\n";

foreach my $reference (sort {$a cmp $b} keys %coverages) {

	my $totalCov;
	my $totalCovSequence;
	my $readsSequence;

	foreach my $position (sort {$b <=> $a} keys %{$coverages{$reference}}) {

		print COVERAGES $reference , "\t", $position, "\t", $coverages{$reference}{$position}, "\n";

		if($coverages{$reference}{$position} > $max) {$max = $coverages{$reference}{$position}}

		$totalCov += $coverages{$reference}{$position};

		++$totalCovSequence;

	} #end foreach

	my $averageSeq = $totalCov/ ($totalCovSequence-1);

	$totalCount += $totalCov;
	$numPositions += $totalCovSequence;

	$sequenceStats{$reference}{'averageSeq'} = $averageSeq;
	$sequenceStats{$reference}{'totalCovSequence'} = $totalCovSequence;

} #end counting references

close(COVERAGES);

my $average = $totalCount/ ($numPositions-1);

print STDERR "Done printing coverage file, will calculate statistics\n\n";

open(GENERAL, ">$prefix.generalStats.txt") or die "Can't open file $sam.generalStats.txt\n";

open(SEQUENCE, ">$prefix.statsPerSequence.txt") or die "Can't open file $sam.statsPerSequence.txt\n";
print SEQUENCE "#Sequence\tCoveredPositions\tAverageCoverage\tNumMappedReads\tStandardDeviation\n";

print GENERAL "Number of map pairs of reads:	", $mapReads, "\n";
print GENERAL "Number of unmap pairs of reads:	", $unmapReads/2, "\n\n";

print GENERAL "Total number of positions with reads: $numPositions\n";
print GENERAL "Average number of reads per positions: ", sprintf("%.3f", $average), "\n";
print GENERAL "Maximun number of reads for a single positions: $max\n";

foreach my $reference (sort {$a cmp $b} keys %coverages) {

	my $totalSqDevScaffold = 0;

	foreach my $position (sort {$a < $b} keys %{$coverages{$reference}}) {

		$totalSqDev += ($average - $coverages{$reference}{$position}) ** 2;

		$totalSqDevScaffold += ($sequenceStats{$reference}{'averageSeq'} - $coverages{$reference}{$position}) ** 2;

	} #end foreach

	my $stDevScaffold = ($totalSqDevScaffold / ($sequenceStats{$reference}{'totalCovSequence'} - 1) ) ** 0.5;

	print SEQUENCE $reference, "\t", sprintf("%.3f", $sequenceStats{$reference}{'averageSeq'}), "\t", $sequenceStats{$reference}{'totalCovSequence'}, "\t", $sequenceStats{$reference}{'NumMappedReads'}, "\t", , sprintf("%.4f", $stDevScaffold), "\n";

} #end foreach

my $stDev = ($totalSqDev / ($numPositions-1) ) ** 0.5;

print GENERAL "Standard deviation of the coverage across all the scaffolds: ", sprintf("%.4f", $stDev), "\n";

print STDERR "Completed!!\n\n";

close(COVERAGES);
exit;














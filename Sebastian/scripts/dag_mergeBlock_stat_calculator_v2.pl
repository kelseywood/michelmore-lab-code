#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Two are arguments are needed, please inpu them.\n";
	print "Dag merged synteny blocks file from SynMap\n";
        print "Outfile file.\n";

        exit 0;

} #end if

open(MATCHES, $ARGV[0]) or die "Can't open file $ARGV[0]\n";

my %distribution;
my $totalNumBlock = 0;
my $totalNumMatches = 0;
my @blocks;

my $firstLine = <MATCHES>;

while (my $line = <MATCHES>) {

	chomp($line);

        my @header = split("\t", $line);

        if ($line =~ /^#\d/) {

		my @blockInfo = split("\t", $line);

		my $numMatches = $blockInfo[5];

		++$distribution{$numMatches};

		push(@blocks, $numMatches);

		$totalNumMatches += $numMatches;

		++$totalNumBlock;

	} #end if is block header

} #end while

open(RESULTS, ">$ARGV[1].txt");

my $averMatchesperBlock = $totalNumMatches/$totalNumBlock;

print RESULTS "SourceFile\t$ARGV[0]\n\n";

print RESULTS "Total number of blocks\t$totalNumBlock\n";
print RESULTS "Total number of maches\t$totalNumMatches\n";
print RESULTS "Average number of matches per block\t", sprintf("%.2f", $averMatchesperBlock), "\n";

my $squareDif = 0;

foreach my $block (@blocks) {

	$squareDif += ($averMatchesperBlock - $block) * ($averMatchesperBlock - $block);

} #end foreach blocks

my $var = $squareDif / ($totalNumBlock - 1);

my $stdev = sqrt($var);

print RESULTS "Variance of the number of matches per block\t", sprintf("%.3f", $var), "\n";
print RESULTS "Standard deviation of the number of matches per block\t", sprintf("%.3f", $stdev), "\n";

print RESULTS "\nBlock sizes\n";

print RESULTS join(",", @blocks);

print RESULTS "\n\nDistribution of matches per block\n";
print RESULTS "BlockSize\tNumberofBlocks\n";

foreach my $key (sort { $a <=> $b } keys (%distribution) ) {

	print RESULTS $key, "\t", $distribution{$key}, "\n";

} #end forech

close(MATCHES);
close(RESULTS);

exit;

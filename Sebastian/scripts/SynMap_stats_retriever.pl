#!/usr/bin/perl

use strict; use warnings;

if (!defined($ARGV[0])) {

	print "No arguments imputted\n";
	print "Please input the directory where the data it's located\n";
	die "Missing arguments\n";

} #end if no arguments

opendir(DIR, $ARGV[0]);

open(RESULTS, ">$ARGV[0].summaryStats.txt");

while(my $file = readdir(DIR) ) {

	if( $file eq "." || $file eq "..") {next}

	if(-d "$ARGV[0]/$file") {

		print $file, "\n";

		opendir(SUBDIR, "$ARGV[0]/$file");

		my $totalhits = 0;
		my $blocks = 0;
		my $syntenicHits = 0;
		my $blocksMean = 0;
		my $blocksVar = 0;
		my $blocksSD = 0;

		while(my $subFile = readdir(SUBDIR) ) {

#			print "dag_mergeBlock_stat_calculator.pl $ARGV[0]/$file/$subFile $ARGV[0]/$file/$subFile.stats\n";

			if($subFile =~ /.filtered.txt/) {

				open(ALLHITS, "$ARGV[0]/$file/$subFile");

				while (my $allhitsline = <ALLHITS>) {

					++$totalhits;

				} #end of all hits files

				close(ALLHITS);

			} elsif ($subFile =~ /.aligncoords.Dm[0-9].ma[0-9].gcoords.txt$/) {

				if( system("dag_mergeBlock_stat_calculator.pl $ARGV[0]/$file/$subFile $ARGV[0]/$file/$subFile.stats") != 0) {

					print "Could not run dag_mergeBlock_stat_calculator.pl\n";

				} #end if running stats calculator					

				open(STATS, "$ARGV[0]/$file/$subFile.stats.txt");

				my $firstLine = <STATS>;
				my $secondLine = <STATS>;
				my $thirdLine = <STATS>;
				my $fourthLine = <STATS>;
				my $fithLine = <STATS>;
				my $sixthLine = <STATS>;
				my $seventhLine = <STATS>;

				chomp($thirdLine);
				chomp($fourthLine);
				chomp($fithLine);
				chomp($sixthLine);
				chomp($seventhLine);
				
				my @numBlocks = split("\t", $thirdLine);
				my @numSyntenicMatches = split("\t", $fourthLine);
				my @averMatchesPerBlock = split("\t", $fithLine);
				my @variance = split("\t", $sixthLine);
				my @sd = split("\t", $seventhLine);

				$blocks = $numBlocks[1];
				$syntenicHits = $numSyntenicMatches[1];
				$blocksMean = $averMatchesPerBlock[1];
				$blocksVar = $variance[1];
				$blocksSD = $sd[1];

			} #end if require file

		} #end while subdir

		print RESULTS $file, "\t", $totalhits, "\t", $blocks, "\t", $syntenicHits, "\t", $blocksMean, "\t", $blocksVar, "\t", $blocksSD, "\n";

	} #end if it's a dir

} #end while

close(RESULTS);

exit;

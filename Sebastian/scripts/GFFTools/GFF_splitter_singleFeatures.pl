#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

###################################################
#
#
#	GFF3 splitter
#
#
#
###################################################

# This version of the splitter it's design to deal with features that only have one line
# Features such as repeats that don't contain splice information or any parent/child relations

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to change the ID (ID is the first of the features).\n";
        print "Translation Table with sequence names to change.\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <Splitting ranges>\n";
        exit 0;

} #end if

my %fragments;

open(TRANSLATIONTABLE, $ARGV[1]) or die "Can't open $ARGV[1]";

while(my $codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        my @data = split("\t", $codeTranslation);

        $fragments{$data[0]}{$data[1]}{'start'} = $data[2];
        $fragments{$data[0]}{$data[1]}{'end'} = $data[3];

} #end while

open(GFF, $ARGV[0]) or die "Can't open $ARGV[0]\n";
open(UNCHANGED, ">$ARGV[0].unchanged.gff3");
open(SPLITTED, ">$ARGV[0].splitted.gff3");
open(SEGMENTED, ">$ARGV[0].segmented.gff3");
open(ERRORS, ">$ARGV[0].errors.txt");

while ( my $line = <GFF>) {

	if ($line =~ m/^#/) {print SPLITTED $line; next}

	chop($line);

        my @GFF_line = split ("\t", $line);

	if(!defined($fragments{$GFF_line[0]})) {print UNCHANGED $line, "\n"; next}

	my $reference;
	my $start;
	my $end;

	my @pieces;

	foreach my $fragment (keys %{$fragments{$GFF_line[0]}} ) {

		if ( ($GFF_line[3] >= $fragments{$GFF_line[0]}{$fragment}{'start'} ) && ($GFF_line[4] <= $fragments{$GFF_line[0]}{$fragment}{'end'}) ) {

			$reference = $fragment;
			$start = $GFF_line[3]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1);
			$end = $GFF_line[4]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1);

		} elsif ( ($GFF_line[3] < $fragments{$GFF_line[0]}{$fragment}{'start'} ) && ($GFF_line[4] > $fragments{$GFF_line[0]}{$fragment}{'start'}) ) {

			if($GFF_line[4] > $fragments{$GFF_line[0]}{$fragment}{'end'}) {

				push(@pieces, [$fragment, "1", $fragments{$GFF_line[0]}{$fragment}{'end'}-($fragments{$GFF_line[0]}{$fragment}{'start'}-1)] );

			} else {

				push(@pieces, [$fragment, "1", $GFF_line[4]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1)] );

			} #end else


		} elsif  ( ($GFF_line[3] < $fragments{$GFF_line[0]}{$fragment}{'end'} ) && ($GFF_line[4] > $fragments{$GFF_line[0]}{$fragment}{'end'}) ) {

			if($GFF_line[3] < $fragments{$GFF_line[0]}{$fragment}{'start'}) {

				push(@pieces, [$fragment, "1", $fragments{$GFF_line[0]}{$fragment}{'end'}-($fragments{$GFF_line[0]}{$fragment}{'start'}-1)] );

			} else {

				push(@pieces, [$fragment, $GFF_line[3]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1), $fragments{$GFF_line[0]}{$fragment}{'end'}-($fragments{$GFF_line[0]}{$fragment}{'start'}-1)] );

			} #end else

		} #end else


	} #end foreach

	if (defined($reference) ) { print SPLITTED $reference, "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $start, "\t", $end, "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n"}
	else {

		if (!defined($pieces[0])) { print ERRORS "Didn't found proper locations for $GFF_line[8]\n"}

		for(my $i = 0; $i < scalar(@pieces); ++$i) {
			my @linetoPrint = @GFF_line;
			my $subs = $linetoPrint[8] =~ s/ID=((A|B|C|D|E|F|G||H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z|a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z|0|1|2|3|4|5|6|7|8|9|_|-)*);{1}/ID=$1_$i;/;
			$linetoPrint[8] =~ s/Name=((A|B|C|D|E|F|G||H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z|a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z|0|1|2|3|4|5|6|7|8|9|_|-)*);{1}/Name=$1_$i;/;

			if($subs == 0) { print ERRORS "Couldn't find an ID to replace, please verify the characters on $linetoPrint[8], multiple features may end with the same ID\n"}

			$linetoPrint[0] = $pieces[$i][0];
			$linetoPrint[3] = $pieces[$i][1];
			$linetoPrint[4] = $pieces[$i][2];
			print SEGMENTED join("\t", @linetoPrint), "\n";
		} #end for each piece
	}



} #end while

close(GFF);

exit;

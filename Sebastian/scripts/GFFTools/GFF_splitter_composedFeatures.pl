#!/usr/bin/perl
use strict;
use warnings;

###################################################
#
#
#	GFF3 splitter
#
#
#
###################################################

# This version of the splitter it's design to deal with compose features
# Features such as gene models or splice aligmnents that contains parent/child relations

if (@ARGV < 2) {

	print "Two are arguments are needed, please inpu them.\n";
        print "Translation Table with sequence names to change.\n\n";
	print "At least one GFF File were to change the ID (ID is the first of the fields).\n";
	print "perl Fasta_renamer.pl <Splitting ranges> <GFFFile1> <GFFFile2> ... <GFFFileN>\n";
        exit 0;

} #end if

my %fragments;

my $splittingRanges = shift(@ARGV);

open(SPLITTINGRANGES, $splittingRanges) or die "Can't open $splittingRanges";

while(my $codeTranslation = <SPLITTINGRANGES>) {

	chop($codeTranslation);

        my @data = split("\t", $codeTranslation);

        $fragments{$data[0]}{$data[1]}{'start'} = $data[2];
        $fragments{$data[0]}{$data[1]}{'end'} = $data[3];

} #end while

my %segmentedIDs;

foreach my $file (@ARGV) {

	open(GFF, $file) or die "Can't open $file\n";
	open(UNCHANGED, ">$file.unchanged.gff3");
	open(SPLITTED, ">$file.splitted.gff3");
	open(SEGMENTED, ">$file.segmented.gff3");
	#open(ERRORS, ">$file.errors.txt");

	while ( my $line = <GFF>) {

		if ($line =~ m/^#/) {print UNCHANGED $line; next}

		chop($line);

		my @GFF_line = split ("\t", $line);

		if(!defined($fragments{$GFF_line[0]})) {print UNCHANGED $line, "\n"; next}

		my $reference;
		my $start;
		my $end;

		my $id;
		my $parent;

	
		my @features = split (";", $GFF_line[8]);

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "ID") { $id = $value[1] }
			elsif ($value[0] eq "Parent") { $parent = $value[1] } #end elsif

		} #end finding id/parent

		foreach my $fragment (keys %{$fragments{$GFF_line[0]}} ) {

			if ( ($GFF_line[3] >= $fragments{$GFF_line[0]}{$fragment}{'start'} ) && ($GFF_line[4] <= $fragments{$GFF_line[0]}{$fragment}{'end'}) ) {

				$reference = $fragment;
				$start = $GFF_line[3]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1);
				$end = $GFF_line[4]-($fragments{$GFF_line[0]}{$fragment}{'start'}-1);

			} # end if matches

		} #end foreach

		if(!defined($parent)) { $parent = $id }
		if(!defined($id)) { $id = $parent }

		if ( defined($reference) && !defined($segmentedIDs{$parent}) ) {
			print SPLITTED $reference, "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $start, "\t", $end, "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";
		} else {
			print SEGMENTED $line, "\n";
			$segmentedIDs{$id} = "0";


		}

	} #end while

	close(GFF);

} #end looping through files

exit;

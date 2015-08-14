#!/usr/bin/perl
use strict;
use warnings;

###################################################
#
#
#	GFF3 ID/Name/Parent Changer
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "Usage:\n";
	print "GFF_Score_changer.pl <Scores file> <Column Score> <GFFFile> <Type>\n\n";
	print " Scores file: tab delimited file which contains the score. First column should be the ID/Name/Parent of the gff feature.\n";
	print " Column score: which column it's the score located in the file (starting from column 1).\n";
	print " GFFFile: GFF file where to change the score.\n";
	print " Type: feature type to which change the score.\n\n";
        exit 0;

} #end if

my %scores = ();

open(SCORES, $ARGV[0]);

my $scoreColumn = $ARGV[1];

--$scoreColumn;

my $scoreCount;

while(my $score = <SCORES>) {

	chop($score);

        my @data = split("\t", $score);

        $scores{$data[0]} = $data[1];

	++$scoreCount;

} #end while

open(GFF, $ARGV[2]);

my $filename = $ARGV[2];

$filename =~ s/\.[^.]+$//;

open(RESULTS, ">$filename.reScored.gff");

my $type = $ARGV[3];

my $reScoredCount;

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	if ($GFF_line[2] ne $type) {

		print RESULTS $line, "\n";

	} else {

		my @features = split (";", $GFF_line[8]);

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "ID") {

				if(defined($scores{$value[1]})) { 
					$GFF_line[5] = $scores{$value[1]};
					++$reScoredCount;
					last;
				 } #end if defined

			} elsif ($value[0] eq "Name") {

				if(defined($scores{$value[1]})) { 
					$GFF_line[5] = $scores{$value[1]};
					++$reScoredCount;
					last;
				 } #end if defined

			} elsif ($value[0] eq "Parent") {

				if(defined($scores{$value[1]})) { 
					$GFF_line[5] = $scores{$value[1]};
					++$reScoredCount;
					last;
				 } #end if defined

			} #end else

		} #end foreach

		print RESULTS join("\t", @GFF_line);
		
		print RESULTS "\n";

	} #end else

} #end while

close(GFF);
close(RESULTS);

print "$scoreCount scores were readed it from the file\n";
print "Features which score was updated $reScoredCount\n";


exit;

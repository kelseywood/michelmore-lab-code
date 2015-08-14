#!/usr/bin/perl

###################################################
#
#
#	GFF3 Reference Sequence Changer
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to change the reference sequence.\n";
        print "Translation Table with sequence size to change.\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <TranslationTable>\n";
        exit 0;

} #end if

my %scaffolds = ();

open(TRANSLATIONTABLE, $ARGV[1]);

while($codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        @data = split("\t", $codeTranslation);

        $scaffolds{$data[0]} = $data[1];

} #end while

open(GFF, $ARGV[0]);
open(RESULTS, ">$ARGV[0].rename.gff");

while ($line = <GFF>) {

	chop($line);

        @GFF_line = split ("\t", $line);

	print RESULTS $scaffolds{$GFF_line[0]}, "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";

        if ($scaffolds{$GFF_line[0]} eq "") {

        	print "Error: Scaffold not found for ", $GFF_line[0],"\n";

        }

} #end while

close(GFF);
close(RESULTS);

exit;

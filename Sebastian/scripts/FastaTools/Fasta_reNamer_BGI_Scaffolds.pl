#!/usr/bin/perl

###################################################
#
#
#	Fasta Sequence Retriever
#
#
#
###################################################

#
# This script was prepared to split fasta files into
# smaller files containing a specific number of
# sequences
#

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "Fasta file with the sequences to look for.\n";
        print "Translation Table with sequence size to change.\n\n";
	print "perl Fasta_renamer.pl <fastaFile> <prefix>\n";
        exit 0;

} #end if

my %scaffolds = ();

open(TRANSLATIONTABLE, $ARGV[1]);

while($codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        @data = split("\t", $codeTranslation);

        $scaffolds{$data[0]} = [$data[1], $data[2]];

} #end while

open(FASTA, $ARGV[0]);
open(RESULTS, ">$ARGV[0].rename.txt");
open(FAILS, ">$ARGV[0].fails.txt");

while ($line = <FASTA>) {

	if ((substr($line,0,1)) eq ">") {

		chop($line);

                @scaffold_header = split (" ", $line);

        	$scaffold = substr($scaffold_header[0],1,100);

                print RESULTS ">", $scaffolds{$scaffold}[0], " ", $scaffolds{$scaffold}[1], "\n";

                if ($scaffolds{$scaffold}[0] eq "") {

                	print FAILS $scaffold, "\n";

                }

        } else {

                print RESULTS $line;

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;
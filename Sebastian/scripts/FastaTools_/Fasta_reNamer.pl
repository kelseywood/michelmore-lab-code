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

my %newNames = ();

open(TRANSLATIONTABLE, $ARGV[1]);

while($codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        @data = split("\t", $codeTranslation);

        $newNames{$data[0]} = [$data[1]];

} #end while

open(FASTA, $ARGV[0]);
open(RESULTS, ">$ARGV[0].rename.txt");
open(FAILS, ">$ARGV[0].fails.txt");

while ($line = <FASTA>) {

	if ($line =~ /^>/) {

		chomp($line);

                @header = split (" ", $line);

		$name = substr(shift(@header), 1, 300);
		
                print RESULTS ">", $newNames{$scaffold}[0], " ", join(" ", @header), "\n";

                if ($newNames{$name}[0] eq "") {

                	print FAILS $name, "\n";

                }

        } else {

                print RESULTS $line;

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;

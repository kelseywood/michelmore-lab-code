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
	print "Text file with list of GFF's to relocated (one per line), can be get with 'ls -1 *.gff'.\n";
        print "Translation Table with sequence size to change.\n\n";
	print "perl Fasta_renamer.pl <GFFList> <TranslationTable>\n";
        exit 0;

} #end if

my %scaffolds = ();

open(TRANSLATIONTABLE, $ARGV[1]);

while($codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        @data = split("\t", $codeTranslation);

        $scaffolds{$data[0]} = $data[1];

} #end while

open(LIST, $ARGV[0]);

while (my $file = <LIST>) {

	chomp($file);

	open(GFF, $file);
	open(RESULTS, ">$file.reLocated.gff");

	while ($line = <GFF>) {

		chomp($line);

		@GFF_line = split ("\t", $line);

		print RESULTS $scaffolds{$GFF_line[0]}, "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";

		if ($scaffolds{$GFF_line[0]} eq "") {

			print "Error: Scaffold not found for ", $GFF_line[0],"\n";

		}

	} #end while

	print "Done with $file\n\n";

	close(GFF);
	close(RESULTS);

} #end while files

close(LIST);

exit;

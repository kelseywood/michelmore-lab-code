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
        print "Translation Table with sequence size to change.\n\n";
	print "At least one GFF File were to change the reference sequence.\n";
	print "perl Fasta_renamer.pl <TranslationTable> <GFFFile1> <GFFFile2> ... <GFFFileN>\n";
        exit 0;

} #end if

my %scaffolds = ();

my $translationTable = shift(@ARGV);

open(TRANSLATIONTABLE, $translationTable);

while($codeTranslation = <TRANSLATIONTABLE>) {

	chop($codeTranslation);

        @data = split("\t", $codeTranslation);

        $scaffolds{$data[0]} = $data[1];

} #end while


foreach my $file (@ARGV) {

	print "Working on $file\n\n";

	open(GFF, $file);

	my $base = $file;

	$base =~ s{.*/}{};
	$base =~ s{\.[^.]+$}{};

	open(RESULTS, ">$base.reLocated.gff3");
	open(ERRORS, ">$base.reLocated.errors.txt");

	while ($line = <GFF>) {

		chop($line);

		if($line =~ /^#/) { print RESULTS $line }

		@GFF_line = split ("\t", $line);

		print RESULTS $scaffolds{$GFF_line[0]}, "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";

		if ($scaffolds{$GFF_line[0]} eq "") {

			print ERRORS "Error: Scaffold not found for ", $GFF_line[0],"\n";

		}

	} #end while

	close(GFF);
	close(RESULTS);
	close(ERRORS);

} #end looping through files

exit;

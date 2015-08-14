#!/usr/bin/perl

###################################################
#
#
#	Fasta Sequence retriever
#
#
#
###################################################

#
# This script was prepared to split fasta files into
# smaller files containing a specific number of
# sequences
#

if (@ARGV < 2) {

	print "Three are arguments are needed, please inpu them.\n";
	print "Fasta file with the sequences\n";
        print "Text file with the contigs to remove.\n";
        print "Outfile prefix.\n";

        exit 0;

} #end if

my $outputFile = 1;
my $toOutput = "";
my %contigs = "";

open(CONTIGS, $ARGV[1]);

while($contig = <CONTIGS>) {

	chop($contig);

        $contigs{$contig} = "1";

} #end while

open(FASTA, $ARGV[0]);
open(RESULTS, ">$ARGV[2].fasta");

while ($line = <FASTA>) {

	if ((substr($line,0,1)) eq ">") {

		chop($line);

        	$headerLine = substr($line,1,100);

                @header = split(" ", $headerLine);

                if($contigs{$header[0]} eq 1) {

			$toOutput = "FALSE";

                } else {

			print RESULTS $line;
                        print RESULTS "\n";

                        $toOutput = "TRUE";

                } #end else

        } else {

		if($toOutput eq "TRUE") {

                       print RESULTS $line;

                } #end else

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;

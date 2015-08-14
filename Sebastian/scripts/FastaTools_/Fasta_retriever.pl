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
	print "Fasta file where to retrieve the sequence\n";
        print "Text file with the contigs to retrieve.\n";
        print "Outfile prefix.\n";

        exit 0;

} #end if

my $outputFile = 1;
my $toOutput = "";
my %contigs = "";

open(CONTIGS, $ARGV[1]);

while($contig = <CONTIGS>) {

	chomp($contig);

	my @contigLine = split("\t", $contig);

        $contigs{$contigLine[0]} = "1";

} #end while

open(FASTA, $ARGV[0]);
open(RESULTS, ">$ARGV[2].fasta");

while ($line = <FASTA>) {

	if ((substr($line,0,1)) eq ">") {

		chop($line);

        	$headerLine = substr($line,1,100);

                @header = split(" ", $headerLine);

                if($contigs{$header[0]} eq 1) {

			delete($contigs{$header[0]});

			print RESULTS $line;
                        print RESULTS "\n";

                        $toOutput = "TRUE";

                } else {

			$toOutput = "FALSE";

                } #end else

        } else {

		if($toOutput eq "TRUE") {

                       print RESULTS $line;

                } #end else

        } #end else

} #end while

close(FASTA);
close(RESULTS);

foreach $key (%contigs) {

	print "$key not found\n";

} #end foreach

exit;

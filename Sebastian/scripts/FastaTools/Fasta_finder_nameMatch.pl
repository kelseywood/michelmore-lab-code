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
my @contigs;

open(CONTIGS, $ARGV[1]);

while($contig = <CONTIGS>) {

	chomp($contig);

        push (@contigs, $contig);

} #end while

my $searchString = join("|", @contigs);

print $searchString, "\n";

open(FASTA, $ARGV[0]);
open(RESULTS, ">$ARGV[2].fasta");

while ($line = <FASTA>) {

	if ($line =~ /^>/) {

		chomp($line);

                $toOutput = "FALSE";

#                foreach $search (@contigs) {

                	if( $line =~ /($searchString)/) { print RESULTS $line, "\n"; $toOutput = "TRUE" }

#                } #end foreach

        } else {

		if($toOutput eq "TRUE") {

                       print RESULTS $line;

                } #end else

        } #end else

} #end while

close(FASTA);
close(RESULTS);

exit;
#!/usr/bin/perl

######################################################
#                                                    #
#                Exonerate parser                    #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to pares the output from an exonerate analysis with the following options
# --showalignment no --showtargetgff yes --ryo "# --- START OF GENOMIC SEQUENCE DUMP ---\n>%ti_%tab-%tae  Genomic (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tas\n# --- STOP OF GENOMIC SEQUENCE DUMP ---\n# --- START OF CODING SEQUENCE DUMP ---\n>%ti_%tab-%tae  Coding (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tcs\n# --- STOP OF CONDING SEQUENCE DUMP ---\n"
# It will parse the file and create three fasta files and a gff file

use strict;
use warnings;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

my $CodonTable = Bio::Tools::CodonTable->new();

if(!defined($ARGV[0])) {
	print "Please provide an exonerate output\n";
	print "Exonerate should be run with the following options\n";
	print '--showalignment no --showtargetgff yes --ryo "# --- START OF GENOMIC SEQUENCE DUMP ---\n>%ti_%tab-%tae  Genomic (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tas\n# --- STOP OF GENOMIC SEQUENCE DUMP ---\n# --- START OF CODING SEQUENCE DUMP ---\n>%ti_%tab-%tae  Coding (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tcs\n# --- STOP OF CONDING SEQUENCE DUMP ---\n"\n';
	die "Missing file to parse\n";
} #end if not file provided

open(INPUT, $ARGV[0]);

open(GENOMIC, ">$ARGV[0].genomic.fasta");
open(CODING, ">$ARGV[0].coding.fasta");
open(PROTEIN, ">$ARGV[0].protein.fasta");
open(GFF, ">$ARGV[0].gff");

open(STATS, ">$ARGV[0].stats.txt");

print STATS "#Target	Query	Start	End	Strand	Rank	Identity	LengthAlignment\n";

while(my $line = <INPUT>) {

	if ($line =~ /--- START OF GFF DUMP ---/) {

		while(my $gffLine = <INPUT>) {

			if ($gffLine =~ /--- END OF GFF DUMP ---/) {;last}

			print GFF $gffLine;

		} #end while for gff

	} elsif ($line =~ /--- START OF GENOMIC SEQUENCE DUMP ---/) {

		while(my $genomicLine = <INPUT>) {

			if ($genomicLine =~ /--- STOP OF GENOMIC SEQUENCE DUMP ---/) {last}

			print GENOMIC $genomicLine;

		} #end while for genomic sequence

	} elsif ($line =~ /--- START OF CODING SEQUENCE DUMP ---/) {

		my $codingSequence;

		my $header;

		while(my $codingLine = <INPUT>) {

			if ($codingLine =~ /--- STOP OF CONDING SEQUENCE DUMP ---/) {

				$header =~ s/Coding/Protein/;

				print PROTEIN $header;

				my $codonStart = 0;

				my $aa = "gibberish";

				while ($aa ne "") {

					my $codon = substr($codingSequence, $codonStart, 3);

					$aa = $CodonTable->translate_strict($codon);

					print PROTEIN $aa;

					$codonStart += 3;

				} #end while

				print PROTEIN "\n";

				my @info = split(" ", $header);

				print STATS substr($info[0],1,300), "\t", $info[12], "\t", substr($info[2],1,30), "\t", substr($info[4],0,-1), "\t", $info[6], "\t", $info[8], "\t", $info[10], "\t", length($codingSequence), "\n";

				last

			} elsif ( !($codingLine =~ /^>/) && $line ne "") {chomp($codingLine); $codingSequence .= $codingLine

			} elsif ($codingLine =~ /^>/) {$header = $codingLine} #end for elsif

			print CODING $codingLine;

		} #end while for coding sequence

	} # End of looking for printing blocks

} #end while

close(INPUT);
close(GENOMIC);
close(CODING);
close(PROTEIN);
close(GFF);

exit;

#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

my $CodonTable = Bio::Tools::CodonTable->new();

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

if (@ARGV < 3) {

	print "Three are arguments are needed, please inpu them.\n";
	print "Fasta file where to retrieve the sequence\n";
        print "Text file with the coordinates to extract.\n";
        print "Flanking sequence.\n";
        print "Outfile prefix.\n";

        exit 0;

} #end if

my $idx = Bio::Index::Fasta->new(
	                         '-filename' => "$ARGV[0].idx",
	                         '-write_flag' => 1
	                        );

$idx->make_index("$ARGV[0]");

my $toOutput = $ARGV[3];

my $flank = $ARGV[2];

open(CONTIGS, $ARGV[1]);
open(OUTPUT, ">$toOutput.flank.$flank.fasta");

while(my $contig = <CONTIGS>) {

	chomp($contig);

	my @contigLine = split("\t", $contig);

	my $seq = $idx->fetch($contigLine[1]);

	my $start = $contigLine[2]-$flank;

	if($start < 1 ) {$start = 1}

	my $end = $contigLine[3]+$flank;

	if($end > ($seq->length) ) { $end = $seq->length}

	my $sequence = $seq->subseq($start, $end);

	print OUTPUT ">", $contigLine[0], " ", $contigLine[1], ":", $start, "..", $end, "\n";
	print OUTPUT $sequence, "\n";

} #end while


exit;

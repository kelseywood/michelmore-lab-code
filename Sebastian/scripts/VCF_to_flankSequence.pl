#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

###################################################
#
#
#	Fasta VCF to flanking sequence
#
#
#
###################################################

#
# This script was prepared to print a fasta file with the flanking sequence for SNP's using information from a VCF
#

if (@ARGV < 3) {

	print "Three are arguments are needed, please inpu them.\n";
	print "Fasta file where to retrieve the sequence\n";
        print "VCF file with the SNP information\n";
        print "Flanking sequence (number of bases to the sides).\n";
        print "Outfile prefix.\n";
        exit 0;

} #end if

my $idx = Bio::Index::Fasta->new(
	                         '-filename' => "$ARGV[0].idx",
	                         '-write_flag' => 1
	                        );

$idx->make_index("$ARGV[0]");

my $flank = $ARGV[2];

my $output = $ARGV[3];

open(VCF, $ARGV[1]);
open(FASTA, ">$output.flank.$flank.fasta");
open(INFO, ">$output.infoColumn.txt");

while(my $VCFline = <VCF>) {

	if ($VCFline =~ /^\#\#INFO/) {print INFO $VCFline, next}
	if ($VCFline =~ /^\#/) {next}

	chomp($VCFline);

	my @VCFinfo = split("\t", $VCFline);

	my $sequence = $VCFinfo[0];
	my $pos = $VCFinfo[1];
	my $ref = $VCFinfo[3];
	my $alt = $VCFinfo[4];
	my $score = $VCFinfo[5];
	my $info = $VCFinfo[7];

	my $seq = $idx->fetch($sequence);

	my $start = $pos-$flank;
	if($start < 1 ) {$start = 1}

	my $leftFlank = "";
	if($pos > 1 ) { $leftFlank = $seq->subseq($start, ($pos-1)) }

	my $end = $pos+$flank;
	if($end > ($seq->length) ) { $end = $seq->length}

	my $rightFlank = "";
	if($pos < ($seq->length) ) { $rightFlank = $seq->subseq(($pos+1), $end) }

	print FASTA ">", $sequence, "_", $pos, " score=", $score, " info=", $info, "\n";
	print FASTA $leftFlank, "[", $ref, "/", $alt, "]", $rightFlank, "\n";

} #end while

close(VCF);
close(INFO);
close(FASTA);

exit;





















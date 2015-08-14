#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

my $CodonTable = Bio::Tools::CodonTable->new();

###################################################
#
#
#	Fasta multi to one liner
#
#
#
###################################################

#
# This script was prepared to convert multiple line fasta to one liner fastas
#

if (@ARGV < 1) {

	print "One argument is needed, please input them.\n";
        print "Fasta file to translate.\n";

        exit 0;

} #end if

open(FASTA, $ARGV[0]) or die "Can't open file $ARGV[0]\n";
open(RESULTS, ">$ARGV[0].protein.fasta");

my $header = "";
my $sequence;

while (my $line = <FASTA>) {

	chomp($line);

	if ($line =~ /^>/) {

		if ($header eq "") {$header = $line; next};

		print RESULTS $header, "\n";

		my $codonStart = 0;

		my $aa = "gibberish";

		while ($aa ne "") {

			my $codon = substr($sequence, $codonStart, 3);

			$aa = $CodonTable->translate_strict($codon);

			print RESULTS $aa;

			$codonStart += 3;

                } #end while

		print RESULTS "\n";

		$header = $line;
		$sequence = "";

        } else {

		$sequence .= $line;

        } #end else

} #end while

print RESULTS ">", $header, "\n";

my $aa = "gibberish";

my $codonStart = 1;

while ($aa ne "") {

	my $codon = substr($sequence, $codonStart, 3);

	$aa = $CodonTable->translate_strict($codon);

	print RESULTS $aa;

	$codonStart += 3;

} #end while

close(FASTA);
close(RESULTS);

exit;







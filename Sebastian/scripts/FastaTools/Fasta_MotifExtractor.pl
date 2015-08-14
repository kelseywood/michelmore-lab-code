#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 4) {

	print "Please input four arguments\n";
	print "Fasta file where to find the motif\n";
	print "Motif to be found (X residue wildcard is valid)\n";
	print "Number of residues before to retrieve\n";
	print "Number of residues after to retrieve\n";
	print "Usage:\n";
	print "Fasta_MotifExtractor.pl <Fasta File> <Motif> <before> <after>\n\n";
	die "missing arguments\n\n";

} #end 

my $motif = $ARGV[1];
my $before = $ARGV[2];
my $after = $ARGV[3];

my $regex = $motif;

$regex =~ s/X|x/\\S/g;

my $fasta = $ARGV[0];

open(FASTA, $fasta);

while (my $fastaLine = <FASTA>) {

	chomp($fastaLine);

	my $seq = <FASTA>;

	chomp($seq);

	my $baseSeq = $seq;

	$seq =~ /$regex/g;

	if(defined($-[0])) {

		my $start;
		if( $-[0]-$before < 0 ) { $start = 0}
		else { $start = $-[0]-$before}
		my $length = ($+[0]-$-[0])+$before+$after;
		my $motifSeq = substr($baseSeq, $start, $length);

		print $fastaLine." ".$motif." from ".$-[0]." to ".$+[0]."\n";
		print $motifSeq, "\n";

	} #end if match found

#	my $i = 1;
#	while($seq =~ /$regex/g) {
#		my $start;
#		if( $-[0]-$before < 0 ) { $start = 0}
#		else { $start = $-[0]-$before}
#		my $length = ($+[0]-$-[0])+$before+$after;
#		my $motifSeq = substr($baseSeq, $start, $length);
#		print $fastaLine."_".$motif."_".$i." from ".$-[0]." to ".$+[0]."\n";
#		print $motifSeq, "\n";
#		++$i;
#	} #end while

} #end 

close(FASTA);

exit;

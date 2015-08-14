#!/usr/bin/perl

############################################################################
# contigs_to_sim_reads.pl                                                  #
# Cut long reads/contigs into synthetic paired-end reads with either       #
# forward-reverse orientation [pe] or reverse-forward orientation [mate].  #
# Use the provided step size to simulate "fragments".                      #
#                                                                          #
# Joan Wong                                                                #
############################################################################

use strict;
use warnings;

$#ARGV >= 5 || die "Usage: contigs_to_sim_reads.pl <pe|mate> <input.fasta> <frag size> <step size> <read length> <output prefix>  Optional: <offset> <genome size>\n";

my $dir = $ARGV[0];
my $input = $ARGV[1];
my $frag = $ARGV[2];
my $step = $ARGV[3];
my $readlen = $ARGV[4];
my $output = $ARGV[5].".sim_".$dir.$frag;
my $offset = $ARGV[6];
my $genome = $ARGV[7];

if (! $ARGV[6]) {$offset = 0;}
if($readlen > $frag) {die "Please select a read length that is shorter than the fragment size.\n";}
if($dir eq "pe" and $frag >= 1000) {print "\nWarning: Paired-end libraries typically have insert sizes < 1 kb.\n";}
elsif($dir eq "mate" and $frag < 1000) {print "\nNOTE: An insert size of at least 1 kb is recommended for simulating mate pairs.\n";}

open IN, "<$input";
open OUT1, ">$output\_1.fasta";
open OUT2, ">$output\_2.fasta";

my $header = <IN> || die "File cannot be read: $input.\n";	#Advance past the first header.
my $seqs_in = 0;
my $ID = 1;
my $read = "";
my $pos;
my $rev = "";

#Print mate pair "reads" from fragment ends in R-F orientation.
sub matereads {
	$seqs_in++;
	$pos = 0 + $offset;
	while($pos <= length($read)-$frag) {	#Stop when the remaining sequence is shorter than the desired fragment size.
		print OUT1 "\>$output\_frag\_$ID/1\n";	#Print forward reads to output file 1.  Direction = R.
		$rev = substr($read, $pos, $readlen);
		$rev =~ tr/ATGCNXatgcnx/TACGNXtacgnx/;
		$rev = reverse $rev;
		print OUT1 $rev, "\n";
		print OUT2 "\>$output\_frag\_$ID/2\n";	#Print reverse reads to output file 2.	Direction = F.
		print OUT2 substr($read, $pos+$frag-$readlen, $readlen), "\n";
		$pos = $pos+$step;
		$ID++;
	}
	$read = "";	#Clear variable.
}

#Print paired-end "reads" from fragment ends in F-R orientation.
sub pereads {
	$seqs_in++;
	$pos = 0 + $offset;
	while($pos <= length($read)-$frag) {	#Stop when the remaining sequence is shorter than the desired fragment size.
		print OUT1 "\>$output\_frag\_$ID/1\n";	#Print forward reads to output file 1.  Direction = F.
		print OUT1 substr($read, $pos, $readlen), "\n";
		print OUT2 "\>$output\_frag\_$ID/2\n";	#Print reverse reads to output file 2.	Direction = R.
		$rev = substr($read, $pos+$frag-$readlen, $readlen);
		$rev =~ tr/ATGCNXatgcnx/TACGNXtacgnx/;
		$rev = reverse $rev;
		print OUT2 $rev, "\n";
		$pos = $pos+$step;
		$ID++;
	}
	$read = "";	#Clear variable.
}

sub printreads {
	if ($dir eq "mate") {&matereads;}
	else {&pereads;}
}

while(<IN>) {
	chomp;
	if (not />/) {$read .= $_;}	#If not a header line, add it to the sequence variable.
	else {&printreads;}
}
&printreads;	#Process the last sequence in the file.

my $tot = 2*($ID-1)*$readlen;

print "\n=== $input ===\n";
print "Input sequences:   $seqs_in\n";
print "Output read pairs: ",$ID-1,"\n";
print "Total output size: ";
if ($tot < 1000) {print "$tot bp\n";}
elsif ($tot >= 1000 and $tot < 1000000) {printf("%.1f kb\n", $tot/1000);}
elsif ($tot >= 1000000 and $tot < 1000000000) {printf("%.1f Mb\n", $tot/1000000);}
else {printf("%.1f Gb\n", $tot/1000000000)}
if (defined $genome) {printf("Genome coverage:   %.1f\n", $tot/$genome);}
print "--> $output\_1.fasta\n";
print "--> $output\_2.fasta\n";

close IN;
close OUT1;
close OUT2;

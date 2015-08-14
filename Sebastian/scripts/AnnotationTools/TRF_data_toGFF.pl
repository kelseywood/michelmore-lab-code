#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) {
	print "Please provide an TRF .dat file\n";
	die "Missing arguments\n";
} #end if not defined argument

## Script design to convert an tRNAscan ss file into a gff and a fasta file

my $datFile = $ARGV[0];
my $base = $ARGV[0];
$base =~ s{\.[^.]+$}{};

my %SeqStr;

open(DAT, $datFile);

my $sequence;
my $TRindex = 1;

open(GFF, ">$base.gff3");

while(my $line = <DAT>) {

	chomp($line);

	if ($line =~ /^Sequence/) { $sequence = $line
				  
	} elsif ($line =~ /^\d/) {

# 1543 1582 19 2.2 19 86 4 55 70 10 7 12 1.35 AAAAATAAGTTCAAAACAA AAAAATAAGTCAAAACAAAAAAATCAGTTCAAAAGAAAAA

		my @repeatInfo = split(" ", $line);
	
		my @reference = split(" ", $sequence);
		my $start = $repeatInfo[0];
		my $end = $repeatInfo[1];
		my $score = $repeatInfo[7];

		my $PeriodSize = $repeatInfo[2];
		my $CopyNumber = $repeatInfo[3];
		my $PercentMatches = $repeatInfo[5];
		my $PercentIndels = $repeatInfo[6];
		my $Consensus = $repeatInfo[13];
		
		print GFF $reference[1], "\tTRF\tTandemRepeat\t", $start, "\t", $end,"\t", $score, "\t+\t.\t";
		print GFF "ID=TR", $TRindex, ";Name=TR", $TRindex, ";PeriodSize=", $PeriodSize, ";CopyNumber=", $CopyNumber, ";PercentMatches=", $PercentMatches, ";PercentIndels=", $PercentIndels, ";Consensus=", $Consensus, ";";
		print GFF "\n";

		++$TRindex

	} #end if is it's repeat

} #end reading ss file

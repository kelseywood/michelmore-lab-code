#!/usr/bin/perl

use strict;
use warnings;

if(!defined($ARGV[0])) { die "please input a fasta file\n\n"}

open(FASTAINPUT, $ARGV[0]) or die "Can't open $ARGV[0]\n\n";
open(OUTPUT, ">$ARGV[0].cleaned.fasta");

my $sequenceName;
my $sequence = "";

while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\r//;

	if( $line =~ /^>/) {

		if (!defined($sequenceName)) {
		        $sequenceName = $line
		} else {

#			print $sequenceName, "\n";

			my $nonACGT = $sequence;

			$nonACGT =~ s/(a|c|g|t|n|A|C|G|T|N)//g;

			if( length($nonACGT) == 0) {

			print OUTPUT $sequenceName, "\n";
			print OUTPUT $sequence, "\n";
			$sequenceName = $line;
			$sequence = "";

			} else {

			$sequenceName = $line;
			$sequence = "";

		        } #end else
			
		}

	} else {

		$sequence .= $line;

	} #end else

} #end while

			my $nonACGT = $sequence;

			$nonACGT =~ s/(a|c|g|t|n|A|C|G|T|N)//g;

			if( length($nonACGT) == 0) {

			print OUTPUT $sequenceName, "\n";
			print OUTPUT $sequence, "\n";

			} #end if print
			

close(FASTAINPUT);
close(OUTPUT);

exit;


__END__

#for ACGT counting

my $aCount += ($line =~ tr/aA/aA/);
my $tCount += ($line =~ tr/tT/tT/);
my $gCount += ($line =~ tr/gG/gG/);
my $cCount += ($line =~ tr/cC/cC/);
my $nCount += ($line =~ tr/nN/nN/);




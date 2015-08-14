#!/usr/bin/perl

use strict;
use warnings;

if (!(defined($ARGV[0])) ) {

	print "Please input Fasta file where to remove duplicates\n";
	die "Missing arguments\n";

} #end if missing arguments

open(FASTAINPUT, $ARGV[0]);

my $header;
my $sequence = "";

my %sequences;

while(my $line = <FASTAINPUT>) {

	if ($line =~ /^>/) {

		if( defined($header)) {

		  $sequences{$sequence}{$header} = "";


                } #end else

	        $header = $line;
		$sequence = "";

        } else {

		chomp($line);
		$sequence .= $line;

        } #end else

} #end while

$sequences{$sequence}{$header} = "";

close(FASTAINPUT);

open(OUTPUT, ">$ARGV[0]_withoutDuplicates.fasta");

foreach my $seq (keys %sequences) {

  my @headers = keys %{$sequences{$seq}};

  print OUTPUT $headers[0], " ", scalar(@headers), "\n";
  print OUTPUT $seq, "\n";

} #end foreach seq

close(OUTPUT);

exit;


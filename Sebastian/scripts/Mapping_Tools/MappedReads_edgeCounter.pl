#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[0])) {
	print "missing arguments\n";
	die;
} #end if no arguments


open(MPILEUP, $ARGV[0]);

open(COUNTER, ">$ARGV[0].edgesCounter.txt");

while(my $line = <MPILEUP>) {

	chomp($line);

	my @data = split("\t", $line);

	my $numberStarts = 0;
	my $numberEnds = 0;

	$numberStarts = ($data[4] =~ s/\^//g);
	$numberEnds = ($data[4] =~ s/\$//g);

	print COUNTER $data[0], "\t", $data[1], "\t", $numberStarts, "\t", $numberEnds, "\t", ($numberStarts + $numberEnds), "\n";

} #end reading genotypes

close(MPILEUP);

exit;

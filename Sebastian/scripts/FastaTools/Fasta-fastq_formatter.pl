#!/usr/bin/perl
use warnings;
use strict;

if(!defined($ARGV[0])) { die "Missing arguments, please provide a fasta file\n"}

#open DP File
open(INPUT, $ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # Can't open fasta file to convert

my $file = $ARGV[0];

$file =~ s{\.[^.]+$}{};

open(RESULTS, ">$file.fastq") || die "Cannot open file results file";

my $header;
my $sequence = "";

while (my $line = <INPUT>) {

	chop ($line);

	if($line =~ m/^>/ && defined($header) ) {

		my @quality = ();
		my $length = length($sequence);

	   	for(my $i = 0; $i < $length; ++$i) {
			push (@quality, "I");
		} #end for
		
		print RESULTS "@", $header, "\n";
		print RESULTS $sequence, "\n";
		print RESULTS "+\n";
		print RESULTS @quality, "\n";

		$header = substr($line,1,400);

		$sequence = "";

	} elsif (defined($header) && ($line ne "") ) { $sequence .= $line

	} elsif(!defined($header) ) { $header = substr($line,1,400)

	} #end elsif


} #end while


my @quality;
my $length = length($sequence);

for(my $i = 0; $i < $length; ++$i) {
	push (@quality, "I");
} #end for
		
print RESULTS "@", $header, "\n";
print RESULTS $sequence, "\n";
print RESULTS "+", $header, "\n";
print RESULTS @quality, "\n";

close(INPUT);
close(RESULTS);

exit;

#!/usr/bin/perl

use warnings;use strict;

if(!defined($ARGV[0])) { die "Missing arguments, please provide a fasta file\n"}

#open DP File
open(INPUT, $ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # Can't open fasta file to convert

my $file = $ARGV[0];

$file =~ s{\.[^.]+$}{};


my $quality = "I";

my $prevHeaderLine = ();
my @prevQuality = ();

open(RESULTS, ">$file.fastq") || die "Cannot open file results file";

my $line = <INPUT>;
chop ($line);
print RESULTS $line, "\n";

$prevHeaderLine = substr($line,1,100);

$line = <INPUT>;
chop ($line);
print RESULTS $line, "\n";

my $length1 = length($line);

for(my $i=0; $i < $length1; ++$i) {
	push (@prevQuality, $quality);
} #end for

while ($line = <INPUT>) {

	chop ($line);

	my @data = split("", $line);

	if($data[0] eq ">") {

		print RESULTS "@", $prevHeaderLine, "\n";
		print RESULTS @prevQuality, "\n";
		print RESULTS $line, "\n";

		$prevHeaderLine = substr($line,1,100);

	} elsif ($data[0] ne "") {

		print RESULTS $line, "\n";

		my $length = length($line);

		@prevQuality = "";

		for(my $i=0; $i < $length; ++$i) {
			push (@prevQuality, $quality);
		} #end for

	} #end ifelse

} #end while

print RESULTS "@", $prevHeaderLine, "\n";
print RESULTS @prevQuality, "\n";
print RESULTS $line, "\n";


close(INPUT);
close(RESULTS);

exit;

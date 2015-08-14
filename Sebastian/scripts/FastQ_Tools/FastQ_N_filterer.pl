#!/usr/bin/perl

use strict; use warnings;

if( !defined($ARGV[1]) ) {
	print "Please provide a fastq file and a maximum number of N's allowed per read\n";
	print "FastQ_N_filterer <Fastq file> <Max N per read>\n";
	die "missing argumetns\n"
}

open(FASTQ, $ARGV[0]);

my $file = $ARGV[0];

my $maxN = $ARGV[1];

$file =~ s{\.[^.]+$}{};

open(RESULTS, ">$file.readsWithNremoved.fastq");

while(my $header = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	chomp($sequence);

	my $length = length($sequence);

	my $numN = $sequence =~ tr/Nn/Nn/;

	if ($numN < $maxN ) {
		print RESULTS $header;
		print RESULTS $sequence, "\n";
		print RESULTS $header2;
		print RESULTS $qual;
	} #end if print

} #end while


close(FASTQ);
close(RESULTS);
exit;

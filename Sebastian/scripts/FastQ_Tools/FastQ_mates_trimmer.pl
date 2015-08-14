#!/usr/bin/perl

use strict; use warnings;

if( !defined($ARGV[1]) ) {
	print "Please provide a read size to trim and at least one fastq file\n";
	print "FastQ_mates_trimmer <read size> <Fastq file1> <Fastq file2> ...  <Fastq fileN>\n";
	die "missing argumetns\n"
}

my $length = shift(@ARGV);

foreach my $file (@ARGV) {

	open(FASTQ, $file);

	$file =~ s{\.[^.]+$}{};

	open(TRIM, ">$file.trimmedReadsto$length.fastq");

	while(my $header = <FASTQ>) {

		my $sequence = <FASTQ>;
		my $header2 = <FASTQ>;
		my $qual = <FASTQ>;

		if (!defined($qual)) {die "Error in format of fastq file missing quality\n"}

		chomp($sequence);
		chomp($qual);

		my $newSeq = substr($sequence, 0, $length);
		my $newqual = substr($qual, 0, $length);

		print TRIM $header;
		print TRIM $newSeq, "\n";
		print TRIM $header2;
		print TRIM $newqual, "\n";

	} #end while

	close(FASTQ);
	close(TRIM);

} #end looping through files

exit;

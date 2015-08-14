#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Two arguments are need it please input them.\n";
	print "AllPaths final_summary file.\n";
	print "Output prefix for GFF file.\n";
	print "Usage: .\n";
	die;

} #end if

my $file = $ARGV[0];
my $prefix = $ARGV[1];

open(INPUT, $file) or die "Can't open file $file\n";

open(OUTPUT, ">$prefix.gff");

my $scaffold;

my $start = 1;

while(my $line = <INPUT>) {

	if ($line =~ /^scaffold/) {

		chomp($line);

		my @data = split(" ", $line);

		$scaffold = join("_", ($data[0], $data[1]));

	} elsif ( $line =~ /^(1|2|3|4|5|6|7|9|0)/ ) {

		chomp($line);

		$line =~ s/\(//g;
		$line =~ s/\)//g;

		my @data = split(" ", $line);

		$start = 1;

		print OUTPUT $scaffold, "\t", "AllPath", "\t", "Contig", "\t", $start, "\t", $data[3], "\t", ".", "\t", ".", "\t", ".", "\t", "ID=Contig_", $data[0], ";Name=Contig_", $data[0], ";Length=", $data[3], "\n"; 		

		$start += $data[3];
	
	} elsif ( $line =~ /^ -- / ) {
		
		chomp($line);

		$line =~ s/\(//g;
		$line =~ s/\)//g;

		my @data = split(" ", $line);

		$start += $data[1];

		print OUTPUT $scaffold, "\t", "AllPath", "\t", "Contig", "\t", $start, "\t", $start+$data[8]-1, "\t", ".", "\t", ".", "\t", ".", "\t", "ID=Contig_", $data[5], ";Name=Contig_", $data[5], ";Length=", $data[8], "\n"; 		

		$start += $data[8];

	} #end

} #end reading file






















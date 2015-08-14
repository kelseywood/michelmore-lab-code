#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Two arguments are need it please input them.\n";
	print "File to convert.\n";
	print "Source to pass to the GFF file.\n";
	print "Ussage: SlidingWindow_Sumarizer.pl <File> <WindownSize> <Column>.\n";
	die;

} #end if

my $source = $ARGV[1];

open(TABLE, $ARGV[0]);

open(GFF, ">$ARGV[0].gff");

while(my $line = <TABLE>) {

	chomp($line);

	my @data = split("\t", $line);

	print GFF $data[0], "\t", $source, "\tmatch\t", $data[1], "\t", $data[2], "\t", $data[3], "\t.\t.\tID=", $data[0],"-", $data[1],"_", $data[2],";Name=", $data[0],"-", $data[1],"_", $data[2],";\n";

} #end while

close(TABLE);
close(GFF);

exit;

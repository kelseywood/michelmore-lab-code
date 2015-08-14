#!/usr/bin/perl
use strict; use warnings;

if(@ARGV < 1) {

	print "Missing arguments, please input fasta file and prefix for read name\n";

	die "Missing arguments\n";

} #end if arguments included


open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_reName.txt");

my $base = $ARGV[1];

my $readID = 1;

while(my $line = <FASTAINPUT>) {

	chomp($line);

        my $sequence = <FASTAINPUT>;

        chomp($sequence);

        print OUTPUT ">",$base,"_",$readID, "\n";
        print OUTPUT $sequence, "\n";

	++$readID;

} #end while

close(FASTAINPUT);
close(OUTPUT);

exit;

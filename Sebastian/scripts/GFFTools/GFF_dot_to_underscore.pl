#!/usr/bin/perl
use strict; use warnings;
my @chrs = qw(0 1 2 3 4 5 6 7 8 9 A);

open(GFFINPUT, $ARGV[0]);
open(GFFOUTPUT, ">$ARGV[0].newName.gff3");

while(my $line = <GFFINPUT>) {

	foreach my $chr (@chrs) {

		my $orgName = "Lsat.1.v3.g.$chr\.";

		my $newName = "Lsat_1_v3_g_$chr\_";

		$line =~ s/$orgName/$newName/g;

	} #end foreach

	print GFFOUTPUT $line;

} #end while

close(GFFOUTPUT);
close(GFFINPUT);

exit;



















#!/usr/bin/perl
use strict;
use warnings;

if (!defined($ARGV[0])) {

	print "Error, not groups.txt input was entered. Please enter a orthoMCL output file for analysis\n";
	die "Missing arguments\n";

} #end if

open(GROUPS, $ARGV[0]);
open(LIST, ">$ARGV[0].list.txt");


while(my $line = <GROUPS>) {

	chomp($line);

	my @group = split(" ", $line);

	my $groupName = shift(@group);

	$groupName =~ s/://;

	foreach my $member (@group) {

		print LIST $groupName, "\t", $member, "\n";

	} #end foreach member

} #end while file

close(GROUPS);
close(LIST);

exit;

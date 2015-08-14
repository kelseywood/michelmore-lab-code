#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[0])) {
	print "missing arguments\n";
	print "Loc file\n";

	die;
} #end if no arguments

my %markers;

open(LOC, $ARGV[0]);
open(CLEANED, ">$ARGV[0].duplicatedRemoved.loc");
open(BINS, ">$ARGV[0].duplicatedMarkers.txt");

my $header = "true";

while(my $line = <LOC>) {

	if($header eq "true") {

		print CLEANED $line;

		if($line =~ /^;/) {$header = "false"}

	} else {

        	chomp($line);

		my @data = split("\t", $line);

		my $marker = shift(@data);

		my $haplotype = join("\t", @data);

		push(@{$markers{$haplotype}}, $marker);

#                print $haplotype, "\n";

	} #end if header

} #end reading genotypes

close(LOC);

print "Number of unique haplotypes found ", scalar(keys %markers), "\n";

foreach my $haplotype (keys %markers) {

      	my @group = @{$markers{$haplotype}};

	my $head = shift(@group);

        print CLEANED $head, "\t", $haplotype, "\n";

        foreach my $node (@group) { print BINS $head, "\t", $node, "\n"}

} #end

exit;

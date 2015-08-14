#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	print "missing arguments\n";
	print "Loc file, tree clust file and prefix for output files required\n";
	die;
} #end if no arguments

my %markers;

open(LOC, $ARGV[0]);

my $header = "true";
my @individualNames;

while(my $line = <LOC>) {

	if($header eq "true") {

		if($line =~ /^;/) {
			$header = "false";
                        @individualNames = split("\t", $line);
                        shift(@individualNames);
		} #end if

	} else {

		my @data = split("\t", $line);

		my $missings = $line =~ tr/-//;

		my $name = shift(@data);

		$markers{$name}{'missings'} = $missings;
		$markers{$name}{'genotypes'} = join("\t", @data);

	} #end if header

} #end reading genotypes

close(LOC);

my %bins;

open(CLUST, $ARGV[1]);

while (my $node = <CLUST>) {

	my @info = split("\t", $node);

	push( @{$bins{$info[16]}{'markers'}}, $info[25]);

        if( !defined($bins{$info[16]}{'head'}) || ( $markers{$bins{$info[16]}{'head'}}{'missings'} > $markers{$info[25]}{'missings'} ) ) {

	        $bins{$info[16]}{'head'} = $info[25];

        } #end if defining head

} #end reading map

close(CLUST);

open(BINS, ">$ARGV[2].bins.txt");
open(REPRESENTATIVES, ">$ARGV[2].representatives.loc");

print REPRESENTATIVES ";";

foreach my $individual( @individualNames ) {

	print REPRESENTATIVES "\t", $individual;

} #end foreach

my $numGroups = 0;

foreach my $group (keys %bins) {

	++$numGroups;

	my $head = $bins{$group}{'head'};

        print REPRESENTATIVES $head, "\t", $markers{$head}{'genotypes'};

        foreach my $marker (@{$bins{$group}{'markers'}} ) {

        	print BINS $bins{$group}{'head'}, "\t", $group, "\t", $marker, "\n";

        } #end foreach printing markers

} #end foreach printing bins

close(BINS);
close(REPRESENTATIVES);

print "Number of bins found $numGroups\n";

exit;

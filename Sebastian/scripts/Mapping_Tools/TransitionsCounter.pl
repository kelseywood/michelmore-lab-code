#!/usr/bin/perl

use strict;
use warnings;

if ( !defined($ARGV[1])) {
	print "missing arguments\n";
	die;
} #end if no arguments


open(LOC, $ARGV[0]);

open(TRANSITIONS, ">$ARGV[0].transitions.txt");
open(SUSPICIOUS, ">$ARGV[0].putativeBreakpoints.$ARGV[1]orMoreTransitions.txt");

my $header = "true";
my @currentState;
my $curScaffold = "";
my $numInds;

while(my $line = <LOC>) {

	if($header eq "true") {

        	$header = "false";
		my @headerData = split("\t", $line);

                shift(@headerData);
                shift(@headerData);
#                shift(@headerData);
#                shift(@headerData);

                foreach my $ind (@headerData) { push(@currentState, "-")}

	} else {

        	chomp($line);

                $line =~ s/(H|h|U|u)/-/g;

		my @data = split("\t", $line);

                my $markerScaffold = shift(@data);
                my $markerStart = shift(@data);
#                my $markerEnd = shift(@data);
#                my $markerSNPs = shift(@data);



                my $transitions = 0;

                my $indID = 0;

                foreach my $datapoint (@data) {

                	if( $curScaffold ne $markerScaffold ) { $currentState[$indID] = $datapoint}

                	elsif ($datapoint eq "-") { }

                        elsif ($currentState[$indID] eq "-") { $currentState[$indID] = $datapoint }

                        elsif ($currentState[$indID] eq $datapoint) { }

                        elsif ($currentState[$indID] ne $datapoint) { $currentState[$indID] = $datapoint; ++$transitions }

                        ++$indID;

                } #end loading hash

                $curScaffold = $markerScaffold;

#                print TRANSITIONS join("\t", ($markerScaffold, $markerStart, $markerEnd, $markerSNPs, $transitions)), "\n";
#                if($transitions >= $ARGV[1]) { print SUSPICIOUS join("\t", ($markerScaffold, $markerStart, $markerEnd, $markerSNPs, $transitions)), "\n" }

                print TRANSITIONS join("\t", ($markerScaffold, $markerStart, $transitions)), "\n";
                if($transitions >= $ARGV[1]) { print SUSPICIOUS join("\t", ($markerScaffold, $markerStart, $transitions)), "\n" }

	} #end if header

} #end reading genotypes

close(LOC);

exit;

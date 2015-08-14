#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[1])) {

	print "Please input a fasta file to trim and length of seed\n";

	die "Missing input file\n";

} #end if

my $forward = $ARGV[0];

my $seed = $ARGV[1];

open(FORWARD, $forward) or die "Can't open $forward\n";

my %data;

	print "Processing $forward file.\n";

	my $lineCounter = 0;

        while(my $lineForward = <FORWARD>) {

		chomp($lineForward);

        	my $firstChar = substr($lineForward, 0,1);

                if ($firstChar eq "@") {

			my $progress = $lineCounter/100000;

			if (int($progress) eq $progress) {

				print "$lineCounter reads has been analyzed.\n";

			} #end if

			my $header = substr($lineForward, 0,-2);

                        my $sequenceForward = <FORWARD>;

			my $length = length($sequenceForward);

                        my $sequence = substr($sequenceForward,0,$seed)."~".substr($sequenceForward,($length-$seed),$seed);

                        $data{$sequence} = $header;

			++$lineCounter;

#			print $sequence, "\n";

		} #end if 

        } #end while

	print $lineCounter, "\n";

	close(FORWARD);

	open(OUTPUT, ">output.txt");

	print "Loading unique hash done, selecting unique IDs.\n";

	my %uniqReads;

	foreach my $uniqSequence (keys(%data)) {

		print OUTPUT $data{$uniqSequence}, "\n";

		$uniqReads{$data{$uniqSequence}} = 1;

	} #end foreach

	print scalar(keys(%data) ), "\n";

	open(FORWARD, $forward);

	open(FORWARDUNIQ, ">$forward.uniqReads.fq");

	print "Printing into unique $forward file.\n";
	
	while(my $lineForward = <FORWARD>) {

		chomp($lineForward);

        	my $firstChar = substr($lineForward, 0,1);

                if ($firstChar eq "@") {

			my $header = substr($lineForward, 0,-2);

			if(defined($uniqReads{$header})) {

				print FORWARDUNIQ $lineForward, "\n";

				my $sequence = <FORWARD>;
				chomp($sequence);
				print FORWARDUNIQ $sequence, "\n";

				my $secondHeader = <FORWARD>;
				chomp($secondHeader);
				print FORWARDUNIQ $secondHeader, "\n";

				my $quality = <FORWARD>;
				chomp($quality);		
				print FORWARDUNIQ $quality, "\n";

			} #end if
		} #end if
	} #end while

exit;


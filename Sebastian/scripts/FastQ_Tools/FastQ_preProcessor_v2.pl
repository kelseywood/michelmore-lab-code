#!/usr/bin/perl
use strict;
use warnings;

my @files;

open(INPUT, $ARGV[0]);

open(LOG, ">$ARGV[0].log.txt");

while(my $forward = <INPUT>) {

	chomp($forward);

        my $reverse = <INPUT>;

        chomp($reverse);

	open(FORWARD, $forward);

        open(REVERSE, $reverse);

        my %data;

	print "Processing $forward and $reverse files.\n";
	print  LOG "Processing $forward and $reverse files.\n";

	my $lineCounter = 0;

        while(my $lineForward = <FORWARD>) {

		chomp($lineForward);

        	my $lineReverse = <REVERSE>;

		chomp($lineReverse);

        	my $firstChar = substr($lineForward, 0,1);

                if ($firstChar eq "@") {

			my $progress = $lineCounter/100000;

			if (int($progress) eq $progress) {

				print "$lineCounter reads has been analyzed.\n";

			} #end if

			my $header = substr($lineForward, 0,-2);

                        my $sequenceForward = <FORWARD>;

                        my $sequenceReverse = <REVERSE>;

                        my $sequence = substr($sequenceForward,1,40).substr($sequenceReverse,1,40);

                        $data{$sequence} = $header;

			++$lineCounter;

		} #end if 

        } #end while

	print $lineCounter, "\n";

	print LOG "Total number of reads $lineCounter\n";

	close(FORWARD);
	close(REVERSE);

	open(OUTPUT, ">output.txt");

	print "Loading unique hash done, selecting unique IDs.\n";

	my %uniqReads;

	foreach my $uniqSequence (keys(%data)) {

		print OUTPUT $data{$uniqSequence}, "\n";

		$uniqReads{$data{$uniqSequence}} = 1;

	} #end foreach

	print scalar(keys(%data) ), "\n";

	print LOG "Number of uniq reads ", scalar(keys(%data) ), "\n";

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

	open(REVERSE, $reverse);
	open(REVERSEUNIQ, ">$reverse.uniqReads.fq");
	
	print "Printing into unique $reverse file.\n";

	while(my $lineReverse = <REVERSE>) {

		chomp($lineReverse);

        	my $firstChar = substr($lineReverse, 0,1);

                if ($firstChar eq "@") {

			my $header = substr($lineReverse, 0,-2);

			if(defined($uniqReads{$header})) {

				print REVERSEUNIQ $lineReverse, "\n";

				my $sequence = <REVERSE>;
				chomp($sequence);
				print REVERSEUNIQ $sequence, "\n";

				my $secondHeader = <REVERSE>;
				chomp($secondHeader);
				print REVERSEUNIQ $secondHeader, "\n";

				my $quality = <REVERSE>;
				chomp($quality);
				print REVERSEUNIQ $quality, "\n";

			} #end if
		} #end if
	} #end while


} #end while for input files

exit;


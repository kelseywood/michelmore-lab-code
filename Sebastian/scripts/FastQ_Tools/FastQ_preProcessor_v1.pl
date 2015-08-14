#!/usr/bin/perl
use strict;
use warnings;

my @files;

open(INPUT, $ARGV[0]);

while(my $forward = <INPUT>) {

	chomp($forward);

        my $reverse = <INPUT>;

        chomp($reverse);

	open(FORWARD, $forward);

        open(REVERSE, $reverse);

        my %data;

	print "Processing $forward and $reverse files.\n";

	my $lineCounter = 0;

        while(my $lineForward = <FORWARD>) {

		my $progress = $lineCounter/100000;

		if (int($progress) eq $progress) {

			print "$lineCounter reads has been analyzed.\n";

		} #end if
	
		chomp($lineForward);

        	my $lineReverse = <REVERSE>;

		chomp($lineReverse);

        	my $firstChar = substr($lineForward, 0,1);

                if ($firstChar eq "@") {

			my $header = substr($lineForward, 0,-2);

                        my $sequenceForward = <FORWARD>;

                        my $sequenceReverse = <REVERSE>;

                        my $sequence = substr($sequenceForward,1,40).substr($sequenceReverse,1,40);

                        my $headerQualForward = <FORWARD>;

                        my $headerQualReverse = <REVERSE>;

                        my $qualForward = <FORWARD>;

                        my $qualReverse = <REVERSE>;

                        my $qual = substr($qualForward,1,40).substr($qualReverse,1,40);

                        my $qualValue = 0;

                        foreach my $character (split("", $qual) ) {

                        	$qualValue += ord($character);

                        } #end foreach qualities

                        $data{$sequence}{$qualValue} = $header;

		} #end if 

		++$lineCounter;

        } #end while

	print $lineCounter, "\n";

	close(FORWARD);
	close(REVERSE);

	open(OUTPUT, ">output.txt");

	print "Loading unique hash done, selecting unique ID's.\n";

	my %uniqReads;

	foreach my $uniqSequence (keys(%data)) {

		my @quals = sort (keys %{$data{$uniqSequence}} );
#		my @quals = keys %{$data{$uniqSequence}};

		my $maxQual = pop(@quals);

		print OUTPUT $data{$uniqSequence}{$maxQual}, "\n";

		$uniqReads{$data{$uniqSequence}{$maxQual}} = 1;

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


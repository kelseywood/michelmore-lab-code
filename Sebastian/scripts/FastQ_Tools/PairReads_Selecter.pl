#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV < 1) {

	print "No input value entered.\n";
	print "Please input a fofn (File Of File Names), list of fastq files, one file per line (Paired FastQ files should be one after the other)/\n";
	die "Missing Arguments";

} #end if

my @files;

open(INPUT, $ARGV[0]);

open(LOG, ">$ARGV[0].pairSelecter.log.txt");

open(PAIREDFILES, ">$ARGV[0].PairedFiles.txt");
open(SINGLEFILES, ">$ARGV[0].SingleFiles.txt");

while(my $forward = <INPUT>) {

	chomp($forward);

        my $reverse = <INPUT>;

        chomp($reverse);

	print "Processing $forward and $reverse files.\n";
	print  LOG "Processing $forward and $reverse files.\n";

	my %reads = ();

	my $lineCounter = 0;

	open(FORWARD, $forward);

        while(my $lineForward = <FORWARD>) {

#       	my $firstChar = substr($lineForward, 0,1);

		if ($lineForward =~ /^@/) {

#                if ($firstChar eq "@") {

			chomp($lineForward);

			my $name = substr($lineForward, 1, -2);

#			print $name, "\n";

			$reads{$name} = "S";

			++$lineCounter;

		} #end if 

        } #end while

	print "Total number of forward reads $lineCounter\n";

	print LOG "Total number of forward reads $lineCounter\n";

	close(FORWARD);



	$lineCounter = 0;

        open(REVERSE, $reverse);

        while(my $lineReverse = <REVERSE>) {

#       	my $firstChar = substr($lineReverse, 0,1);

		if ($lineReverse =~ /^@/) {

#                if ($firstChar eq "@") {

			chomp($lineReverse);

			my $name = substr($lineReverse, 1, -2);

			if( defined($reads{$name}) ) {

				$reads{$name} = "P";

			} else {

				$reads{$name} = "S";

			} #end else

			++$lineCounter;

		} #end if 

        } #end while

	print "Total number of reverse reads $lineCounter\n";

	print LOG "Total number of reverse reads $lineCounter\n";

	close(REVERSE);


	open(FORWARD, $forward);

	open(FORWARDPAIRED, ">$forward.pairedReads.fq");
	open(FORWARDSINGLE, ">$forward.singleReads.fq");

	print PAIREDFILES "$forward.pairedReads.fq","\n";
	print SINGLEFILES "$forward.singleReads.fq","\n";

	my $forwardPaired = 0;
	my $forwardSingle = 0;

	print "Printing into files for $forward.\n";
	
	while(my $lineForward = <FORWARD>) {

#       	my $firstChar = substr($lineForward, 0,1);

		if ($lineForward =~ /^@/) {

#                if ($firstChar eq "@") {

			chomp($lineForward);

			my $header = substr($lineForward, 1,-2);

#			print $header, "\n";

			if($reads{$header} eq "P" ) {

				++$forwardPaired;

				print FORWARDPAIRED $lineForward, "\n";

				my $sequence = <FORWARD>;
				chomp($sequence);
				print FORWARDPAIRED $sequence, "\n";

				my $secondHeader = <FORWARD>;
				chomp($secondHeader);
				print FORWARDPAIRED $secondHeader, "\n";

				my $quality = <FORWARD>;
				chomp($quality);		
				print FORWARDPAIRED $quality, "\n";

			} else {

				++$forwardSingle;

				print FORWARDSINGLE $lineForward, "\n";

				my $sequence = <FORWARD>;
				chomp($sequence);
				print FORWARDSINGLE $sequence, "\n";

				my $secondHeader = <FORWARD>;
				chomp($secondHeader);
				print FORWARDSINGLE $secondHeader, "\n";

				my $quality = <FORWARD>;
				chomp($quality);		
				print FORWARDSINGLE $quality, "\n";

			} #end else

		} #end if

	} #end while


	open(REVERSE, $reverse);

	open(REVERSEPAIRED, ">$reverse.pairedReads.fq");
	open(REVERSESINGLE, ">$reverse.singleReads.fq");
	
	print PAIREDFILES "$reverse.pairedReads.fq","\n";
	print SINGLEFILES "$reverse.singleReads.fq","\n";

	my $reversePaired = 0;
	my $reverseSingle = 0;

	print "Printing into files for $reverse.\n";

	while(my $lineReverse = <REVERSE>) {

#       	my $firstChar = substr($lineReverse, 0,1);

		if ($lineReverse =~ /^@/) {

#                if ($firstChar eq "@") {

			chomp($lineReverse);

			my $header = substr($lineReverse, 1,-2);

			if($reads{$header} eq "P" ) {

				++$reversePaired;

				print REVERSEPAIRED $lineReverse, "\n";

				my $sequence = <REVERSE>;
				chomp($sequence);
				print REVERSEPAIRED $sequence, "\n";

				my $secondHeader = <REVERSE>;
				chomp($secondHeader);
				print REVERSEPAIRED $secondHeader, "\n";

				my $quality = <REVERSE>;
				chomp($quality);		
				print REVERSEPAIRED $quality, "\n";

			} else {

				++$reverseSingle;

				print REVERSESINGLE $lineReverse, "\n";

				my $sequence = <REVERSE>;
				chomp($sequence);
				print REVERSESINGLE $sequence, "\n";

				my $secondHeader = <REVERSE>;
				chomp($secondHeader);
				print REVERSESINGLE $secondHeader, "\n";

				my $quality = <REVERSE>;
				chomp($quality);		
				print REVERSESINGLE $quality, "\n";

			} #end else

		} #end if

	} #end while

	print LOG "file\tPairedReads\tSingleReads\n";
	print LOG "$forward\t$forwardPaired\t$forwardSingle\n";
	print LOG "$reverse\t$reversePaired\t$reverseSingle\n";
	print LOG "\n\n";

} #end while for input files

exit;


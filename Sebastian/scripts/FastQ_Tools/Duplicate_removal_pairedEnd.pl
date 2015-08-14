#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[2])) {

	print "Please input a pair of fastq files and length of seed\n";

	die "Missing input file\n";

} #end if

my $forward = $ARGV[0];

my $reverse = $ARGV[1];

my $seed = $ARGV[2];

open(FORWARD, $forward) or die "Can't open $forward\n";
open(REVERSE, $reverse) or die "Can't open $reverse\n";

my %data;

print "Processing $forward file.\n";

my $lineCounter = 0;

while(my $lineForward = <FORWARD>) {
    
    chomp($lineForward);

    my $progress = $lineCounter/100000;

    if (int($progress) eq $progress) { print "$lineCounter reads has been analyzed.\n" } #end if

    my $header = substr($lineForward, 0,-2);

    my $sequenceForward = <FORWARD>;
    my $secondHeaderForward = <FORWARD>;
    my $qualForward = <FORWARD>;
    
    my $headerReverse = <REVERSE>;
    my $sequenceReverse = <REVERSE>;
    my $secondHeaderReverse = <REVERSE>;
    my $qualReverse = <REVERSE>;
    
    my $sequence = substr($sequenceForward,0,$seed)."~".substr($sequenceReverse,0,$seed);

    $data{$sequence} = $header;

    ++$lineCounter;

} #end while

print $lineCounter, "\n";

close(FORWARD);
close(REVERSE);

print "Loading unique hash done, selecting unique IDs.\n";

open(OUTPUT, ">output.txt");


my %uniqReads;

foreach my $uniqSequence (keys(%data)) {

    print OUTPUT $data{$uniqSequence}, "\n";

    $uniqReads{$data{$uniqSequence}} = 1;

} #end foreach

print scalar(keys(%data) ), "\n";

open(FORWARD, $forward);
open(REVERSE, $reverse);

open(FORWARDUNIQ, ">$forward.uniqReads.fq");
open(REVERSEUNIQ, ">$reverse.uniqReads.fq");

print "Printing into unique $forward file.\n";
	
while(my $lineForward = <FORWARD>) {

    chomp($lineForward);

    my $header = substr($lineForward, 0,-2);

    my $sequenceForward = <FORWARD>;
    my $secondHeaderForward = <FORWARD>;
    my $qualForward = <FORWARD>;
    
    my $headerReverse = <REVERSE>;
    my $sequenceReverse = <REVERSE>;
    my $secondHeaderReverse = <REVERSE>;
    my $qualReverse = <REVERSE>;
        
    if(defined($uniqReads{$header})) {

	print FORWARDUNIQ $lineForward, "\n";
	print FORWARDUNIQ $sequenceForward;
	print FORWARDUNIQ $secondHeaderForward;
	print FORWARDUNIQ $qualForward;

	print REVERSEUNIQ $headerReverse;
	print REVERSEUNIQ $sequenceReverse;
	print REVERSEUNIQ $secondHeaderReverse;
	print REVERSEUNIQ $qualReverse;

    } #end if

} #end while

close(FORWARD);
close(REVERSE);
close(FORWARDUNIQ);
close(REVERSEUNIQ);

exit;


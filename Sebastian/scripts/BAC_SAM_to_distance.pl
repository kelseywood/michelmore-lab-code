#!/usr/bin/perl
use strict;
use warnings;

#

if(!defined($ARGV[2])) {

  print "Please provide two SAM files (Forward and Reverse) and output file\n";
  die "Missing arguments";
} #end argument check

my %forwReads;
my %refSizes;

open(FORWARD, $ARGV[0]) or die "Can't open forward file $ARGV[0]\n";

while (my $forwAlignment = <FORWARD>) {

  if($forwAlignment =~ /^@/) {

    my @headInfo = split("\t", $forwAlignment);
    $headInfo[2] =~ s/LN://;
    $headInfo[1] =~ s/SN://;

    $refSizes{$headInfo[1]} = $headInfo[2];

    next;

  } #getting the size of the references from the sam header

  my @alignmentData = split("\t", $forwAlignment);

  if($alignmentData[2] eq "*") {next}

  my $start = $alignmentData[3];
  my $end = $start + length($alignmentData[9]);

  my $readName = $alignmentData[0];
  $readName =~ s/\.b1//;

  push(@{$forwReads{$readName}}, [$alignmentData[2], $alignmentData[1], $start, $end]);

} #end reading forward file

my %revReads;

open(REVERSE, $ARGV[1]) or die "Can't open reverse file $ARGV[1]\n";

while (my $revAlignment = <REVERSE>) {

  if($revAlignment =~ /^@/) {next}

  my @alignmentData = split("\t", $revAlignment);

  if($alignmentData[2] eq "*") {next}

  my $start = $alignmentData[3];
  my $end = $start + length($alignmentData[9]);

  my $readName = $alignmentData[0];
  $readName =~ s/\.g1//;

 push(@{$revReads{$readName}}, [$alignmentData[2], $alignmentData[1], $start, $end]);

} #end reading reverse file

open(OUTPUT, ">$ARGV[2]");
open(DATA, ">$ARGV[2].data.txt");

my $maxInsert = 150000;
my $meanInsert = 100000;

foreach my $read (keys %forwReads) {

  if(!defined($revReads{$read})) {next}

  if( scalar(@{$forwReads{$read}}) == 1 && scalar(@{$revReads{$read}}) == 1 ) {

   my @forwRead = @{$forwReads{$read}};
   my @revRead = @{$revReads{$read}};

   print STDERR $forwRead[0][0], "\t", $revRead[0][0], "\n";

   my $distance = "-";
   my $senses = "-";
   my $stat = "blank";

   my $forwRef = $forwRead[0][0];
   my $forwFlag = $forwRead[0][1];
   my $forwStart = $forwRead[0][2];
   my $leftForwTail = $forwRead[0][2];
   my $rightForwTail = $refSizes{$forwRead[0][0]} - $forwRead[0][3];

   my $revRef = $revRead[0][0];
   my $revFlag = $revRead[0][1];
   my $revStart = $revRead[0][2];
   my $leftRevTail = $revRead[0][2];
   my $rightRevTail = $refSizes{$revRead[0][0]} - $revRead[0][3];

   if ( ($leftForwTail > $maxInsert || $rightForwTail > $maxInsert) && !($leftForwTail > $maxInsert && $rightForwTail > $maxInsert) && ($leftRevTail > $maxInsert || $rightRevTail > $maxInsert) && !($leftRevTail > $maxInsert && $rightRevTail > $maxInsert) ) {
#   if( $refSizes{$forwRef} > $maxInsert*2 && $refSizes{$revRef} > $maxInsert*2 ) {

      $stat = "DoubleAnchor";

      if( $forwStart < $maxInsert && $revStart > $maxInsert && $rightRevTail < $maxInsert) {

	if( $forwFlag == 16 && $revFlag == 0) {$senses = "-/-"; $distance = $meanInsert - $leftForwTail - $rightRevTail }

      } elsif( $forwStart < $maxInsert && $revStart < $maxInsert) {
	
	if( $forwFlag == 16 && $revFlag == 16) {$senses = "-/+"; $distance = $meanInsert - $leftForwTail - $leftRevTail }

      } elsif( $forwStart > $maxInsert && $revStart > $maxInsert && $rightForwTail < $maxInsert && $rightRevTail < $maxInsert) {
	
	if( $forwFlag == 0 && $revFlag == 0) {$senses = "+/-"; $distance = $meanInsert - $rightForwTail - $rightRevTail }

      } elsif( $forwStart > $maxInsert && $revStart < $maxInsert && $rightForwTail < $maxInsert) {
	
 	if( $forwFlag == 0 && $revFlag == 16) {$senses = "+/+"; $distance = $meanInsert - $rightForwTail - $leftRevTail }

      } #end check which end are they

    } elsif ( ($leftForwTail > $maxInsert || $rightForwTail > $maxInsert) && !($leftForwTail > $maxInsert && $rightForwTail > $maxInsert) && ($leftRevTail < $maxInsert && $rightRevTail < $maxInsert) ) {
#    } elsif( $refSizes{$forwRef} > $maxInsert*2 && $refSizes{$revRef} < $maxInsert*2 ) {

      $stat = "SingleAnchor";

      if($forwStart < $maxInsert) {

	if($forwFlag == 16) {
	  if($revFlag == 0) { $senses = "-/-"; $distance = $meanInsert - $leftForwTail - $rightRevTail }
	  elsif ($revFlag == 16) { $senses = "-/+"; $distance = $meanInsert - $leftForwTail - $leftRevTail } #end of orienting rev read
	} #end verifying orientation

      } elsif($forwStart > $maxInsert && $rightForwTail < $maxInsert) {

	if($forwFlag == 0) {
	  if($revFlag == 0) { $senses = "+/-"; $distance = $meanInsert - $rightForwTail - $rightRevTail }
	  elsif ($revFlag == 16) { $senses = "+/+"; $distance = $meanInsert - $rightForwTail - $leftRevTail } #end of orienting rev read
	} #end verifying orientation

      } #end anchoring large scaffold

    } elsif ( ( $leftForwTail < $maxInsert && $rightForwTail < $maxInsert) && ($leftRevTail > $maxInsert || $rightRevTail > $maxInsert) && !($leftRevTail > $maxInsert && $rightRevTail > $maxInsert) ) {
#    } elsif( $refSizes{$forwRef} < $maxInsert*2 && $refSizes{$revRef} > $maxInsert*2 ) {

      $stat = "SingleAnchor";

      if($revStart < $maxInsert) {

	if($revFlag == 16) {
	  if($forwFlag == 0) { $senses = "+/+"; $distance = $meanInsert - $rightForwTail - $leftRevTail }
	  elsif ($forwFlag == 16) { $senses = "-/+"; $distance = $meanInsert - $leftForwTail - $leftRevTail } #end of orienting rev read
	} #end verifying orientation

      } elsif($revStart > $maxInsert && $rightRevTail < $maxInsert) {

	if($revFlag == 0) {
	  if($forwFlag == 0) { $senses = "+/-"; $distance = $meanInsert - $rightForwTail - $rightRevTail }
	  elsif ($forwFlag == 16) { $senses = "-/-"; $distance = $meanInsert - $leftForwTail - $rightRevTail } #end of orienting rev read
	} #end verifying orientation

      } #end anchoring large scaffold

    } elsif ( ( $leftForwTail < $maxInsert && $rightForwTail < $maxInsert) && ($leftRevTail < $maxInsert && $rightRevTail < $maxInsert) ) {
#    } elsif( $refSizes{$forwRef} < $maxInsert*2 && $refSizes{$revRef} < $maxInsert*2 ) {

      $stat = "NoAnchor";

      if($forwFlag == 16) {
	if($revFlag == 0) { $senses = "-/-"; $distance = $meanInsert - $leftForwTail - $rightRevTail }
	elsif ($revFlag == 16) { $senses = "-/+"; $distance = $meanInsert - $leftForwTail - $leftRevTail } #end of orienting rev read
      } elsif($forwFlag == 0) {
	if($revFlag == 0) { $senses = "+/-"; $distance = $meanInsert - $rightForwTail - $rightRevTail }
	elsif ($revFlag == 16) { $senses = "+/+"; $distance = $meanInsert - $rightForwTail - $leftRevTail } #end of orienting rev read
      } #end verifying orientation

   } #end elsif

    print OUTPUT $forwRef, "\t", $revRef, "\t", $distance, "\t", $senses, "\t", $stat, "\t", $read, "\n";
    print DATA $read, "\t", $forwRef, "\t", $forwFlag, "\t", $forwStart, "\t", $leftForwTail, "\t", $rightForwTail, "\t", $revRef, "\t", $revFlag, "\t", $revStart, "\t", $leftRevTail, "\t", $rightRevTail, "\n";

  } #end if not multiple mapping and both ends are mapped

}# end analysing each read

close(OUTPUT);

exit;


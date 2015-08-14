#!/usr/bin/perl
use POSIX;
use strict;
use warnings;

if (@ARGV < 2) {

	print "Two are arguments are needed, please inpu them.\n";
        print "Fasta file and maximun sequence size.\n";

        exit 0;

} #end if

my $input = $ARGV[0];
my $maxSize = $ARGV[1];

open(FASTAINPUT, $input);
open(OUTPUT, ">$input.splitted_max$maxSize.txt");

my $sequenceName;
my $sequenceSize = 0;
my $sequence;

while(my $line = <FASTAINPUT>) {

	chomp($line);

	$line =~ s/\r//;

	if( $line =~ /^>/) {

		if ($sequenceSize > $maxSize) {

			my $pieceSize = ceil($sequenceSize/ceil($sequenceSize/$maxSize));

			my $start = 0;

			my $piece = 1;

			while ($start <= $sequenceSize) {

			        my $subSequence = substr($sequence,$start,$pieceSize);

				if($subSequence eq "") {last}

				print OUTPUT $sequenceName, "_", $piece, "\n";
				print OUTPUT $subSequence, "\n";

				$start += $pieceSize;
				++$piece;

			} #end while printing piece

			$sequenceSize = 0;
			$sequence = "";

		} elsif ($sequenceSize > 0) {

			print OUTPUT $sequenceName, "\n";
			print OUTPUT $sequence, "\n";

		} #end if sequence to long

		my @header = split(" ", $line);

		$sequenceName = $header[0];	

	} else {

		$sequenceSize += length($line);
		$sequence .= $line;

	} #end else

} #end while

if ($sequenceSize > $maxSize) {

	my $pieceSize = ceil($sequenceSize/ceil($sequenceSize/$maxSize));

#	print $sequenceSize, "\n";

#	print $pieceSize, "\n";

	my $start = 0;

	my $piece = 1;

	while ($start <= $sequenceSize) {

	        my $subSequence = substr($sequence,$start,$pieceSize);

		print OUTPUT $sequenceName, "_", $piece, "\n";
		print OUTPUT $subSequence, "\n";

		$start += $pieceSize;
		++$piece;

	} #end while printing piece

} elsif ($sequenceSize > 0) {

	print OUTPUT $sequenceName, "\n";
	print OUTPUT $sequence, "\n";

} #end if sequence to long

close(FASTAINPUT);
close(OUTPUT);

exit;


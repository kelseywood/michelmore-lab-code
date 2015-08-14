#!/usr/bin/perl
use strict; use warnings;
use List::Util qw(sum);

if (@ARGV < 2) {

	print "Three arguments are need it please input them.\n";
	print "File to summarize.\n";
	print "Windown size.\n";
	print "Allowed gap size.\n";
	print "Column to summarize (first column equal to 0).\n";
	print "Ussage: SlidingWindow_Sumarizer.pl <File> <WindownSize> <AllowSpace> <Column>.\n";
	die;

} #end if

my $file = $ARGV[0];
my $windowSize = $ARGV[1];
my $allowSpace = $ARGV[2];
my $columnToSumarize = $ARGV[3];
my $currentScaffold = "";
my $startPosition = 0;
my $prevPosition = 0;

open(FILE, $file);
open(OUTPUT, ">$file.Window_$windowSize.AllowSpace_$allowSpace.txt");

my @data = [];

my $baseLine = <FILE>;
my @line = split("\t", $baseLine);
push(@data, $line[$columnToSumarize]);
$currentScaffold = $line[0];
$startPosition = $line[1];
$prevPosition = $line[1];

while($baseLine = <FILE>) {

	my @line = split("\t", $baseLine);

	if ($currentScaffold eq $line[0]) {

		push(@data, $line[$columnToSumarize]);

		if (((scalar @data) >= $windowSize) || (abs($line[1] - $prevPosition) > $allowSpace)) {

			shift(@data);

			my $num = scalar(@data);
			my $sum = sum(0, @data);
			my $average = $sum/(scalar @data);
			my $sqsum = 0;
			foreach my $x (@data) { $sqsum += ($average - $x) ** 2 }
			my $var = $sqsum/$num;
			my $stand_dev = sqrt ($var);

			print OUTPUT $currentScaffold, "\t", $startPosition, "\t", $prevPosition, "\t", sprintf("%.2f", $average), "\t", sprintf("%.2f", $stand_dev), "\n";

			@data = [];

			$startPosition = $line[1];

		} #end if

		$prevPosition = $line[1];

	} else {


		shift(@data);

		my $num = scalar(@data);
		my $sum = sum(0, @data);
		my $average = $sum/(scalar @data);
		my $sqsum = 0;
		foreach my $x (@data) { $sqsum += ($average - $x) ** 2 }
		my $var = $sqsum/$num;
		my $stand_dev = sqrt ($var);

		print OUTPUT $currentScaffold, "\t", $startPosition, "\t", $prevPosition, "\t", sprintf("%.2f", $average), "\t", sprintf("%.2f", $stand_dev), "\n";

		@data = [];

		$startPosition = $line[1];

		push(@data, $line[$columnToSumarize]); 

		$currentScaffold = $line[0];

		$startPosition = $line[1];

		$prevPosition = $line[1];

	} #end else

} #end while

close(FILE);
close(OUTPUT);

exit;



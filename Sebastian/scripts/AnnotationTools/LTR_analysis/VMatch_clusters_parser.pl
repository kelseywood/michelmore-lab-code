#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[1])) {
	print "Missing inputs, please provide an output file and one VMatch output file\n";
	die
} #end else

my $output = shift(@ARGV);

my %sequences;
my @files;

foreach my $file (@ARGV) {

	my $fileBase = $file;

	$fileBase =~ s{.*/}{};      # removes path  
	$fileBase =~ s{\.[^.]+$}{}; # removes extension

	open(INPUT, $file) or print STDERR "Can't open $file\n";

	my $cluster;

	push(@files, $fileBase);

	while(my $line = <INPUT>) {

		if($line =~ /^#/) {next}

		elsif($line =~ /^\d/) {

			chomp($line);

			$line =~ s/://;

			$cluster = $line;

		} else {

			chomp($line);

			$line =~ s/\s//g;

			$sequences{$line}{$fileBase} = $cluster;

		} #end else

	} #end while reading clusters file

	close(INPUT);

} #end of looping through files

print "Finish loading cluster data\n\n";

open(OUTPUT, ">$output");

print OUTPUT "Sequence\t", join("\t", @files), "\n";

foreach my $sequence (sort keys (%sequences)) {

	print OUTPUT $sequence;

	foreach my $file (@files) {

		print OUTPUT "\t";

		if( defined($sequences{$sequence}{$file}) ) { print OUTPUT $sequences{$sequence}{$file} }
		else { print OUTPUT "NA"}

	} #end foreach $file

	print OUTPUT "\n";

} #end printing

print "Done printing done, all good\n";

close(OUTPUT);

exit;



				

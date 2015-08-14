#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 3) {

	print "Usage: Input text file with IDs to be retain, GFF file where to look for the data and output gff file name.\n";
	print "GFF_filterer_byType.pl IDsList.txt In.gff Out.gff\n";
	die "Missing Arguments";

} #end if

open(IDSFILE, $ARGV[0]) or die "Can't open file $ARGV[0]\n";;

my %IDs;

while( my $id = <IDSFILE>) {

	chomp($id);

	$IDs{$id} = 0;

} #end foreach

my $in = $ARGV[1];
my $out = $ARGV[2];

open(GFF, $in) or die "Can't open $in\n";
open(OUTPUT, ">$out") or die "Can't open $out\n";

while(my $gffLine = <GFF>) {

	chomp($gffLine);

	my $firstChar = substr($gffLine,0,1);

	my @fields = split("\t", $gffLine);

	my $toprint = "no";

	if($firstChar ne "#" && $firstChar ne "" && defined($fields[2])) {
	
		my @features = split (";", $fields[8]);

		foreach my $feature (@features) {

			my @value = split("=", $feature);

		        if($value[0] eq "ID") {
				if(defined($IDs{$value[1]})  ) {
					$toprint = "yes";
					++$IDs{$value[1]};
				} #end if to print

		        } elsif ($value[0] eq "Name") {
				if(defined($IDs{$value[1]}) ) {
					$toprint = "yes";
					++$IDs{$value[1]};
				} #end if to print

		        } elsif ($value[0] eq "Parent") {
				if(defined($IDs{$value[1]}) ) {
					$toprint = "yes";
					++$IDs{$value[1]};
				} #end if to print

		        } #end elsif

		} #end foreach

		if ($toprint eq "yes") {

			print OUTPUT $gffLine, "\n";

		} #end else

	} #end if

} #end while

print STDERR "Missing Id's\n\n";

foreach my $key ( keys %IDs) { if($IDs{$key} == 0 ) { print STDERR $key, "\n" } }

close(GFF);
close(OUTPUT);


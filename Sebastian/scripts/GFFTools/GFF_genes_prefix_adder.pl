#!/usr/bin/perl
use warnings;
use strict;

###################################################
#
#
#	GFF3 ID/Name/Parent Changer
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "Two are arguments are needed, please inpu them.\n";
	print "GFF File were to add the prefix.\n";
        print "Prefix to add to all IDs, names and parents (do not need a separator, its added automatically on the script).\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <prefix>\n";
        exit 0;

} #end if

my $input = $ARGV[0];
my $prefix = $ARGV[1];

open(GFF, $input);

$input =~ s{.*/}{};
$input =~ s{\.[^.]+$}{};

open(RESULTS, ">$input.withPrefix.gff");


while (my $line = <GFF>) {

	if($line =~ /^#/) { print RESULTS $line }

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	my $ID;
	my $name;
	my $alias;
	my $parent;

	my @otherFeatures;

	print RESULTS $GFF_line[0], "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t";

	foreach my $feature (@features) {

		my @value = split("=", $feature);

		print RESULTS $value[0];

	        if($value[0] eq "ID") {

			print RESULTS "=", $prefix, "-", $value[1], ";";

	        } elsif ($value[0] eq "Name") {

			print RESULTS "=", $prefix, "-", $value[1], ";";

	        } elsif ($value[0] eq "Parent") {

			$value[1] =~ s/,/,$prefix-/g;

			print RESULTS "=", $prefix, "-", $value[1], ";";

	        } else {

			print RESULTS "=", $value[1], ";";

	        } #end else

	} #end foreach

	print RESULTS "\n";
	
} #end while

close(GFF);
close(RESULTS);

exit;

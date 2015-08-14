#!/usr/bin/perl
use strict; use warnings;
###################################################
#
#
#	GFF3 Orphan feature remover
#
#
#
###################################################

#

if (@ARGV < 1) {

	print "One argument is needed, please input it.\n";
	print "GFF File were to remove orphan features.\n";
	print "perl GFF_orphan_remover.pl <GFFFile>\n";
        exit 0;

} #end if

my %ids;

open(GFF, $ARGV[0]);
open(RESULTS, ">$ARGV[0].withoutOrphans.gff");
open(ORPHANS, ">$ARGV[0].Orphans.gff");

while (my $line = <GFF>) {

	chomp($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	my $orphan = "FALSE";

	foreach my $feature (@features) {

		my @value = split("=", $feature);

	        if($value[0] eq "ID") {
			$ids{$value[1]} = 1;
			if(defined($ids{$value[1]})) {last}

	        } elsif ($value[0] eq "Parent") {

			if(defined($value[1])) {

				if(defined($ids{$value[1]})) {

					$orphan = "FALSE";

				} else {

					$orphan = "TRUE";

				} #end else

			} else {

				$orphan = "TRUE";

			} #end else

	        } #end elsif

	} #end foreach

	if( $orphan eq "FALSE" ) {
		print RESULTS $GFF_line[0], "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";

	} else {

		print ORPHANS $GFF_line[0], "\t", $GFF_line[1], "\t", $GFF_line[2], "\t", $GFF_line[3], "\t", $GFF_line[4], "\t", $GFF_line[5], "\t", $GFF_line[6], "\t", $GFF_line[7], "\t", $GFF_line[8], "\n";

	} #end else printing

} #end while

close(GFF);
close(RESULTS);

exit;

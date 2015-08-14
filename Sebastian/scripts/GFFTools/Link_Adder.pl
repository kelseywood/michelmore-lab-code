#!/usr/bin/perl
use warnings;

#Open current GFF file and temporal output GFF				
open(GFF, $ARGV[0]) or print "Can't open file $ARGV[1]";

while($GFFLine = <GFF>) {

	chomp($GFFLine);

	@fields = split("\t", $GFFLine);

	if (($fields[1] eq "blastx" && $fields[2] eq "protein_match") || ($fields[1] eq "protein2genome" && $fields[2] eq "protein_match") || ($fields[1] eq "est2genome" && $fields[2] eq "expressed_sequence_match") || ($fields[1] =~ m/^GeneWise/ && $fields[2] eq "mRNA")) {

		@attributes = split(";", $fields[8]);

		foreach my $attribute (@attributes) {

			$type = substr($attribute,0,4);

			if($type eq "Name") {

				$name = substr($attribute, 5, 100);

				$start = substr($name,0,2);

				my $link;

				if($start eq "gi") {

					@parts = split("\\|", $name);

					$link = $parts[1];

				} elsif ($start eq "AT" || $start eq "Gl") {

					@parts = split("\\.", $name);

					$link = $parts[0];


				} elsif ($start eq "GS" || $start eq "PG" || $start eq "Tc") {

					@parts = split("-", $name);

					$link = $parts[0];

				} elsif ($start eq "La") {

					$link = $name;

				} else {

					$link = "";

				} #end if

				print $GFFLine, "Note=", $link, ";\n";

			} #end if

		} #end foreach

	} else {

		print $GFFLine, "\n";

	} #end else

} #end while

close(GFF);

exit;

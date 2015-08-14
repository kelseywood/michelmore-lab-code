#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 2) {

	print "Usage: GFF file where to look for the data and output text file name.\n";
	print "GFF_interproFormat2domainTable.pl In.gff Out.txt\n";
	die "Missing Arguments";

} #end if

my $in = $ARGV[0];
my $out = $ARGV[1];

open(GFF, $in) or die "Can't open $in\n";
open(OUTPUT, ">$out") or die "Can't open $out\n";

while(my $gffLine = <GFF>) {

	chomp($gffLine);

        if($gffLine =~ m/>/) {last}

	if( !($gffLine =~ m/#/) ) {

		my @fields = split("\t", $gffLine);

            	my $gene = $fields[0];
                my $source = $fields[1];
                my $start = $fields[3];
                my $end = $fields[4];
                my $score = $fields[5];
                my $strand = $fields[6];
                my $domain = "";
                my $signature_desc = "";
                my $date = "";
                my $Ontology_term = "";
                my $Dbxref = "";

                if($source eq '.') {next}

		my @features = split (";", $fields[8]);

		my  $name = "";

		foreach my $feature (@features) {

			my @value = split("=", $feature);

			if($value[0] eq "Name") {

				$domain = $value[1];

			} elsif($value[0] eq "signature_desc") {

				$signature_desc = $value[1];

			} elsif($value[0] eq "date") {

				$date = $value[1];

			} elsif($value[0] eq "Ontology_term") {

				$Ontology_term = $value[1];

			} elsif($value[0] eq "Dbxref") {

				$Dbxref = $value[1];

		        } #end if require value

		} #end foreach

                my @geneModels = split("\\|", $gene);

                foreach my $geneModel (@geneModels) {

	                print OUTPUT $geneModel, "\t",
	                        $domain, "\t",
	                        $source, "\t",
	                        $start, "\t",
	                        $end, "\t",
	                        $score, "\t",
	                        $strand, "\t",
	                        $signature_desc, "\t",
	                        $date, "\t",
	                        $Dbxref, "\n";

                } #end of gene models

	} #end if

} #end while

close(GFF);
close(OUTPUT);

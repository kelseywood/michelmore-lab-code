#!/usr/bin/perl

######################################################
#                                                    #
#          Maker qsub annotation sumarizer           #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to gather all the results from independent maker runs into single gff file per source
# Part of the main pipeline "Maker_qsub_annotator.pl"

use strict;
use warnings;

my $gffprefix = $ARGV[0];
my $gffDir = $ARGV[1];

##### Generate complete GFF #####

print STDERR "Consolidating annotation info by source\n";

if(-e "$gffprefix*.gff") {

	system("rm $gffprefix*.gff");

} #end if file exist delete then

opendir(GFFFILES, $gffDir);

while(my $gff = readdir(GFFFILES) ) {

	if($gff ne "." && $gff ne "..") {

		if($gff =~ /.gff$/) {

			my %annotation;

			open(GFF, "GFF_files/$gff");

			while(my $gffLine = <GFF>) {

				if($gffLine =~ /^[[:alnum:]]/) {

					my @fields = split("\t", $gffLine);

				        if($fields[1] =~ /^[[:alnum:]]/) {

				        	push(@{$annotation{$fields[1]}}, $gffLine);

				        } #end if source starts with a letter

				} elsif ($gffLine =~ /^>/) {

					last;

				} #end elsif line

			} #end while

			close(GFF);

			foreach my $source (keys(%annotation)) {

				open(SOURCE, ">>$gffprefix.maker.output.$source.gff");

				foreach my $line (@{$annotation{$source}}) { print SOURCE $line }

				close(SOURCE);

			} #end foreach

		} #end if its a gff

	} #end if its a real file

} #end while gff files

exit;

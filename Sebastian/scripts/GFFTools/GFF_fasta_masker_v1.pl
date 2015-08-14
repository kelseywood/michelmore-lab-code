#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Bio::Seq;
use Bio::Index::Fasta;

if (@ARGV < 4) {

	print "Usage: print input GFF file for what to mask, fasta input file, fasta output file, masking type.\n";
	print "GFF_masker.pl <regions.gff> <IN.fasta> <OUT.fasta> <HARD/SOFT>\n";

	die;

} #end if

my $maskingType = $ARGV[3];

if ( ($maskingType ne "HARD") && ($maskingType ne "SOFT") ) {

	print "Please specify HARD or SOFT masking with those specific words.\n";
        die "Non recognize masking type\n";

} #end masking type wrong

my %maskingRegions;

open(GFF, $ARGV[0]);

while(my $gffLine = <GFF>) {

	chomp($gffLine);

	my @fields = split("\t", $gffLine);

	$maskingRegions{$fields[0]}{$fields[3]} = $fields[4];

} #end while

foreach my $scaffold ( keys %maskingRegions) {

	foreach my $start (sort {$a<=>$b} keys %{$maskingRegions{$scaffold}}) {

		my $end = $maskingRegions{$scaffold}{$start};

		foreach my $otherStart (sort {$a<=>$b} keys %{$maskingRegions{$scaffold}}) {
		
			if( ( $end > $otherStart ) && ( $maskingRegions{$scaffold}{$otherStart} > $start ) ) {

				$maskingRegions{$scaffold}{$start} = $maskingRegions{$scaffold}{$otherStart};
				delete($maskingRegions{$scaffold}{$otherStart});

			} #end if overlap
		} #end second foreach
	} #end first foreach
} #end foreach

#foreach my $scaffold ( keys %maskingRegions) {
#	foreach my $start (sort {$a<=>$b} keys %{$maskingRegions{$scaffold}}) {
#	print $scaffold, "\t", $start, "\t", $maskingRegions{$scaffold}{$start}, "\n";
#	} #end first foreach
#} #end foreach



close(GFF);

my $idx = Bio::Index::Fasta->new(
                                 -filename => "$ARGV[1].idx",
                                 -write_flag => 1
                                );
$idx->make_index("$ARGV[1]");

open(FASTA, $ARGV[1]) or die "can't open fasta input\n";;

open(NEWFASTA, ">$ARGV[2]") or die "can't open fasta output\n";

while(my $fastaLine = <FASTA>) {

	chomp($fastaLine);

        my $firstChar = substr($fastaLine,0,1);

        if($firstChar eq ">") {

                my @header = split(" ", $fastaLine);

                my $seqName = substr($header[0], 1, 60);

                my @name = split("_", $seqName);

		my $seq = $idx->fetch($seqName);

                my $length = $seq->length;

#		print "Cleaning $seqName ";

                if(defined($maskingRegions{$seqName})) {

                	print NEWFASTA ">", $seqName, " ", $length, "\n";

                        my $prevEnd = 0;

                        foreach my $region (sort {$a<=>$b} keys %{$maskingRegions{$seqName}}) {

#				print "masking $region..$maskingRegions{$seqName}{$region}\n";
#				print "masking $region..$prevEnd\n";

				if ( $region > $prevEnd ) {

					if($region ne 1) {

	                                	my $sequence1 = $seq->subseq(($prevEnd+1), ($region-1));

        	                                print NEWFASTA $sequence1;

					} #end if

                                } else {

					$region = $prevEnd + 1;

				} #ends

                                my $sequence2 = $seq->subseq($region, $maskingRegions{$seqName}{$region});

                                if($maskingType eq "HARD") {

	                                $sequence2 =~ tr/ACGT/NNNN/;

                                } else {

	                                $sequence2 =~ tr/ACGT/acgt/;

                                } #end masking

                                print NEWFASTA $sequence2;

                                $prevEnd = $maskingRegions{$seqName}{$region};

                        } #end foreach masking regions

                        if ($prevEnd ne $length) {

                        	my $sequence3 = $seq->subseq(($prevEnd+1), $length);

                                print NEWFASTA $sequence3;

                        } #end if print tail

#			print "\n";

                } else {

                	print NEWFASTA ">", $seqName, " ", $length, "\n";

                        my $sequence = $seq->subseq(1, $length);

                        print NEWFASTA $sequence;

                } #end else

                print NEWFASTA "\n";

        } #end if sequence header

} #end while fasta file

close (FASTA);
close (NEWFASTA);

exit;

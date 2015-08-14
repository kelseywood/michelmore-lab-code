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

                my $version = substr($name[2], 1, 1);

                ++$version;

                my $newName = $name[0]."_".$name[1]."_v".$version."_".$name[3]."_".$name[4]."_".$name[5];

		my $seq = $idx->fetch($seqName);

                my $length = $seq->length;

                if(defined($maskingRegions{$seqName})) {

                	print NEWFASTA ">", $newName, " ", $length, "\n";

                        my $prevEnd = 0;

                        foreach my $region (sort {$a<=>$b} keys %{$maskingRegions{$seqName}}) {

                        	if($region ne 1) {

                                	my $sequence1 = $seq->subseq(($prevEnd+1), ($region-1));

                                        print NEWFASTA $sequence1;

                                } #end if

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

                } else {

                	print NEWFASTA ">", $newName, " ", $length, "\n";

                        my $sequence = $seq->subseq(1, $length);

                        print NEWFASTA $sequence;

                } #end else

                print NEWFASTA "\n";

        } #end if sequence header

} #end while fasta file

close (FASTA);
close (NEWFASTA);

exit;

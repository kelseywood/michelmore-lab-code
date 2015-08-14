#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;

###################################################
#
#
#	GFF3 Genomic Sequence Extractor
#
#
#
###################################################

#

if (@ARGV < 2) {

	print "Three are arguments are needed, please inpu them.\n";
        print "Fasta reference file where to take the sequence from.\n";
	print "GFF File were to take position information.\n";
	print "Size of the flanking region out of the ORF.\n";
	print "Extract intron, Yes o No.\n\n";
	print "perl GFF_mRNA_GenomicSequence_Extractor.pl <FastaFile> <GFFFile> <Flanking size> <Yes or No>\n";
        exit 0;

} #end if

my $idx = Bio::Index::Fasta->new(
                                 -filename => "$ARGV[0].idx",
                                 -write_flag => 1
                                );
$idx->make_index("$ARGV[0]");

my $flanking_span = $ARGV[2];

my %mRNAs;

my $countmRNAs =1;

my $getIntron = $ARGV[3];

open(GFF, $ARGV[1]);

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	foreach my $feature (@features) {

#		print $feature, "\n";

		my @value = split("=", $feature);

	        if ( ($value[0] eq "Name") && ($GFF_line[2] eq "mRNA") ) {

			$mRNAs{$value[1]}{'data'} = [$GFF_line[0],$GFF_line[1],$GFF_line[6],$GFF_line[8]];

			++$countmRNAs;

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "CDS") ) {

			$mRNAs{$value[1]}{$GFF_line[3]} = $GFF_line[4];

	        } #end else

	} #end foreach

} #end while

print "$countmRNAs were loaded and will be use for sequence extraction\n\n";

close(GFF);

if($getIntron eq "Yes") {

	open(FASTA, ">$ARGV[1].flank$flanking_span.withIntrons.fasta");

} elsif ($getIntron eq "No") {

	open(FASTA, ">$ARGV[1].flank$flanking_span.withoutIntrons.fasta");

} else {

	die "Please specifi if you want introns or not with 'Yes' or 'No'\n";

} #end else


foreach my $mRNA (keys %mRNAs) {

#	print "Analyzing $mRNA\n";

	my @mRNA_info = @{$mRNAs{$mRNA}{'data'}};

	delete $mRNAs{$mRNA}{'data'};

	print FASTA ">$mRNA $mRNA_info[1]";

	if ($mRNA_info[2] eq "+") {

		print FASTA " (+ strand)\n";

		my $firstExon = "true";

		my $prevEnd = 0;

		foreach my $cds (sort {$a<=>$b} keys %{$mRNAs{$mRNA}} ) {

			if ( ($firstExon eq "true") && ($cds > 1) ) {

				my $seq = $idx->fetch($mRNA_info[0]);

				my $start = $cds-$flanking_span;

				if($start < 1 ) {

					$start = 1;

				} #else is start below
			
				my $end = $cds-1;

				my $left_flank = lc($seq->subseq($start, $end));

				print FASTA $left_flank;

			} elsif ($firstExon eq "false" && ($getIntron eq "Yes") ) {

				my $start = $prevEnd + 1;
			
				my $end = $cds - 1;

				my $seq = $idx->fetch($mRNA_info[0]);
				my $intron = lc($seq->subseq($start, $end));
				
				print FASTA $intron;

			} #end print flank and intron

			my $start = $cds;
			
			my $end = $mRNAs{$mRNA}{$cds};

			my $seq = $idx->fetch($mRNA_info[0]);
			my $exon = uc($seq->subseq($start, $end));
				
			print FASTA $exon;

			$firstExon = "false";

			$prevEnd = $mRNAs{$mRNA}{$cds};

		} #end print CDS

		my $seq = $idx->fetch($mRNA_info[0]);

		if($prevEnd < $seq->length) {

			my $start = $prevEnd+1;

			my $end = $prevEnd+$flanking_span;

			if($end > $seq->length ) {

				$end = $seq->length;

			} #end this

			my $right_flank = lc($seq->subseq($start, $end));

			print FASTA $right_flank;

		} #end if end of sequence

		print FASTA "\n";
		
	} else {

		print FASTA " (- strand)\n";

		my $firstExon = "true";

		my $prevEnd = 0;

		my $prevStart = 0;

		foreach my $cds (sort { $mRNAs{$mRNA}{$b} <=> $mRNAs{$mRNA}{$a} } keys %{$mRNAs{$mRNA}} ) {

			if ($firstExon eq "true") {

				my $seq = $idx->fetch($mRNA_info[0]);

				if($mRNAs{$mRNA}{$cds} < $seq->length) {

					my $end = $mRNAs{$mRNA}{$cds}+$flanking_span;

					if($end > $seq->length ) {

						$end = $seq->length;

					} #end this
			
					my $start = $mRNAs{$mRNA}{$cds}+1;

					my $left_flank = reverse(lc($seq->subseq($start, $end)));

					$left_flank =~ tr/acgt/tgca/;

					print FASTA $left_flank;

				} #end if end of sequence

			} elsif ( ($firstExon eq "false") && ($getIntron eq "Yes") ) {

				my $end = $prevStart - 1;
			
				my $start = $mRNAs{$mRNA}{$cds} + 1;

				my $seq = $idx->fetch($mRNA_info[0]);
				my $intron = reverse(lc($seq->subseq($start, $end)));

				$intron =~ tr/acgt/tgca/;
				
				print FASTA $intron;

			} #end print flank and intron

			my $start = $cds;
			
			my $end = $mRNAs{$mRNA}{$cds};

			my $seq = $idx->fetch($mRNA_info[0]);
			my $exon = reverse(uc($seq->subseq($start, $end)));

			$exon =~ tr/ACGT/TGCA/;
				
			print FASTA $exon;

			$firstExon = "false";

			$prevStart = $cds;

			$prevEnd = $mRNAs{$mRNA}{$cds};

		} #end print CDS

		if ($prevStart > 1) {

			my $seq = $idx->fetch($mRNA_info[0]);

			my $start = $prevStart-$flanking_span;

			my $end = $prevStart-1;

			if ($start < 1) {

				$start = 1;

			} #end if end after

			my $right_flank = reverse(lc($seq->subseq($start, $end)));

			$right_flank =~ tr/acgt/tgca/;

			print FASTA $right_flank;

		} #end if

		print FASTA "\n";
		
	} #end print plus and minus strands

} #end printing mRNAs

close(FASTA);

exit;

#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

my $CodonTable = Bio::Tools::CodonTable->new();

###################################################
#
#
#	GFF3 Genomic Sequence Extractor
#
#
#
###################################################

#####################################
#
# Diff from v1
# add flaking reginos
#

if (@ARGV < 2) {

	print "Two are arguments are needed, please input them.\n";
        print "Fasta reference file where to take the sequence from.\n";
	print "GFF File were to take position information.\n";
	print "perl GFF_mRNA_GenomicSequence_Extractor.pl <FastaFile> <GFFFile>\n";
        exit 0;

} #end if

#my $idx;

#if(!(-f "$ARGV[0].idx") ) {

	print "Building fasta index\n\n";

	my $idx = Bio::Index::Fasta->new(
		                         '-filename' => "$ARGV[0].idx",
		                         '-write_flag' => 1
		                        );

#} #end if index not present

$idx->make_index("$ARGV[0]");

my %mRNAs;

open(GFF, $ARGV[1]);

while (my $line = <GFF>) {

	if ($line =~ "^#") {next}

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	foreach my $feature (@features) {

		my @value = split("=", $feature);

	        if ( ($value[0] eq "ID") && ($GFF_line[2] eq "mRNA") ) { $mRNAs{$value[1]}{'data'} = [$GFF_line[0],$GFF_line[1],$GFF_line[6],$GFF_line[8],$GFF_line[3],$GFF_line[4]] }

	        elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "exon") ) { $mRNAs{$value[1]}{'exon'}{$GFF_line[3]} = $GFF_line[4] }

	} #end foreach

} #end while

my $countmRNAs = scalar(keys %mRNAs);

print "$countmRNAs mRNA's were loaded and will be use for sequence extraction\n\n";

close(GFF);

open(TRANSCRIPT, ">$ARGV[1].transcript.fasta");
open(ERRORS, ">$ARGV[1].errors.txt");

my $analyzedmRNA = 0;

foreach my $mRNA (sort {$a <=> $b} keys %mRNAs) {

	if( ($analyzedmRNA/500) =~ m/^\d+$/) { print "Done with $analyzedmRNA mRNA's of $countmRNAs\n" } #end if integer

	++$analyzedmRNA;

	if(!defined(@{$mRNAs{$mRNA}{'data'}})) {print ERRORS "Empty dataset for $mRNA\n"; next}

	# Get mRNA info and seq location
	my @mRNA_info = @{$mRNAs{$mRNA}{'data'}};
	delete $mRNAs{$mRNA}{'data'};
	my $seq = $idx->fetch($mRNA_info[0]);

        if(defined($seq)) {

		# Start sequence storage arrays
	        my @transcript;

		# Get sorted start coordinates for CDS's and exons
		my @ExonStarts = sort {$a<=>$b} keys %{$mRNAs{$mRNA}{'exon'}} ;

		# For monoexonic genes
		if(scalar(@ExonStarts) == 1) {

			#Get exon
			my $exonSeq = uc($seq->subseq($ExonStarts[0], $mRNAs{$mRNA}{'exon'}{$ExonStarts[0]}));
			push(@transcript, $exonSeq);

		# For non-monoexonic genes
		} else {

#			print STDERR $mRNA, "\n";

			# Get first CDS
			my $firstExon = shift(@ExonStarts);
			my $firstExonSeq = uc($seq->subseq($firstExon, $mRNAs{$mRNA}{'exon'}{$firstExon}));
			push(@transcript, $firstExonSeq);

			# Get information out of exons
		        foreach my $exon (@ExonStarts) {

			        my $exonSeq = uc($seq->subseq($exon, $mRNAs{$mRNA}{'exon'}{$exon}));

			        #Print exon
			        push(@transcript, $exonSeq);

		        } #end print exon

		} #end else for monoexonic gene models

		# Print sequences into respective files
		my $transcriptSequence;

                ###  sequence is on plus ###
	        if ($mRNA_info[2] eq "+") { $transcriptSequence = join("",@transcript)

                ###  sequence is on minus ###
	        } else {

			$transcriptSequence = reverse(join("",@transcript));
			$transcriptSequence =~ tr/ACGTacgt/TGCAtgca/;

	        } #end print plus and minus strands

	        print TRANSCRIPT ">$mRNA transcript ".length($transcriptSequence)."bp $mRNA_info[0]:$mRNA_info[4]..$mRNA_info[5] ($mRNA_info[2] strand)\n".$transcriptSequence."\n";

        #sequence doesn't found on index
        } else {

        	print ERRORS "Sequence for $mRNA_info[0] not found in fasta file for mRNA $mRNA\n";

        } #end else

} #end printing mRNAs

exit;







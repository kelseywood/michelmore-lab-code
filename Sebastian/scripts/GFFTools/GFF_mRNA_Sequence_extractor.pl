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

#

if (@ARGV < 1) {

	print "Two are arguments are needed, please input them.\n";
        print "Fasta reference file where to take the sequence from.\n";
	print "GFF File were to take position information.\n";
	print "perl GFF_mRNA_GenomicSequence_Extractor.pl <FastaFile> <GFFFile>\n";
        exit 0;

} #end if

my $idx = Bio::Index::Fasta->new(
                                 '-filename' => "$ARGV[0].idx",
                                 '-write_flag' => 1
                                );
$idx->make_index("$ARGV[0]");

my %mRNAs;

my $countmRNAs = 0;

open(GFF, $ARGV[1]);

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	foreach my $feature (@features) {

		my @value = split("=", $feature);

	        if ( ($value[0] eq "Name") && ($GFF_line[2] eq "mRNA") ) {

			$mRNAs{$value[1]}{'data'} = [$GFF_line[0],$GFF_line[1],$GFF_line[6],$GFF_line[8]];

			++$countmRNAs;

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "UTR_5") ) {

			$mRNAs{$value[1]}{'UTR_5'}{$GFF_line[3]} = $GFF_line[4];

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "UTR_3") ) {

			$mRNAs{$value[1]}{'UTR_3'}{$GFF_line[3]} = $GFF_line[4];

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "CDS") ) {

			$mRNAs{$value[1]}{'CDS'}{$GFF_line[3]} = $GFF_line[4];

	        } #end else

	} #end foreach

} #end while

print "$countmRNAs mRNA's were loaded and will be use for sequence extraction\n\n";

close(GFF);


open(GENOMIC, ">$ARGV[1].genomic.fasta");
open(TRANSCRIPT, ">$ARGV[1].transcript.fasta");
open(CDS, ">$ARGV[1].CDS.fasta");
open(PROTEIN, ">$ARGV[1].protein.fasta");
open(ERRORS, ">$ARGV[1].errors.txt");

#open(GENOMIC2, ">$ARGV[1].genomic2.fasta");
#open(TRANSCRIPT2, ">$ARGV[1].transcript2.fasta");
#open(CDS2, ">$ARGV[1].CDS2.fasta");
#open(PROTEIN2, ">$ARGV[1].protein2.fasta");

my $analyzedmRNA = 0;

foreach my $mRNA (keys %mRNAs) {

	if( ($analyzedmRNA/500) =~ m/^\d+$/) {

		print "Done with $analyzedmRNA mRNA's of $countmRNAs\n";

	} #end if integer

	++$analyzedmRNA;

	if(defined(@{$mRNAs{$mRNA}{'data'}})) {

		my @mRNA_info = @{$mRNAs{$mRNA}{'data'}};

		delete $mRNAs{$mRNA}{'data'};

		my $seq = $idx->fetch($mRNA_info[0]);

		if(defined($seq)) {

	#	        print GENOMIC2 ">", $mRNA, " genomic\n";
	#	        print TRANSCRIPT2 ">", $mRNA, " transcript\n";
	#	        print CDS2 ">", $mRNA, " CDS\n";
	#	        print PROTEIN2 ">", $mRNA, " peptide\n";

			my @genomic;
			my @transcript;
			my @CDS;
			my @protein;

			my $remaining_genomic = 100;
			my $remaining_transcript = 100;
			my $remaining_cds = 100;
			my $remaining_protein = 0;

			my $size_genomic = 0;
			my $size_transcript = 0;
			my $size_cds = 0;
			my $size_protein = 0;

		        #sequence is on plus
			if ($mRNA_info[2] eq "+") {

			        my $firstExon = "true";

			        my $prevEnd = 0;

			        my $codon = "";

			        # print 5' UTR

			        if(defined($mRNAs{$mRNA}{'UTR_5'}) ) {

			                foreach my $utr (sort { $mRNAs{$mRNA}{'UTR_5'}{$b} <=> $mRNAs{$mRNA}{'UTR_5'}{$a} } keys %{$mRNAs{$mRNA}{'UTR_5'}} ) {

			                        my $start = $utr;

			                        my $end = $mRNAs{$mRNA}{'UTR_5'}{$utr};
			                        
			                        my $utrSeq = lc($seq->subseq($start, $end));

			                        push(@genomic, $utrSeq);
			                        push(@transcript, $utrSeq);

	#	                                print GENOMIC2 $utrSeq;
	#	                                print TRANSCRIPT2 $utrSeq;

			                } #end utr foreach

			        } #end printing 5' UTR

				my $codingSequence = "";

			        foreach my $cds (sort {$a<=>$b} keys %{$mRNAs{$mRNA}{'CDS'}} ) {

			                #Print intron

			                if ($firstExon eq "false" ) {

			                        my $start = $prevEnd + 1;

			                        my $end = $cds - 1;

			                        my $intron = lc($seq->subseq($start, $end));

			                        my $start_fragment = 0;

			                        push(@genomic, $intron);

	#	                                print GENOMIC2 $intron;

			                } #end print flank and intron

			                my $start = $cds;

			                my $end = $mRNAs{$mRNA}{'CDS'}{$cds};

			                my $exon = uc($seq->subseq($start, $end));

			                #Print exon

			                push(@genomic, $exon);
			                push(@transcript, $exon);
			                push(@CDS, $exon);

	#	                        print GENOMIC2 $exon;
	#	                        print TRANSCRIPT2 $exon;
	#	                        print CDS2 $exon;

			                $codingSequence = $codingSequence.$exon;

			                $firstExon = "false";

			                $prevEnd = $mRNAs{$mRNA}{'CDS'}{$cds};

			        } #end print CDS

		                #print protein

		                my $codonStart = 0;

		                $codon = substr($codingSequence, $codonStart, 3);
		                my $aa = $CodonTable->translate_strict($codon);

		                while ($aa ne "") {

		                        push(@protein, $aa);

	#                                print PROTEIN2 $aa;

		                        ++$remaining_protein;
		                        ++$size_protein;

		                        $codonStart += 3;

		                        $codon = substr($codingSequence, $codonStart, 3);

		                        $aa = $CodonTable->translate_strict($codon);

		                } #end while


			        if(defined($mRNAs{$mRNA}{'UTR_3'}) ) {

			                foreach my $utr (sort { $mRNAs{$mRNA}{'UTR_3'}{$b} <=> $mRNAs{$mRNA}{'UTR_3'}{$a} } keys %{$mRNAs{$mRNA}{'UTR_3'}} ) {

			                        my $start = $utr;

			                        my $end = $mRNAs{$mRNA}{'UTR_3'}{$utr};

			                        my $utrseq = reverse(lc($seq->subseq($start, $end)));

			                        $utrseq =~ tr/ACGT/TGCA/;

			                        push(@genomic, $utrseq);
			                        push(@transcript, $utrseq);

	#	                                print GENOMIC2 $utrseq;
	#	                                print TRANSCRIPT2 $utrseq;

			                } #end utr foreach

			        } #end printing 3' UTR

		        #sequence is on minus
			} else {

			        my $firstExon = "true";

			        my $prevStart = 0;

			        my $codon = "";

			        # print 5' UTR

			       if(defined($mRNAs{$mRNA}{'UTR_5'}) ) {

			                foreach my $utr (sort { $mRNAs{$mRNA}{'UTR_5'}{$b} <=> $mRNAs{$mRNA}{'UTR_5'}{$a} } keys %{$mRNAs{$mRNA}{'UTR_5'}} ) {

			                        my $start = $utr;

			                        my $end = $mRNAs{$mRNA}{'UTR_5'}{$utr};

			                        my $utrSeq = reverse(lc($seq->subseq($start, $end)));

			                        $utrSeq =~ tr/ACGT/TGCA/;

			                        push(@genomic, $utrSeq);
			                        push(@transcript, $utrSeq);

	#	                                print GENOMIC2 $utrSeq;
	#	                                print TRANSCRIPT2 $utrSeq;

			                } #end utr foreach

			        } #end printing 5' UTR

				my $codingSequence = "";

			        foreach my $cds (sort { $mRNAs{$mRNA}{'CDS'}{$b} <=> $mRNAs{$mRNA}{'CDS'}{$a} } keys %{$mRNAs{$mRNA}{'CDS'}} ) {

			                if ( $firstExon eq "false" ) {

			                        my $end = $prevStart - 1;

			                        my $start = $mRNAs{$mRNA}{'CDS'}{$cds} + 1;

			                        my $intron = reverse(lc($seq->subseq($start, $end)));

			                        $intron =~ tr/acgt/tgca/;

			                        push(@genomic, $intron);

	#	                                print GENOMIC2 $intron;

			                } #end print flank and intron

			                my $start = $cds;

			                my $end = $mRNAs{$mRNA}{'CDS'}{$cds};

			                my $exon = reverse(uc($seq->subseq($start, $end)));

			                $exon =~ tr/ACGT/TGCA/;

			                #Print exon

			                push(@genomic, $exon);
			                push(@transcript, $exon);
			                push(@CDS, $exon);

	#	                        print GENOMIC2 $exon;
	#	                        print TRANSCRIPT2 $exon;
	#	                        print CDS2 $exon;

			                $codingSequence = $codingSequence.$exon;

			                $firstExon = "false";

			                $prevStart = $cds;

			        } #end print CDS

		                #print protein

		                my $codonStart = 0;

		                $codon = substr($codingSequence, $codonStart, 3);
		                my $aa = $CodonTable->translate_strict($codon);

		                while ($aa ne "") {

		                        push(@protein, $aa);

	#                                print PROTEIN $aa;

		                        ++$remaining_protein;
		                        ++$size_protein;

		                        $codonStart += 3;

		                        $codon = substr($codingSequence, $codonStart, 3);

		                        $aa = $CodonTable->translate_strict($codon);

		                } #end while

			        # print 3' UTR genomic

			        if(defined($mRNAs{$mRNA}{'UTR_3'}) ) {

			                foreach my $utr (sort { $mRNAs{$mRNA}{'UTR_3'}{$b} <=> $mRNAs{$mRNA}{'UTR_3'}{$a} } keys %{$mRNAs{$mRNA}{'UTR_3'}} ) {

			                        my $start = $utr;

			                        my $end = $mRNAs{$mRNA}{'UTR_3'}{$utr};

			                        my $utrSeq = reverse(lc($seq->subseq($start, $end)));

			                        $utrSeq =~ tr/ACGT/TGCA/;

			                        #Print utr

			                        push(@genomic, $utrSeq);
			                        push(@transcript, $utrSeq);

	#	                                print GENOMIC2 $utrSeq;
	#	                                print TRANSCRIPT2 $utrSeq;

			                } #end utr foreach

			        } #end printing 3' UTR genomic

			} #end print plus and minus strands

	#	        print GENOMIC2 "\n";
	#	        print TRANSCRIPT2 "\n";
	#	        print CDS2 "\n";
	#	        print PROTEIN2 "\n";

			print GENOMIC ">$mRNA genomic ".length(join("",@genomic))."bp reference_sequence $mRNA_info[0] ($mRNA_info[2] strand)\n".join("",@genomic)."\n";
			print TRANSCRIPT ">$mRNA transcript ".length(join("",@transcript))."bp reference_sequence $mRNA_info[0] ($mRNA_info[2] strand)\n".join("",@transcript)."\n";
			print CDS ">$mRNA CDS ".length(join("",@CDS))."bp reference_sequence $mRNA_info[0] ($mRNA_info[2] strand)\n".join("",@CDS)."\n";
			print PROTEIN ">$mRNA peptide ".length(join("",@protein))."aa reference_sequence $mRNA_info[0] ($mRNA_info[2] strand)\n".join("",@protein)."\n";

		#sequence doesn't found on index
		} else {

			print "Sequence for $mRNA_info[0] not found in fasta file for mRNA $mRNA\n";
			print ERRORS "Sequence for $mRNA_info[0] not found in fasta file for mRNA $mRNA\n";

		} #end else

	} else {

		print ERRORS "Undefined parent for mRNA $mRNA\n";

	} #end else undefined data

} #end printing mRNAs

exit;

#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

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

my $idx = Bio::Index::Fasta->new(
                                 '-filename' => "$ARGV[0].idx",
                                 '-write_flag' => 1
                                );
$idx->make_index("$ARGV[0]");

my %features;

my $count = 0;

open(GFF, $ARGV[1]);

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	$GFF_line[8] =~ s/;//;
	
	$features{$GFF_line[8]} = [$GFF_line[0],$GFF_line[1],$GFF_line[3],$GFF_line[4],$GFF_line[5],$GFF_line[6]];

	++$count;

} #end while

print "$count features's were loaded and will be use for sequence extraction\n\n";

close(GFF);


open(FASTA, ">$ARGV[1].fasta");
open(ERRORS, ">$ARGV[1].errors.txt");

my $analyzedfeature = 0;

foreach my $feature (keys %features) {

	if( ($analyzedfeature/1000) =~ m/^\d+$/) {

		print "Done with $analyzedfeature repeats's of $count\n";

	} #end if integer

	++$analyzedfeature;

	my @feature_info = @{$features{$feature}};

	delete $features{$feature};

	my $seq = $idx->fetch($feature_info[0]);

        if(defined($seq)) {

		my $start = $feature_info[2];
		my $end = $feature_info[3];

		print FASTA ">", $feature_info[0], ":", $start, "..", $end, " ", $feature, " source=", $feature_info[1], " score=", $feature_info[4]," (strand ", $feature_info[5]. ")\n";

                #sequence is on plus
	        if ($feature_info[5] eq "+") {

			my $featureSeq = lc($seq->subseq($start, $end));
			print FASTA $featureSeq, "\n";

                #sequence is on minus
	        } else {
			
			my $featureSeq = reverse(lc($seq->subseq($start, $end)));
			$featureSeq =~ tr/ACGTacgt/TGCAtgca/;

			print FASTA $featureSeq, "\n";

	        } #end print plus and minus strands

        #sequence doesn't found on index
        } else {

        	print "Sequence for $feature_info[0] not found in fasta file for repeat $feature\n";
        	print ERRORS "Sequence for $feature_info[0] not found in fasta file for repeat $feature\n";

        } #end else

} #end printing mRNAs

exit;

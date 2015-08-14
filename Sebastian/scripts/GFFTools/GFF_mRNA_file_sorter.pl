#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;

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
	print "perl GFF_mRNA_GenomicSequence_Extractor.pl <InputGFF> <OutputGFF>\n";
        exit 0;

} #end if

my %genes;

open(GFF, $ARGV[0]);

while (my $line = <GFF>) {

	if($line =~ /^#/) { next }

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	my $name;

	foreach my $feature (@features) {
		if ($feature =~ m/ID=/) {
			my @value = split("=", $feature);
			$name = $value[1];
			last;
	        } elsif ( $feature =~ m/Parent/) {
			my @value = split("=", $feature);
			$name = $value[1];
	        } #end else
	} #end foreach

	my @parent = split("\\.", $name);

	$genes{$parent[0]}{$name}{$GFF_line[2]}{$GFF_line[3]} = $line;

} #end while

print STDERR scalar(keys %genes)," genes's were loaded and will be sorted\n\n";

close(GFF);

my @typeList = ("exon","CDS");

open(RESORTED, ">$ARGV[1].gff3") or die "can't open file $ARGV[1].gff3";
print RESORTED "##gff-version 3\n";

foreach my $gene (sort keys %genes) {

	if (defined($genes{$gene}{$gene}{"gene"})) {
        	foreach my $key (keys %{$genes{$gene}{$gene}{"gene"}} ) { print RESORTED $genes{$gene}{$gene}{"gene"}{$key}, "\n" }
	} else { print "Missing gene tag for $gene\n"}

	delete $genes{$gene}{$gene};

	foreach my $mRNA (sort keys %{$genes{$gene}}) {

        	my @mRNAkey = keys %{$genes{$gene}{$mRNA}{"mRNA"}};

		print RESORTED $genes{$gene}{$mRNA}{"mRNA"}{$mRNAkey[0]}, "\n";

		my @mRNAInfo = split("\t",$genes{$gene}{$mRNA}{"mRNA"}{$mRNAkey[0]});

		foreach my $type (@typeList) {

                	if($mRNAInfo[6] eq "+") {
	                        foreach my $line (sort {$a <=> $b} keys %{$genes{$gene}{$mRNA}{$type}}) {
	                                print RESORTED $genes{$gene}{$mRNA}{$type}{$line}, "\n";
	                        } #end foreach line
                        } else {
 	                        foreach my $line (sort {$b <=> $a} keys %{$genes{$gene}{$mRNA}{$type}}) {
	                                print RESORTED $genes{$gene}{$mRNA}{$type}{$line}, "\n";
	                        } #end foreach line
			} #end else

		} #end foreach type
	} #end foreach mRNA
} #end printing genes

exit;

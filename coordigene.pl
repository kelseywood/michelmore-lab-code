#!/usr/bin/perl
# coordigene.pl by kelsey
use strict; use warnings;
use FAlite;

#gets genomic sequence given genome (FASTA format), FASTA identifier of genes and their start/stop coordinates

die "usage: coordigene.pl <genome> <tab delimited file w contig ID and coords>" unless @ARGV == 2;

#file format: contig start stop

my %seqs;
open(my $genome, "<$ARGV[0]") or die "error reading $ARGV[0] for reading";
        my $fasta = new FAlite($genome);
        while (my $entry = $fasta->nextEntry) {
		    my ($id) = $entry->def =~ /^>(\S+)/;
		    $seqs{$id} = $entry->seq;
	    }
close($genome);
	    
open(my $coord, "<$ARGV[1]") or die "error reading $ARGV[1] for reading";
while (<$coord>) {
    my ($contig, $start, $stop) = split;
    my @ids = keys(%seqs);
    if ($contig ~~ @ids){
        my $seq = $seqs{$contig};
        if ($start < $stop){
            my $orf = substr($seq, ($start-1), ($stop-$start+1)); #substr(string, start, length)
            print ">$contig\_$start-$stop\_fwd\n";
            print "$orf\n";
        }
        else {
            my $orf = substr($seq, ($stop-1), ($start-$stop+1));
            my $rev = reverse($orf);    #reverse sequence
            $rev =~ tr/a-z/A-Z/;        #all uppercase
            $rev =~ tr/ATCGYRN/TAGCRYN/;  #complement sequence
            print ">$contig\_$stop-$start\_rev\n";
            print "$rev\n";
        }
    }
}
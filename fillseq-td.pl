#!/usr/bin/perl
# fillseq.pl by kelsey
use strict; use warnings;
use FAlite;

#this script will take a stockholm MSA from HMMer and fill out the rest of the sequence to be a given length
#requires a FASTA of the full sequences

die "usage: fillseq.pl <multiple seq alignment> <FASTA file with full seqs> <seq length>" unless @ARGV == 3;
my $length = $ARGV[2];

#get gene names and position of motif from msa file
open(my $msa, "<$ARGV[0]") or die "error reading $ARGV[0] for reading";
while (<$msa>){
    my $line = $_;
    next if $line =~ /^\s*#/;  #ignore lines w hashtags from stockholm
    next unless /\S/;      #ignore lines that are blank spaces
    if (/^\s*\/\//) {}     #end of file (ends with // in Stockholm)
    else{
        my ($des, $motif) = split;  #splits description and sequence by tab
        fill($des, $motif);         #send to subroutine fill
    }
}

#subroutine to fill out sequence for each motif
sub fill{
    my ($des, $motif) = @_;
    my ($name, $pos) = $des =~ /(\w+)\/(\d+)/; #extract the gene name & pos from description -(\w+ = one or more letters, numbers or _; \/ = forward slash; \d+ digits)
    
    open(my $genome, "<$ARGV[1]") or die "error reading $ARGV[1] for reading";
    my $fasta = new FAlite($genome);
    while (my $entry = $fasta->nextEntry) {
        my $def = $entry->def;            # get definition line from object
	    my $seq = $entry->seq;            # get sequence from object
	    if ($def =~ /$name/){
            #print "$name\t$pos\n";
            my $dashes = count_dashes($motif);
            my $start = $pos-$dashes-1;
            if ($start < 0){$start = 0} #if there are more dashes than amino acids in the protein (ie the motif is at the n-term)
            my $filledseq = substr($seq, $start, $length);
            print "$des\t";
            print "$filledseq\n";
        }
    }
}

#subroutine to count dashes (actually only leading dashes)
sub count_dashes{
    my ($motif) = @_;
    my $count;
    if ($motif =~ /^\w/){$count = 0} #if motif starts with a character instead of a dash then dash count = 0
    else{
        my ($dashes) = $motif =~ /(^[\-\.]+)\w+/; #number of leading dashes OR dots
        $count = ($dashes =~ tr/[\-\.]/\-/);    #use transliterate to count dashes/dots (replace with -)
    }
    return $count;
}

   
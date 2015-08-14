#!/usr/bin/perl
# randoseqset_krb.pl by Keith Bradnam
# gets random sequences from one FASTA file, using target lengths from another FASTA file
use strict; 
use warnings;
use FAlite;


die "usage: $0 <FASTA file to imitate> <FASTA of sequences>" unless @ARGV == 2;
my ($seqset, $genome) = @ARGV;


############################################################################################################
# Read genome FASTA file and save data to a hash of arrays
# hash key is sequence length, value is array of sequences (and headers) with that length
############################################################################################################

my %seqs;
open(my $in, "<", $genome) or die "error reading $genome\n";
my $fasta = new FAlite($in);
while (my $entry = $fasta->nextEntry) {
    my $def = $entry->def;
    my $seq = $entry->seq;
    my $len = length($seq);
    push(@{$seqs{$len}}, "$def\n$seq");
}
close($in);



############################################################################################################
# Generate list of lengths from target FASTA file
############################################################################################################

my @target_seq_lengths = extract_target_seq_lengths($seqset);
my $n = @target_seq_lengths;

for (my $i = 1; $i <= $n; $i++){
    #warn "Choosing random sequence from genome $i of $n\n";

    # select a random length from all target lengths
    my $target_length = $target_seq_lengths[rand($n)];
    
    # how many sequences in the genome are there at this length?
    my $possible_matches = @{$seqs{$target_length}};
    #print "$possible_matches\n";
    
    # choose one sequence of that length randomly
    my $random_seq = ${$seqs{$target_length}}[rand($possible_matches)];
    print "$random_seq\n";
    my $len = length($random_seq);
    #print "$len\n";
}



exit;

####################################
#
#     S U B R O U T I N E S
#
####################################


sub extract_target_seq_lengths{
    my ($file) = @_;
    my @lengths;

    open(my $in, "<", $file) or die "error reading $file\n";
    my $fasta = new FAlite($in);
    while (my $entry = $fasta->nextEntry) {
        my $def = $entry->def;
        my $seq = $entry->seq;
        my $len = length($seq);
        push(@lengths, $len);
    }
    close($in);
    return(@lengths);
}

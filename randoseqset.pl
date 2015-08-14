#!/usr/bin/perl
# randoseqset.pl by kelsey
use strict; use warnings;
use FAlite;

#gets random sequences from a FASTA file
#tries to match sequence length distribution as original file (not perfect)
#currently using mean as lower limit and maximum value as upper limit
#input = fasta and number of sequences you want to pull out
#output = random sequences, same number as in the file you want to mimic, within your length dist.

die "usage: $0 <FASTA file to imitate> <FASTA of sequences> <outfile name>" unless @ARGV == 3;
my $seqset = $ARGV[0];
my $genome = $ARGV[1];
my $outname = $ARGV[2];

#get sequence range
my ($min, $max) = seq_range($seqset);

#read genome and save to hash
my %seqs;
open(my $gen, "<$genome") or die "error reading $ARGV[0] for reading";
        my $fasta = new FAlite($gen);
        while (my $entry = $fasta->nextEntry) {
		    my ($id) = $entry->def =~ /^>(\S+)/;
		    $seqs{$id} = $entry->seq;
	    }
close($genome);

#only pick sequences within the same size range
#make new hash for sequences within the target range
my %targetseqs;
foreach my $key (keys %seqs){
    my $seq = $seqs{$key};
    my $len = length($seq);
    if ($min < $len && $len < $max){
       $targetseqs{$key} = $seq;
       #print "$len\n";
       #print "\>$key\n";
       #print "$seq\n";
    }
}

print "Target sequence range is between the mean to the maximum sequence size\n";
print "Sequence length range target: $min - $max\n";   #print this to stdout

#count number of sequences in FASTAs
my $settotal = seq_counter($seqset);
my @keys = keys(%targetseqs);
my $hashtotal = @keys;

print "Sequences in genome within this size range = $hashtotal\n";
print "Total number of sequences to generate: $settotal\n";

#generate array of random numbers
my @rando;
for (my $i = 0; $i < $settotal; $i++){
    push @rando, int(rand($hashtotal));
}

#pick random sequences from the array of sequences of the target size

open(my $out, ">$outname") or die "error creating random fasta";

my $count; #for sequence number
foreach my $target (keys %targetseqs){
    $count++;
    if ($count ~~ @rando){  #if count matches ANY of the numbers in random number list
        my $targetseq = $seqs{$target};
        #print $out ">$target\n";
        #print $out "$targetseqs{$target}\n";
        #my $len = length($targetseq);
        #print $out "$len\n";
        #print "\>$key\n";
        #print "$seq\n";
    }
}
close($out);

print "Random sequences saved to $outname\n";


####(subroutines)####

sub seq_counter{
	my ($file) = @_;
	my $total = 0;
	open(IN, "$file") or die "error reading $file for reading";
    while (<IN>){
        chomp;
        if (m/>/){
            $total++;
        }
    }
	return $total;	
}

sub seq_range{
    my ($file) = @_;
    my (@lengths) = seq_length($file);
    my $count = @lengths;
    my $sum = 0;
    foreach my $length (@lengths){
        print "$length\n";
        $sum += $length;
        #print "$sum\n";
    }
    my $median = seq_median(@lengths);
    #print "Median = $median\n";
    my $mean = $sum/$count;
    #actual longest sequence
    my @sorted = sort {$a <=> $b} @lengths;
    my $maxseq = $sorted[-1]; #last item of sorted lengths
    my $var = 0;
    foreach my $length (@lengths){
        my $dif = ($mean - $length)**2;
        $var += $dif;
    }
    #print "Mean = $mean\n";
    my $samplevar = $var/$count;
    my $stdev = sqrt($samplevar);
    #print "Standard Dev = $stdev\n";
    my $min = $mean - $stdev;
    my $max = $maxseq;
    return ($min, $maxseq);
}

sub seq_length{
	my ($file) = @_;
	my @lengths;
	open(IN, "$file") or die "error reading $file for reading";
    while (<IN>){
        chomp;
        my $line = $_;
        next if $line =~ m/^>/;
        if ($line =~ m/\S*/){
            my $length = length($line);
            push (@lengths,$length);
        }
    }
	return @lengths;	
}

sub seq_median{
    my (@lengths) = @_;
    my $count = @lengths;
    my @sorted = sort(@lengths);
    my $median = 0;
    if (@sorted % 2 == 1){ #to check if odd or even (if odd -> remainder = 1)
        $median=splice(@sorted,($count/2),1);
    }
    else {
        $median= ((splice(@sorted,($count/2),1)) + splice(@sorted,(($count/2)-1),1))/2;
    } 
    return $median;
}
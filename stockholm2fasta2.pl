#!/usr/bin/perl -w
#modified by kelsey to not print sequence by number of columns

my $gapped = 0;
my $progname = $0;
$progname =~ s/^.*?([^\/]+)$/$1/;
 
my $usage = "Usage: $progname [<Stockholm file(s)>]\n";
$usage .=   "             [-h] print this help message\n";
$usage .=   "             [-g] write gapped FASTA output\n";
$usage .=   "             [-s] sort sequences by name\n";

# parse cmd-line opts
my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
        die $usage;
    } elsif ($arg eq "-g") {
        $gapped = 1;
    } elsif ($arg eq "-s"){
        $sorted = 1;
    } else {
        push @argv, $arg;
    }
}
@ARGV = @argv;
 
my @seqorder = ();
 
my %seq;
while (<>) {
    next unless /\S/;
        next if /^\s*\#/;
    if (/^\s*\/\//) { printseq() }
    else {
        chomp;
        my ($name, $seq) = split;
        $seq =~ s/[\.\-]//g unless $gapped;
        push @seqorder, $name unless exists $seq{$name};
        $seq{$name} .= $seq;
    }
}
printseq();
 
sub printseq {
    
    @seqorder = sort @seqorder if $sorted;
 
    foreach $key (@seqorder) {
            print ">$key\n";
            print "$seq{$key}\n";
    }
 
    %seq = ();
    @seqorder = ();
}
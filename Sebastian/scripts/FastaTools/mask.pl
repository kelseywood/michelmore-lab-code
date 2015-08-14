#!/usr/bin/perl -w
use strict;
use List::Util qw(max);
use Getopt::Long;

#Written by CB

my $fa;
my $gff;
my $mask="N";

# Will mask every contig identified in column 1 at the ranges specified in columns 4 and 5. If two overlapping ranges are specified, the union of the ranges will be masked.
# Gff file is 1 based, i.e. first position in a contig is 1, not 0.
# USAGE: mask.pl --fa=[lidas.fa] --gff=[lidas.gff] (optionally --mask=[mask char, e.g. "N" or "n"])

GetOptions("fa=s"=>\$fa, "gff=s"=>\$gff, "mask=s"=>\$mask);
if(!defined($fa)||!defined($gff)){die "Damn it Lida!! You have to set --fa and --gff."}

# Build a hash of positions to be masked
my %reph;
open(IN, "<$gff");
while(<IN>){
chomp;
my @arr=split();
my $span = $arr[4]-$arr[3]+1;
if($span<0){die "Damn it Lida!! Column 5 must always be greater than column 4 in your GFF."}
my $currspan=0;
if(defined($reph{$arr[0]}{$arr[3]})){$currspan=$reph{$arr[0]}{$arr[3]};}
my $mn = max($span, $currspan);

# Debug:
#print "$arr[0]\t$arr[3]\t$mn\n";

$reph{$arr[0]}{$arr[3]}=$mn;

}
close(IN);

# Read the input, mask, and print to output
open(IN, "<$fa");
my $id=0;
while(<IN>){
chomp;
if(substr($_,0,1)=~m/\>/){
my @arr=split(/ /, $_); 
$id=$arr[0]; $id =~ s/\>//;
print "$_\n";

}else{
my $line = $_;
my $len=length($line);
if(defined($reph{$id})){foreach my $key (keys %{$reph{$id}}){
# Replace that portion of the string with N's
my $start = $key - 1;
my $span = $reph{$id}{$key};
for(my $i=0; $i<$span; $i++){
substr($line, ($start+$i), 1) = $mask;
}
}}

# Debug:
#my $toprint = substr($line, 0, 100);
#print "$toprint\n";

print "$line\n";

}
}
close(IN);




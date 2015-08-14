#!/usr/bin/perl
# 3of5parser.pl by kelsey
use strict; use warnings;

die "usage: 3of5parser.pl <text file with 3 of 5 output> <output name>" unless @ARGV == 2;

#outline:
#get ORF names - always the line after a number followed by a period
#get sequences - starts with 001
#get position of motif - ends with :
#get sequence of motif - below position info
#get name of motif - three lines below position info


open(my $input, "<$ARGV[0]") or die "error reading $ARGV[0] for reading";
open(my $out, ">$ARGV[1]") or die "error creating $ARGV[1]";

my $motifcount = 0; #if there are multiple motifs per protein, then we want to add a few tabs to align with the protein name

while (my $line = <$input> ) {
    if ($line =~ m/^\d+\./){        #digits and then a period
        $motifcount = 0;        #reset counter
        my $nextline = <$input>; #get next line (the name)
        chomp $nextline;         #take off the newline
        print $out "$nextline\t";     #print the name
     }
     elsif($line =~ m/^001/){    #get first row of sequence (starts with 001)
        $line =~ s/[\d+\s]//g;   #get rid of numbers and whitespace
        my $nextline = <$input>; #get second row of sequence
        $nextline =~ s/[\d+\s]//g;
        print $out "$line$nextline\t";
     }
     elsif($line =~ m/^\d+.*\:$/){ #find position info, digits and ends with colon
        $motifcount++; #add one to motif counter
        ($line) = $line =~ m/^(\d+).*\:$/; #take out just the start position
        my $nextline = <$input>;
        chomp $nextline;
        <$input>; #next line
        <$input>; #next line
        my $pattern = <$input>;
        chomp $pattern;
        if($motifcount == 1){
            print $out "$line\t$nextline\t$pattern\n";
        }
        else{
            print $out "\t\t$line\t$nextline\t$pattern\n"; #print tabs if this is not the first motif
        }
     }
}

close($input);
close($out);

    #print if /\d+\./ .. /^\s*$/; #where $/ is the newline character

        #foreach (1..2) { }
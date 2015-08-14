#!/usr/bin/perl

use strict;
use warnings;

if(@ARGV < 1) {

	print "Program usage: at least three arguments required, please input them\n";
        print "First argument is the output filename.\n";
        print "Second argument is the file with the reference information.\n";
        print "Later arguments, count data to merge on the base file (format 'GeneName\tCount').\n";
        exit 0;
}

my $numFiles = scalar @ARGV;

my %posData = ();

open (OUTPUT, ">$ARGV[0]") or die "Can't create $ARGV[0]\n";

print OUTPUT "\n";

open (REFERENCE, $ARGV[1]) or die "Can't open $ARGV[1]\n";

while ( my $referenceLine = <REFERENCE>) { print OUTPUT $referenceLine }

close(OUTPUT);

for(my $file = 2; $file < $numFiles; ++$file) {

	open(INPUT, $ARGV[$file]) or die "Can't open $ARGV[$file]\n";

        my %countData;

        while (my $line = <INPUT>) {

        	chomp($line);

                my @countInfo = split("\t", $line);

                $countData{$countInfo[0]} = $countInfo[1];

        } #end while

       	open(BASEFILE, $ARGV[0]) or die "Can't open file $ARGV[0]\n";

        open(RESULTS, ">TemporaryCounts.txt");

        my $firstLine = <BASEFILE>;

        chomp($firstLine);

        print RESULTS $firstLine, "\t", $ARGV[$file],"\n";

       while (my $baseLine = <BASEFILE>) {

	        chomp($baseLine);

	        my @baseData = split("\t", $baseLine);

                if(defined($countData{$baseData[0]})) {

	                print RESULTS $baseLine, "\t", $countData{$baseData[0]}, "\n";

                } else {

                        print RESULTS $baseLine, "\t0\n";

                } #end else

        } #end while

        close(RESULTS);
        close(BASEFILE);

        if(system("mv TemporaryCounts.txt $ARGV[0]") != 0) {
		print "error executing system command for rename output file file ", $ARGV[$file], "\n";
	}  #end if

} #end for

exit;

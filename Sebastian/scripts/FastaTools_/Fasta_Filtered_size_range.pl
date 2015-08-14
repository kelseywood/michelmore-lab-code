#!/usr/bin/perl


if(@ARGV < 3) {

	print "Not enough arguments where inputed\n";
	print "Please input the needed arguments\n";
	print "Usage; Fasta_Filtered_size_range.pl <FASTA> Min_size Max_size\n";
	exit;

}#end if

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_$ARGV[1]-$ARGV[2].txt");

while($line = <FASTAINPUT>) {

	chomp($line);

        $sequenceName = substr($line,1,200);

        $sequence = <FASTAINPUT>;

        chomp($sequence);

        if( (length($sequence) >= $ARGV[1]) && (length($sequence) <= $ARGV[2]) ) {

	        print OUTPUT ">", $sequenceName, "\n";
	        print OUTPUT $sequence, "\n";

	} #end if

} #end while

close(FASTAINPUT);
close(OUTPUT);

exit;


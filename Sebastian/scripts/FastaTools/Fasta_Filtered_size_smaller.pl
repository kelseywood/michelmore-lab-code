#!/usr/bin/perl

if (!(defined($ARGV[1])) ) {

	print "Please input Fasta file where to filter and size threeshold\n";
	die "Missing arguments\n";

} #end if missing arguments

open(FASTAINPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0]_smallerthan_$ARGV[1].fasta");

my $header;
my $sequence = "";

while(my $line = <FASTAINPUT>) {

	chomp($line);
	$line =~ s/\015//g;

	if ($line =~ /^>/) {

		if( (length($sequence) < $ARGV[1]) && defined($header)) {

			print OUTPUT $header, "\n";
                	print OUTPUT $sequence, "\n";

                } #end else

	        $header = $line;
		$sequence = "";

        } else {

		chomp($line);
		$sequence .= $line;

        } #end else

} #end while

if( (length($sequence) < $ARGV[1]) && defined($header)) {

	print OUTPUT $header, "\n";
        print OUTPUT $sequence, "\n";

} #end else

close(FASTAINPUT);
close(OUTPUT);

exit;


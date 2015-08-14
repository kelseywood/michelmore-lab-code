#!/usr/bin/perl

open(INPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0].stats.txt");

my $prevScaffold = "";
my %data = ();

$line = <INPUT>;
chomp($line);
@info = split("\t", $line);
push(@{$data{$info[1]}}, $info[2]);

$prevScaffold = $info[0];

while($line = <INPUT>) {

	chomp($line);

        @info = split("\t", $line);

        if ($prevScaffold eq $info[0]) {

        	push(@{$data{$info[1]}}, $info[2]);
                $prevScaffold = $info[0];

        } else {

        	print OUTPUT $prevScaffold;

                foreach $LG (keys %data) {

	                my $total = 0;
	                my $numPos = 0;

                        print OUTPUT "\t";
                	print OUTPUT $LG;

                        foreach $pos (@{$data{$LG}}) {
                        	$total += $pos;
                                ++$numPos;
                        } #end foreach

                        my $average = $total/$numPos;

                        my $sqtotal = 0;

                        foreach $pos (@{$data{$LG}}) {
                        	$sqtotal += ($average-$pos)**2;
                        } #end foreach

                        $std = ($sqtotal/$numPos)**0.5;

			print OUTPUT "\t", $numPos;
                        print OUTPUT "\t";
                        printf OUTPUT ("%.2f", $average);
                        print OUTPUT "\t";
                        printf OUTPUT ("%.2f", $std);

                } #end foreach

                print OUTPUT "\n";

	        %data = ();

	        push(@{$data{$info[1]}}, $info[2]);
	        $prevScaffold = $info[0];

          } #end else

} #end while

close(INPUT);
close(OUTPUT);
exit;
#!/usr/bin/perl

use strict; use warnings;

if( !defined($ARGV[0]) ) { print "Please provide a fastq file\n"; die "missing argumetns\n"}

open(FASTQ, $ARGV[0]);

my $file = $ARGV[0];

$file =~ s{\.[^.]+$}{};

open(COUNTS, ">$file.Ncounts.txt");

my $stat = "";

my @Ncounts;

my %distribution;

while(my $header = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	chomp($sequence);

	my $length = length($sequence);

	my $numN = $sequence =~ tr/(N|n)//;

	chomp($header);

	print COUNTS $header, "\t", $numN, "\t", sprintf( "%.2f", ($numN/$length)*100), "\n";

	push(@Ncounts, $numN);

	++$distribution{$numN};

} #end while

open(STATS, ">$file.Nstats.txt");

my $total;

my $max = 0;

foreach my $readNs (@Ncounts) {

	$total += $readNs;

	if($max < $readNs) { $max = $readNs}

} #end if

my $aver = $total/scalar(@Ncounts);

print STATS "Maximum number of N's in a reads: ", $max, "\n"; 
print STATS "Average number of N's per read: ", sprintf( "%.2f", $aver ), "\n\n";
print STATS "Distribution of N's per read\n";
print STATS "Number of N's\tRead count\n";

foreach my $N (sort {$b<=>$a} keys %distribution ) {

	print STATS $N, "\t", $distribution{$N}, "\n";

} #end foreach printing distribution 



close(FASTQ);
close(COUNTS);
exit;
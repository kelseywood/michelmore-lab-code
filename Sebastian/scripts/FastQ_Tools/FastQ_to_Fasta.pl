#!/usr/bin/perl

use strict; use warnings;

if( !defined($ARGV[0])) { die "Missing arguments, please provide a FastQ file to convert/n"}

open(FASTQ, $ARGV[0]);

my $fasta = $ARGV[0];

$fasta =~ s{\.[^.]+$}{};

open(FASTA, ">$fasta.FASTA");

my $stat = "";

while(my $line = <FASTQ>) {

	chomp($line);

	if($line =~ /^@*.\//) {

		$line =~ s/@/>/;

		print FASTA $line, "\n";

		$stat = "sequence";

	} elsif($line =~ /^\+/) {

		$stat = "quality";

	} else {

		if ($stat eq "sequence") {

			print FASTA $line, "\n";

		} #end if its sequence

	} #end if its else

#        $firstElement = substr($line, 0, 1);

#        if ($firstElement eq ">") {

#        	print FASTA $line;
#                print FASTA "\n";

#                $stat = "sequence";

#        } elsif ($firstElement eq "@") {

#        	$stat = "quality";

#        } else {

#        	if ($stat eq "sequence") {

#                	print FASTA $line;
#                        print FASTA "\n";

#                } #end if

#        } #end else

} #end while







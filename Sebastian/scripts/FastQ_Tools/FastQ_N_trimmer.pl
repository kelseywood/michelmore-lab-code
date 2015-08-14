#!/usr/bin/perl

use strict; use warnings;

if( !defined($ARGV[0]) ) {
	print "Please provide a fastq file and a maximum number of N's allowed per read\n";
	print "FastQ_N_filterer <Fastq file>\n";
	die "missing argumetns\n"
}

open(FASTQ, $ARGV[0]);

my $file = $ARGV[0];

$file =~ s{\.[^.]+$}{};

open(CLEAN, ">$file.cleanReads.fastq");
open(BAD, ">$file.badReads.fastq");
open(TRIM, ">$file.trimmedReads.fastq");

while(my $header = <FASTQ>) {

	my $sequence = <FASTQ>;
	my $header2 = <FASTQ>;
	my $qual = <FASTQ>;

	chomp($sequence);

	my $length = length($sequence);

	my $numN = $sequence =~ tr/(N|n)//;

	if ($numN == 0 ) {
		print CLEAN $header;
		print CLEAN $sequence, "\n";
		print CLEAN $header2;
		print CLEAN $qual;
	} elsif ( ($numN/$length) > 0.4 ) {
		print BAD $header;
		print BAD $sequence, "\n";
		print BAD $header2;
		print BAD $qual;
	} else {

		my @seq = split("", $sequence);

#		my @blocks;

		my $lastGoodBase = 1;

#		my $startBlock = 1;

		my $Ncount = 0;

		my $totalN = 0;

		for (my $i = 1; $i < $length; ++$i) {

			if ($seq[$i] =~ /(N|n)/) {

				++$Ncount;
				++$totalN;

			} else {

				if ($Ncount > 3) {

					if ($totalN-$Ncount > 4) {

						print BAD $header;
						print BAD $sequence, "\n";
						print BAD $header2;
						print BAD $qual;

					} else {

						my $trimSeq = substr($sequence, 1, $lastGoodBase);
						my $trimQual = substr($qual, 1, $lastGoodBase);

#						push(@blocks, ($startBlock, $lastGoodBase));
#						$startBlock = $i;

						print TRIM $header;
						print TRIM $trimSeq, "\n";
						print TRIM $header2;
						print TRIM $trimQual, "\n";

					} #end else

					last;

				} #end if big gap

				$lastGoodBase = $i;
				$Ncount = 0;

			} #ene else

		} #end for seq				

	} #end else	

} #end while


close(FASTQ);
close(CLEAN);
close(BAD);
close(TRIM);
exit;

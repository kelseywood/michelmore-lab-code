#!/usr/bin/perl
use strict;
use warnings;

#########################
#
# Diff from v0
#
# Filtering base in bitwise flag, if BWA mem use, flag -M required.
#
#########################
my $sam = $ARGV[0];

if (!defined($sam)) {

	print "Error, not sam input was entered. Please enter at least one SAM file for analysis\n";
	print "If bwa mem used, flag -M required\n";
	print "If other alignerd used, flagging for secundary aligments it's required\n\n";
	print "Usage:\n";
	print "SAM_filterer_v2.pl <SAM FILE> <SAM FILE2> <SAM FILE3> ... <SAM FILEN>\n";
	die "Missing arguments\n";

} #end if

foreach my $sam (@ARGV) {

	open(SAM, $sam) or die "Can't open $sam\n";

	open(SINGLE, ">$sam.filterMappedSingle.sam");
	open(PAIRED, ">$sam.filterMappedPaired.sam");
	open(STATS, ">$sam.stats.txt");

	print "Analyzing $sam\n";

	my $numMultipleMappings = 0;

	my $numMapped = 0;
	my $numUnMapped = 0;

	my $numUnPair = 0;
	my $numUnSing = 0;
	my $numMaPair = 0;
	my $numMaSing = 0;

	while (my $samLine = <SAM>) {

		chomp($samLine);

		if( !($samLine =~ /^@/) ) {

			my @alignment = split("\t", $samLine);

			my $flags = join("|", getBitFlag($alignment[1]));

	#		print $flags, "\n";

			if ($flags =~ m/read_unmapped/) {

				if($flags =~ m/mate_unmapped/) { ++$numUnPair }
				else { ++$numUnSing } #end if mate unmapped
				++$numUnMapped;

			} elsif ($flags =~ m/not_primary_alignment/) {

				++$numMultipleMappings;

			} elsif ($flags =~ m/first_in_pair/) {

				$alignment[0] .= "/1";

				if($flags =~ m/mate_unmapped/) {
					++$numMaSing;
					print SINGLE join("\t", @alignment), "\n";
				} else {
					++$numMaPair;
					print PAIRED join("\t", @alignment), "\n";
				} #end if mate unmapped

				++$numMapped;

			} elsif ($flags =~ m/second_in_pair/) {

				$alignment[0] .= "/2";

				if($flags =~ m/mate_unmapped/) {
					++$numMaSing;
					print SINGLE join("\t", @alignment), "\n";
				} else {
					++$numMaPair;
					print PAIRED join("\t", @alignment), "\n";
				} #end if mate unmapped

				++$numMapped;

			} #end filtering conditions

		


		} #end if its an alignment line

	} #end while
	close(SAM);

	print "Done filtering reads for $sam, will print stats\n";

	my $numTotalReads = $numMapped + $numUnMapped;

	print STATS "\nAlignment Stats for $sam\n\n";

	print STATS "Totals\n";
	print STATS "Total number of reads		  ", $numTotalReads, "\n";
	print STATS "Total number of map reads	  ", $numMapped, "\n";
	print STATS "Total number of un-map reads	  ", $numUnMapped, "\n";
	print STATS "Total number of secundary mappings ", $numMultipleMappings, "\n\n";

	print STATS "Read Counts\n";
	print STATS "Total number of map reads in pairs		", $numMaPair, "\n";
	print STATS "Total number of map reads as broken pairs	", $numMaSing, "\n";
	print STATS "Total number of un-map reads in pairs		", $numUnPair, "\n";
	print STATS "Total number of un-map reads  as broken pairs	", $numUnSing, "\n\n";

	print STATS "Read Percents\n";
	print STATS "Percentage of map reads in pairs		", sprintf("%.2f", ($numMaPair/$numTotalReads)*100 ), "\n";
	print STATS "Percentage of map reads as broken pairs		", sprintf("%.2f", ($numMaSing/$numTotalReads)*100 ), "\n";
	print STATS "Percentage of un-map reads in pairs		", sprintf("%.2f", ($numUnPair/$numTotalReads)*100 ), "\n";
	print STATS "Percentage of un-map reads  as broken pairs	", sprintf("%.2f", ($numUnSing/$numTotalReads)*100 ), "\n\n";

} #end each argv

exit;

####################################
####################################
#
#     END OF SCRIPT
#
####################################
####################################

sub getBitFlag {

	my $bitwise = unpack("B32", pack("N", shift));

	my @bitString = split("", substr(reverse($bitwise),0,12));

	my @bitFlags = ('read_paired',
		'read_mapped_in_proper_pair',
		'read_unmapped',
		'mate_unmapped',
		'read_reverse_strand',
		'mate_reverse_strand',
		'first_in_pair',
		'second_in_pair',
		'not_primary_alignment',
		'read_fails_platform/vendor_quality_checks',
		'read_is_PCR_or_optical_duplicate',
		'supplementary_alignment'
	);

	my @upFlags;

	foreach my $i (0..11) { if($bitString[$i] == 1) {push(@upFlags, $bitFlags[$i])} }

	return(@upFlags);
	
} #end getBitFlag sub







































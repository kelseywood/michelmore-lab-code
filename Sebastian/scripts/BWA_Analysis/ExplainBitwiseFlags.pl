#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) {die "please give me a flag!\n"}

my $flag = $ARGV[0];

my @explanation = getBitFlag($flag);

print join( "\t", @explanation), "\n";

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


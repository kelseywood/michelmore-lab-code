#!/usr/bin/perl

if (@ARGV < 1) {

	print "Two input parameters are needed, please input them.\n";
	print "Usage: SAM_Ranges_Summarizer.pl <SAM file> <Minimmun Phred quality>\n";
	exit;

} #end if


my $phred = $ARGV[1];

open(SAM, $ARGV[0]);
open(RESULTS, ">$ARGV[0]_qualityAbove_$phred.txt");

$map_reads = 0;
$unmap_reads = 0;
$reads_above_phred = 0;
$reads_below_phred = 0;

while($line = <SAM>) {

	if ($line =~ /^@/) {next}

	chomp($line);

	@fields = split("\t", $line);

	if($fields[2] ne "*") {

		if($fields[4] >= $phred) {

			$end = $fields[3] + length($fields[9]);

			print RESULTS $fields[0], "\t", $fields[2], "\t", $fields[4], "\t", $fields[3], "\t", $end, "\n";

			++$reads_above_phred;
			++$map_reads;

		} else {

			++$reads_below_phred;
			++$map_reads;

		} #end else

	} else {

		++$unmap_reads;

	} #end else

} #end while

print "General Stats\n";
print "SAM file $ARGV[0]\n";
print "Minimmun phred quality $phred\n";
print "Total number of reads on SAM file: ", ($map_reads + $unmap_reads), "\n";
print "Total number of unmap reads: $unmap_reads\n";
print "Total number of map reads: $map_reads\n";
print "Reads above $phred quality: $reads_above_phred\n";
print "Reads below $phred quality: $reads_below_phred\n";

exit;





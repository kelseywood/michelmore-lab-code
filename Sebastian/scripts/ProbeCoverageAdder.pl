#!/usr/bin/perl

if (@ARGV < 2) {

	print "Two arguments needed, please input them.\n";
        print "First argument: Number of probes per position (The colums are Contig, Position and number of probes).\n";
        print "Second argument: Table to input add the data (The first two colums are Contig and Position,\n";
        print "Third argument: Number of header lines.\n";
        print "the rest of the colums will be printed has such.\n";
        exit 0;

} #end if

open(INPUT1, $ARGV[0]) || die "Unable to open the file $ARGV[0]"; # First input, file with number of probes per position

my %coverage;

while ($line = <INPUT1>) {

	chop($line);

	@POSInfo = split ("\t", $line);

        $posID = "$POSInfo[0].$POSInfo[1]";

	$coverage{$posID} = $POSInfo[2];

}# end while

close(INPUT1);

print "Done loading the positions probe coverage.\n";

open(INPUT2, $ARGV[1]) || die "Unable to open the file $ARGV[1]"; # Second input, table to add the coverage

open(RESULTS, ">$ARGV[1]_probeCoverage.txt");

#for ($header = 0; $header < $ARGV[2]; ++$header) {
#
#	$line2 = <INPUT2>;
#
#        chop($line2);
#
#       	@headerInfo = split ("\t", $line2);
#
#        print RESULTS $headerInfo[0], "\t", $headerInfo[1], "\tProbeCoverage";
#
#        for($field = 2, $field <= @headerInfo, ++$field) {
#
#        	print RESULTS "\t", $headerInfo[$field];
#
#        } #end for
#
#        print RESULTS "\n";
#
#} #end for

#print "Done printing the file header to the new file.\n";

while($line2 = <INPUT2>) {

	chop($line2);

	my @contigInfo = split ("\t", $line2);

        $posId = "$contigInfo[0].$contigInfo[1]";

        $coverage = $coverage{$posId};

        print RESULTS $contigInfo[0], "\t", $contigInfo[1], "\t", $coverage;

        for($field = 2; $field <= @contigInfo; ++$field) {

        	print RESULTS "\t", $contigInfo[$field];

        } #end for

        print RESULTS "\n";

} #end while

close(INPUT2);
close(RESULTS);

print "Finish printing the position data with the probe coverage to the new file.\n";

exit 0;


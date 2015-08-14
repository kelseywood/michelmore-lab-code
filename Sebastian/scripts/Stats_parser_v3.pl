#!/usr/bin/perl -w
use warnings;

my $Dir = $ARGV[0];
my @myFiles = ();

# open and read DIR
opendir(DIR, $Dir);

while ( defined( my $file = readdir(DIR) ) ) {

	# skip '.' and '..'
	if ( $file ne "." && $file ne ".." ) {

       		push( @myFiles, $file );

        }    # end if
}    # end while

closedir(DIR);

print "\n";
print join("\n", @myFiles);
print "\n";

open (RESULTS, ">Compilation_Stats_Files.txt");

print RESULTS "CompID\tTotalSNP\tTotalSNPContig\tNonDetectedSNPContig\tDetectedSNP\tNonDetectedSNP\tTotalSPP\tTotalSPPContig\tNonDetectedSPPContig\tFalseSPP\tCorrectSPP\n";

foreach $file (@myFiles) {

	open(FILE, "$Dir/$file");

	$line1 = <FILE>;
        chop($line1);

        @line1_2 = split(" ", $line1);
        $compID = pop(@line1_2);

        print RESULTS $compID, "\t";

	$line2 = <FILE>;

	$line3 = <FILE>;
	chop($line3);

        @line3_2 = split(" ", $line3);
        $totalSPPs = pop(@line3_2);

        print RESULTS $totalSPPs, "\t";

	$line4 = <FILE>;
 	chop($line4);

        @line4_2 = split(" ", $line4);
        $totalContigSPPs = pop(@line4_2);

        print RESULTS $totalContigSPPs, "\t";

	$line5 = <FILE>;
 	chop($line5);

        @line5_2 = split(" ", $line5);
        $missContigSPPs = pop(@line5_2);

        print RESULTS $missContigSPPs, "\t";

	$line6 = <FILE>;
 	chop($line6);

        @line6_2 = split(" ", $line6);
        $nonDetectedSPPs = pop(@line6_2);

        print RESULTS $nonDetectedSPPs, "\t";

	$line7 = <FILE>;
 	chop($line7);

        @line7_2 = split(" ", $line7);
        $detectedSPPs = pop(@line7_2);

        print RESULTS $detectedSPPs, "\t";

        print RESULTS "\n";

        close(FILE);

} #end foreach

close(RESULTS);

exit;
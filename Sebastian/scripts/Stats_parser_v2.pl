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
        chop($line2);

	$line3 = <FILE>;
	chop($line3);

        @line3_2 = split(" ", $line3);
        $totalSNPs = pop(@line3_2);

        print RESULTS $totalSNPs, "\t";

	$line4 = <FILE>;
 	chop($line4);

        @line4_2 = split(" ", $line4);
        $totalContigSNPs = pop(@line4_2);

        print RESULTS $totalContigSNPs, "\t";

	$line5 = <FILE>;
 	chop($line5);

        @line5_2 = split(" ", $line5);
        $missContigSNPs = pop(@line5_2);

        print RESULTS $missContigSNPs, "\t";

	$line6 = <FILE>;
 	chop($line6);

        @line6_2 = split(" ", $line6);
        $detectedSNPs = pop(@line6_2);

        print RESULTS $detectedSNPs, "\t";

	$line7 = <FILE>;
 	chop($line7);

        @line7_2 = split(" ", $line7);
        $nonDetectedSNPs = pop(@line7_2);

        print RESULTS $nonDetectedSNPs, "\t";


	$line8 = <FILE>;
	chop($line8);

	$line9 = <FILE>;
	chop($line9);

        @line9_2 = split(" ", $line9);
        $totalSPPs = pop(@line9_2);

        print RESULTS $totalSPPs, "\t";

	$line10 = <FILE>;
 	chop($line10);

        @line10_2 = split(" ", $line10);
        $totalContigSPPs = pop(@line10_2);

        print RESULTS $totalContigSPPs, "\t";

	$line11 = <FILE>;
 	chop($line11);

        @line11_2 = split(" ", $line11);
        $missContigSPPs = pop(@line11_2);

        print RESULTS $missContigSPPs, "\t";

	$line12 = <FILE>;
 	chop($line12);

        @line12_2 = split(" ", $line12);
        $detectedSPPs = pop(@line12_2);

        print RESULTS $detectedSPPs, "\t";

	$line13 = <FILE>;
 	chop($line13);

        @line13_2 = split(" ", $line13);
        $nonDetectedSPPs = pop(@line13_2);

        print RESULTS $nonDetectedSPPs, "\t";

        print RESULTS "\n";

        close(FILE);

} #end foreach

close(RESULTS);

exit;
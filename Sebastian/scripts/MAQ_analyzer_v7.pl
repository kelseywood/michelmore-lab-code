if(scalar(@ARGV) < 4) {

	print "Usage: Four input parameters are needed. Please enter them\n";
	print "First input is the file with the MAQ results combined for both genotypes\n";
	print "Second and Third input are the genotypes\n";
	print "Fourth input is the minimun number of reads per genotype to call a SNP\n";
	exit 0;

}

open(INPUT1, $ARGV[0]) || die "Unable to open file $ARGV[0]"; # First file with the results from the maq

$geno1 = $ARGV[1];
$geno2 = $ARGV[2];

$numReads = $ARGV[3];

open (RESULTS, ">Maq_results_$ARGV[0].txt") || die "Unable to open results file";

open (SNPRESULTS, ">Maq_results_SNPFile_$ARGV[0]_depth$numReads.txt") || die "Unable to open SNP's results file"; 

open (COVERSNPRESULTS, ">Maq_results_CoveredSNP_$ARGV[0]_depth$numReads.txt") || die "Unable to open SNP's results file";

print RESULTS "ContigSal\tPos\t#Probes\tRefSeq\tConsensus$geno1\tDepth$geno1\tConsensus$geno2\tDepth$geno2\tPosStats\n";

print SNPRESULTS "Contig\tPos\t#Probes\tRefSeq\tConsensus$geno1\tDepth$geno1\tA\tC\tG\tT\tConsensus$geno2\tDepth$geno2\tA\tC\tG\tT\n";

print COVERSNPRESULTS "Contig\tPos\t#Probes\tRefSeq\tConsensus$geno1\tDepth$geno1\tA\tC\tG\tT\tConsensus$geno2\tDepth$geno2\tA\tC\tG\tT\n";

my @depthGeno1 = ();
my @depthGeno2 = ();
my $Geno1NoCover = 0;
my $Geno2NoCover = 0;
my $nonValidGeno1 = 0;
my $nonValidGeno2 = 0;
my @contigSNP = ();
my $numSNP = 0;
my $numCoverSNP = 0;
my $numPos = 0;
my $validWOProbesGeno1 = 0;
my $validWOProbesGeno2 = 0;

while($line = <INPUT1>) {

	chop($line);
	@pos = split("\t", $line);

	$posDepthGeno1 = $pos[9];
        $posDepthGeno2 = $pos[19];

	if ($pos[9] > 0) {
		push (@depthGeno1, $pos[9]);
                if ($pos[9] < $numReads) {
                        ++$nonValidGeno1;
                } else {
			if($pos[3] < 1) {
				++$validWOProbesGeno1;
			} #end if
		} #end else
	} else {
		++$Geno1NoCover;
	}

        if ($pos[19] > 0) {
		push (@depthGeno2, $pos[19]);
                if ($pos[19] < $numReads) {
                        ++$nonValidGeno2;
                } else {
                        if($pos[3] < 1) {
                                ++$validWOProbesGeno2;
                        } #end if
                } #end else

        } else {
                ++$Geno2NoCover;
        }

	print RESULTS $pos[1],"\t", $pos[2],"\t", $pos[3],"\t", $pos[7],"\t", $pos[8],"\t", $pos[9],"\t", $pos[18],"\t", $pos[19], "\t";

	if($pos[5] eq $pos[15] && $pos[6] eq $pos[16]) {

		if(@pos[8] ne "N" && @pos[8] ne "" && @pos[18] ne "N" && @pos[18] ne "") {

			if(@pos[8] eq @pos[18]) {
				print RESULTS "Same";
			} else {

				if($pos[9] > $numReads && $pos[19] > $numReads) {

					print SNPRESULTS $pos[1],"\t", $pos[2],"\t", $pos[3],"\t", $pos[7],"\t", $pos[8],"\t", $pos[9],"\t", $pos[10],"\t", $pos[11],"\t", $pos[12],"\t", $pos[13],"\t", $pos[18],"\t", $pos[19], "\t", $pos[20], "\t", $pos[21], "\t", $pos[22], "\t", $pos[23];
					print SNPRESULTS "\n";

					++$numSNP;
					push(@contigSNP, $pos[1]);

	                                if($pos[3] > 1) {

 	                                        print RESULTS "CoverSNP";

		                                print COVERSNPRESULTS $pos[1],"\t", $pos[2],"\t", $pos[3],"\t", $pos[7],"\t", $pos[8],"\t", $pos[9],"\t", $pos[10],"\t", $pos[11],"\t", $pos[12],"\t", $pos[13],"\t", $pos[18],"\t", $pos[19], "\t", $pos[20], "\t", $pos[21], "\t", $pos[22], "\t", $pos[23];
						print COVERSNPRESULTS "\n";

                                       		++$numCoverSNP;

					} else {

                                                print RESULTS "nonCoverSNP";

					} #end else

				}  else {

                                        print RESULTS "lowDepthSNP";

				} #end else

			} #end else

		} elsif (@pos[8] eq "N" || @pos[18] eq "N") {

			print RESULTS "N";
		} else {

			print RESULTS "NA";

		} #end else


	} else {

		print RESULTS "Not Match"

	} #end else

	print RESULTS "\n";

	++$numPos;

} #end while

my $totalDepthGeno1 = 0;

foreach (@depthGeno1) {
	$totalDepthGeno1 += $_;
} 

$averDepthGeno1 = $totalDepthGeno1 / (scalar @depthGeno1);

my $totalDepthSer = 0;

foreach (@depthGeno2) {
        $totalDepthGeno2 += $_;
}

$averDepthGeno2 = $totalDepthGeno2 / (scalar @depthGeno2);

%contigSNPs = map {$_, 1} @contigSNP;

@uniqueContig = keys %contigSNPs;
$numUniqueContig = scalar (@uniqueContig);

print "\nAnalyzing file: $ARGV[0]\n";
print "Depth for alignment : $numReads\n\n";

print "Number of total positions on the file: $numPos\n";
print "Number of non cover positions for $geno1: $Geno1NoCover\n";
print "Number of non cover positions for $geno2: $Geno2NoCover\n";
print "Number of non valid positions for $geno1: $nonValidGeno1\n";
print "Number of non valid positions for $geno2: $nonValidGeno2\n";
print "Number of valid positions without probes for $geno1: $validWOProbesGeno1\n";
print "Number of valid positions without probes for $geno2: $validWOProbesGeno2\n";
print "Average depth for $geno1: $averDepthGeno1\n";
print "Average depth for $geno2: $averDepthGeno2\n\n";

print "Number of contigs with SNP's: $numUniqueContig\n";
print "Number of SNP's found: $numSNP\n";
print "Number of covered SNP's found: $numCoverSNP\n\n";

close(INPUT1);
close(RESULTS);
close(SNPRESULTS);
close(COVERSNPRESULTS);

end;



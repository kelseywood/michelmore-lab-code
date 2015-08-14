
if(@ARGV < 6) {

	print "Program usage: Four arguments needed, please input them.\n";
        print "First Genotype table file (from genotype table generator script).\n";
        print "Second and Third argument are the genotypes to compare (i.e 1 2).\n";
        print "Fourth argument is the minimmun depth per position (10).\n";
        print "Fifth argument is the minimmun allele frequency (heterozygote calling) (0.75).\n";
        print "Sixth argument is the output prefix name.\n";
        print "Seventh argument determine header presence (TRUE or FALSE).\n";
        exit 0;
}

my $minDepth = $ARGV[3];
my $alleleFreq = $ARGV[4];

###########################
#                         #
#   Calculate colums to   #
#        work with        #
#                         #
###########################

my $ConsGeno1 = 9 + (($ARGV[1]-1)*3);
my $colDepthGeno1 = 10 + (($ARGV[1]-1)*3);
my $AlleleCompGeno1 = 11 + (($ARGV[1]-1)*3);

my $ConsGeno2 = 9 + (($ARGV[2]-1)*3);
my $colDepthGeno2 = 10 + (($ARGV[2]-1)*3);
my $AlleleCompGeno2 = 11 + (($ARGV[2]-1)*3);


##################
#                #
#   Open files   #
#                #
##################

open(GENOTABLE,  $ARGV[0]);

open(NONCOVER, ">$ARGV[5]_minAlleleFreq_$alleleFreq.depth_$minDepth.noncover.txt");
open(SAME, ">$ARGV[5]_minAlleleFreq_$alleleFreq.depth_$minDepth.nonSNP.txt");
open(SNP, ">$ARGV[5]_minAlleleFreq_$alleleFreq.depth_$minDepth.SNP.txt");
open(HETS, ">$ARGV[5]_minAlleleFreq_$alleleFreq.depth_$minDepth.heterozygotes.txt");

if($ARGV[6] eq "TRUE") {
	$fileNamesline = <GENOTABLE>;
	$headerline = <GENOTABLE>;
}

while ($snpLine = <GENOTABLE>) {

        ###########################
        #                         #
        #   Get Data from line    #
	#                         #
        ###########################

	chop $snpLine;
        @snpData = split("\t", $snpLine);

	$depthGeno1 = $snpData[$colDepthGeno1];
	$depthGeno2 = $snpData[$colDepthGeno2];
	$alleleGeno1 = $snpData[$ConsGeno1];
	$alleleGeno2 = $snpData[$ConsGeno2];
	@alleleCompGeno1 = split(":", $snpData[$AlleleCompGeno1]);
	@alleleCompGeno2 = split(":", $snpData[$AlleleCompGeno2]);

        $refAllele = $snpData[2];

	$depth = $depthGeno1 + $depthGeno2;


        ###########################
        #                         #
        #    Calculate Allele     #
      	#        Frequency        #
	#                         #
        ###########################

        if ($depthGeno1 ne 0) {

	     	$freqAGeno1 = ($alleleCompGeno1[0]/$depthGeno1);
     		$freqCGeno1 = ($alleleCompGeno1[1]/$depthGeno1);
	     	$freqGGeno1 = ($alleleCompGeno1[2]/$depthGeno1);
     		$freqTGeno1 = ($alleleCompGeno1[3]/$depthGeno1);

        } else {
	        $freqAGeno1 = 0;
     		$freqCGeno1 = 0;
	     	$freqGGeno1 = 0;
     		$freqTGeno1 = 0;
        } #end else

        if ($depthGeno2 ne 0) {

	     	$freqAGeno2 = ($alleleCompGeno2[0]/$depthGeno2);
     		$freqCGeno2 = ($alleleCompGeno2[1]/$depthGeno2);
	     	$freqGGeno2 = ($alleleCompGeno2[2]/$depthGeno2);
     		$freqTGeno2 = ($alleleCompGeno2[3]/$depthGeno2);

        } else {
	        $freqAGeno2 = 0;
     		$freqCGeno2 = 0;
	     	$freqGGeno2 = 0;
     		$freqTGeno2 = 0;
        } #end else


        ###########################
        #                         #
        #    Compare Genotypes    #
	#                         #
        ###########################

        if(($alleleGeno1 && $alleleGeno2) ne "N") {

		if($depthGeno1 > $minDepth && $depthGeno2 > $minDepth) {

        		if ($freqAGeno1 > $alleleFreq || $freqCGeno1 > $alleleFreq || $freqGGeno1 > $alleleFreq || $freqTGeno1 > $alleleFreq) {

                		if ($freqAGeno2 > $alleleFreq || $freqCGeno2 > $alleleFreq || $freqGGeno2 > $alleleFreq || $freqTGeno2 > $alleleFreq) {

                        		if ($alleleGeno1 ne $alleleGeno2) {

                                		print SNP $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";

	                                } else {
        	                        	print SAME $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";
                	                } #end else

                        	} else {
                         		print HETS $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";
	                        } #end else

        	        } else {
                	      	print HETS $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";
	                } #end else

	        } else {
        	 	print NONCOVER $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";
	        } #end else

        } else {
        	print NONCOVER $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $depth, "\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[$AlleleCompGeno1], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[$AlleleCompGeno2], "\n";
        } #end else

} #end while

exit;
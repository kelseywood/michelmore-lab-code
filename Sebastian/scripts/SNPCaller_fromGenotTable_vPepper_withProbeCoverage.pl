
if(@ARGV < 5) {

	print "Program usage: Four arguments needed, please input them.\n";
        print "First Genotype table file (from genotype table generator script).\n";
        print "Second argument is the minimmun depth per position (10).\n";
        print "Third argument is the minimmun allele frequency (heterozygote calling) (0.75).\n";
        print "Fourth argument determine header presence (TRUE or FALSE).\n";
        print "Fifth argument is the output prefix name.\n";
        exit 0;
}

my $minDepth = $ARGV[1];
my $alleleFreq = $ARGV[2];

##################
#                #
#   Open files   #
#                #
##################

open(GENOTABLE,  $ARGV[0]);

open(NONCOVER, ">$ARGV[4]_minAlleleFreq_$alleleFreq.depth_$minDepth.noncover.txt");
open(SAME, ">$ARGV[4]_minAlleleFreq_$alleleFreq.depth_$minDepth.nonSNP.txt");
open(SNP, ">$ARGV[4]_minAlleleFreq_$alleleFreq.depth_$minDepth.SNP.txt");
open(HETS, ">$ARGV[4]_minAlleleFreq_$alleleFreq.depth_$minDepth.heterozygotes.txt");

if($ARGV[3] eq "TRUE") {
	$fileNamesline = <GENOTABLE>;
	$headerline = <GENOTABLE>;
}

print NONCOVER $fileNamesline;
print NONCOVER $headerline;

print SAME $fileNamesline;
print SAME $headerline;

print SNP $fileNamesline;
print SNP $headerline;

print HETS $fileNamesline;
print HETS $headerline;

while ($snpLine = <GENOTABLE>) {

        ###########################
        #                         #
        #   Get Data from line    #
	#                         #
        ###########################

	chop $snpLine;
        @snpData = split("\t", $snpLine);

	$alleleGeno1 = $snpData[10];
	$alleleGeno2 = $snpData[13];
	$alleleGeno3 = $snpData[16];
	@alleleCompGeno1 = split(":", $snpData[12]);
	@alleleCompGeno2 = split(":", $snpData[15]);
	@alleleCompGeno3 = split(":", $snpData[18]);

       	$depthGeno1 = (@alleleCompGeno1[0] + @alleleCompGeno1[1] + @alleleCompGeno1[2] + @alleleCompGeno1[3]);
	$depthGeno2 = (@alleleCompGeno2[0] + @alleleCompGeno2[1] + @alleleCompGeno2[2] + @alleleCompGeno2[3]);
	$depthGeno3 = (@alleleCompGeno3[0] + @alleleCompGeno3[1] + @alleleCompGeno3[2] + @alleleCompGeno3[3]);

        $refAllele = $snpData[3];

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

        if ($depthGeno3 ne 0) {

	     	$freqAGeno3 = ($alleleCompGeno3[0]/$depthGeno3);
     		$freqCGeno3 = ($alleleCompGeno3[1]/$depthGeno3);
	     	$freqGGeno3 = ($alleleCompGeno3[2]/$depthGeno3);
     		$freqTGeno3 = ($alleleCompGeno3[3]/$depthGeno3);

        } else {
	        $freqAGeno3 = 0;
     		$freqCGeno3 = 0;
	     	$freqGGeno3 = 0;
     		$freqTGeno3 = 0;
        } #end else

        if ($freqAGeno1 > $alleleFreq || $freqCGeno1 > $alleleFreq || $freqGGeno1 > $alleleFreq || $freqTGeno1 > $alleleFreq) {

        	$hetGeno1 = FALSE;

        } else {

	        $hetGeno1 = TRUE;

        } #end else

        if ($freqAGeno2 > $alleleFreq || $freqCGeno2 > $alleleFreq || $freqGGeno2 > $alleleFreq || $freqTGeno2 > $alleleFreq) {

        	$hetGeno2 = FALSE;

        } else {

	        $hetGeno2 = TRUE;

        } #end else

        if ($freqAGeno3 > $alleleFreq || $freqCGeno3 > $alleleFreq || $freqGGeno3 > $alleleFreq || $freqTGeno3 > $alleleFreq) {

        	$hetGeno2 = FALSE;

        } else {

	        $hetGeno2 = TRUE;

        } #end else


        ###########################
        #                         #
        #    Compare Genotypes    #
	#                         #
        ###########################

        if( (($alleleGeno1 && $alleleGeno2) ne "N") || (($alleleGeno1 && $alleleGeno3) ne "N") || (($alleleGeno3 && $alleleGeno2) ne "N")) {

		if(($depthGeno1 > $minDepth && $depthGeno2 > $minDepth) || ($depthGeno1 > $minDepth && $depthGeno3 > $minDepth) || ($depthGeno3 > $minDepth && $depthGeno2 > $minDepth)) {

                	if( (($hetGeno1 && $hetGeno2) eq "FALSE") || (($hetGeno1 && $hetGeno3) eq "FALSE") || (($hetGeno3 && $hetGeno2) eq "FALSE")) {

                   		if (($alleleGeno1 ne $alleleGeno2) || ($alleleGeno1 ne $alleleGeno3) || ($alleleGeno3 ne $alleleGeno2)) {

                                	print SNP $snpLine, "\n";

	                        } else {
        		          	print SAME $snpLine, "\n";
                		} #end else

                	} else {

                        	print HETS $snpLine, "\n";

                        } #end else

	        } else {
        	 	print NONCOVER $snpLine, "\n";
	        } #end else

        } else {
        	print NONCOVER $snpLine, "\n";
        } #end else

} #end while

exit;
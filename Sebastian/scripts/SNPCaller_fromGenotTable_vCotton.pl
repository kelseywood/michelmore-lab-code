
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

	$alleleGeno1 = $snpData[9];
	$alleleGeno2 = $snpData[12];
	$alleleGeno3 = $snpData[15];
	$alleleGeno4 = $snpData[18];
	$alleleGeno5 = $snpData[21];


	@alleleCompGeno1 = split(":", $snpData[11]);
	@alleleCompGeno2 = split(":", $snpData[14]);
	@alleleCompGeno3 = split(":", $snpData[17]);
	@alleleCompGeno4 = split(":", $snpData[20]);
	@alleleCompGeno5 = split(":", $snpData[23]);

       	$depthGeno1 = (@alleleCompGeno1[0] + @alleleCompGeno1[1] + @alleleCompGeno1[2] + @alleleCompGeno1[3]);
	$depthGeno2 = (@alleleCompGeno2[0] + @alleleCompGeno2[1] + @alleleCompGeno2[2] + @alleleCompGeno2[3]);
	$depthGeno3 = (@alleleCompGeno3[0] + @alleleCompGeno3[1] + @alleleCompGeno3[2] + @alleleCompGeno3[3]);
	$depthGeno4 = (@alleleCompGeno4[0] + @alleleCompGeno4[1] + @alleleCompGeno4[2] + @alleleCompGeno4[3]);
	$depthGeno5 = (@alleleCompGeno5[0] + @alleleCompGeno5[1] + @alleleCompGeno5[2] + @alleleCompGeno5[3]);

        $refAllele = $snpData[2];

	$depth = $depthGeno1 + $depthGeno2 + $depthGeno3 + $depthGeno4 + $depthGeno5;


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

        if ($depthGeno4 ne 0) {

	     	$freqAGeno4 = ($alleleCompGeno4[0]/$depthGeno4);
     		$freqCGeno4 = ($alleleCompGeno4[1]/$depthGeno4);
	     	$freqGGeno4 = ($alleleCompGeno4[2]/$depthGeno4);
     		$freqTGeno4 = ($alleleCompGeno4[3]/$depthGeno4);

        } else {
	        $freqAGeno4 = 0;
     		$freqCGeno4 = 0;
	     	$freqGGeno4 = 0;
     		$freqTGeno4 = 0;
        } #end else

        if ($depthGeno5 ne 0) {

	     	$freqAGeno5 = ($alleleCompGeno5[0]/$depthGeno5);
     		$freqCGeno5 = ($alleleCompGeno5[1]/$depthGeno5);
	     	$freqGGeno5 = ($alleleCompGeno5[2]/$depthGeno5);
     		$freqTGeno5 = ($alleleCompGeno5[3]/$depthGeno5);

        } else {
	        $freqAGeno5 = 0;
     		$freqCGeno5 = 0;
	     	$freqGGeno5 = 0;
     		$freqTGeno5 = 0;
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

        	$hetGeno3 = FALSE;

        } else {

	        $hetGeno3 = TRUE;

        } #end else

         if ($freqAGeno4 > $alleleFreq || $freqCGeno4 > $alleleFreq || $freqGGeno4 > $alleleFreq || $freqTGeno4 > $alleleFreq) {

        	$hetGeno4 = FALSE;

        } else {

	        $hetGeno4 = TRUE;

        } #end else

         if ($freqAGeno5 > $alleleFreq || $freqCGeno5 > $alleleFreq || $freqGGeno5 > $alleleFreq || $freqTGeno5 > $alleleFreq) {

        	$hetGeno5 = FALSE;

        } else {

	        $hetGeno5 = TRUE;

        } #end else


        ###########################
        #                         #
        #    Compare Genotypes    #
	#                         #
        ###########################

        my $nonNcount = 0;
        if($alleleGeno1 ne "N") { ++$nonNcount;}
        if($alleleGeno2 ne "N") { ++$nonNcount;}
        if($alleleGeno3 ne "N") { ++$nonNcount;}
        if($alleleGeno4 ne "N") { ++$nonNcount;}
        if($alleleGeno5 ne "N") { ++$nonNcount;}

        my $covCount = 0;
        if($depthGeno1 > $minDepth) { ++$covCount;}
        if($depthGeno2 > $minDepth) { ++$covCount;}
        if($depthGeno3 > $minDepth) { ++$covCount;}
        if($depthGeno4 > $minDepth) { ++$covCount;}
        if($depthGeno5 > $minDepth) { ++$covCount;}

	my $hetCount = 0;
        if($hetGeno1 eq "TRUE") { ++$hetCount;}
        if($hetGeno2 eq "TRUE") { ++$hetCount;}
        if($hetGeno3 eq "TRUE") { ++$hetCount;}
        if($hetGeno4 eq "TRUE") { ++$hetCount;}
        if($hetGeno5 eq "TRUE") { ++$hetCount;}

        my %alleles = ();

        $alleles{$alleleGeno1} = "1";
        $alleles{$alleleGeno2} = "2";
        $alleles{$alleleGeno3} = "3";
        $alleles{$alleleGeno4} = "4";
        $alleles{$alleleGeno5} = "5";

        @Allele = keys %alleles, "\n";

        $numAlleles = scalar @Allele;

        if($nonNcount >= 2) {

		if($covCount >= 2) {

                	if($hetCount <= 3) {

                        	if ($numAlleles > 1) {

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
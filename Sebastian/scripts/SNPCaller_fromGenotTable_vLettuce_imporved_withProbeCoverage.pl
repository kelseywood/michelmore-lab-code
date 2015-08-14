#! /usr/bin/perl


if(@ARGV < 6) {

	print "Program usage: Five arguments needed, please input them.\n";
        print "First Genotype table file (from genotype table generator script).\n";
        print "Second argument is the minimmun depth per position (10).\n";
        print "Third argument is the maximum depth per position (although not being used in the algorithm you need to put a number here ).\n";
        print "Fourth argument is the minimmun allele frequency (heterozygote calling) (0.75).\n";
        print "Fifth argument determine header presence (TRUE or FALSE).\n";
        print "Sixth argument is the output prefix name.\n";
        exit 0;
}

$minDepth = $ARGV[1];
$maxDepth = $ARGV[2];
$alleleFreq = $ARGV[3];

##################
#                #
#   Open files   #
#                #
##################

open(GENOTABLE,  $ARGV[0]);

open(NONCOVER, ">$ARGV[5]_minAlleleFreq_$alleleFreq.mindepth_$minDepth.Non_or_High_Cover.txt");
open(SAME, ">$ARGV[5]_minAlleleFreq_$alleleFreq.mindepth_$minDepth.NonSNP.txt");
open(SNPWITHPROBES, ">$ARGV[5]_minAlleleFreq_$alleleFreq.mindepth_$minDepth.SNPwithProbes.txt");
open(SNPWITHOUTPROBES, ">$ARGV[5]_minAlleleFreq_$alleleFreq.mindepth_$minDepth.SNPwithoutProbes.txt");
open(HETS, ">$ARGV[5]_minAlleleFreq_$alleleFreq.mindepth_$minDepth.Heterozygotes.txt");

if($ARGV[4] eq "TRUE") {
	$fileNamesline = <GENOTABLE>;
	$headerline = <GENOTABLE>;
}

print NONCOVER $fileNamesline;
print NONCOVER join ("\t", ('Contig', 'Pos', 'NumProbes', 'Depth', 'Ref', 'Consensus', 'Depth', 'Consensus', 'Depth', 'A:C:G:T', 'Consensus', 'Depth', 'A:C:G:T'));
print NONCOVER "\n";

print SAME $fileNamesline;
print SAME join ("\t", ('Contig', 'Pos', 'NumProbes', 'Depth', 'Ref', 'Consensus', 'Depth', 'Consensus', 'Depth', 'A:C:G:T', 'Consensus', 'Depth', 'A:C:G:T'));
print SAME "\n";

print SNPWITHPROBES $fileNamesline;
print SNPWITHPROBES join ("\t", ('Contig', 'Pos', 'NumProbes', 'Depth', 'Ref', 'Consensus', 'Depth', 'Consensus', 'Depth', 'A:C:G:T', 'Consensus', 'Depth', 'A:C:G:T'));
print SNPWITHPROBES "\n";

print SNPWITHOUTPROBES $fileNamesline;
print SNPWITHOUTPROBES join ("\t", ('Contig', 'Pos', 'NumProbes', 'Depth', 'Ref', 'Consensus', 'Depth', 'Consensus', 'Depth', 'A:C:G:T', 'Consensus', 'Depth', 'A:C:G:T'));
print SNPWITHOUTPROBES "\n";

print HETS $fileNamesline;
print HETS join ("\t", ('Contig', 'Pos', 'NumProbes', 'Depth', 'Ref', 'Consensus', 'Depth', 'Consensus', 'Depth', 'A:C:G:T', 'Consensus', 'Depth', 'A:C:G:T'));
print HETS "\n";

while ($snpLine = <GENOTABLE>) {

        ###########################
        #                         #
        #   Get Data from line    #
	#                         #
        ###########################

	chop $snpLine;
        @snpData = split("\t", $snpLine);

        $numProbes = $snpData[2];
        $allele_ref = $snpData[3];
        $allele_consensus = $snpData[4];
	#read the allele composition for each line
	@alleleCompGeno1 = split(":", $snpData[12]);
	@alleleCompGeno2 = split(":", $snpData[15]);

       	$depthGeno1 = (@alleleCompGeno1[0] + @alleleCompGeno1[1] + @alleleCompGeno1[2] + @alleleCompGeno1[3]);
	$depthGeno2 = (@alleleCompGeno2[0] + @alleleCompGeno2[1] + @alleleCompGeno2[2] + @alleleCompGeno2[3]);

        #make an array of all depths
        @all_depth = ("$depthGeno1", "$depthGeno2");

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

        if ($freqAGeno1 > $alleleFreq) {
        	$hetGeno1 = FALSE;
		$hetGenoCover1 = FALSE;
                $alleleGeno1 = "A";

        } elsif ($freqCGeno1 > $alleleFreq) {
        	$hetGeno1 = FALSE;
		$hetGenoCover1 = FALSE;
                $alleleGeno1 = "C";

        } elsif ($freqGGeno1 > $alleleFreq) {
        	$hetGeno1 = FALSE;
		$hetGenoCover1 = FALSE;
                $alleleGeno1 = "G";

        } elsif ($freqTGeno1 > $alleleFreq) {
        	$hetGeno1 = FALSE;
		$hetGenoCover1 = FALSE;
                $alleleGeno1 = "T";

        } else {
        	$hetGeno1 = TRUE;
	        $alleleGeno1 = $snpData[10];

        	if($depthGeno1 > $minDepth ) {
		        $hetGenoCover1 = TRUE;
                } else {
                	$hetGenoCover1 = FALSE;
                } #end else
        } #end else

        if ($freqAGeno2 > $alleleFreq) {
        	$hetGeno2 = FALSE;
		$hetGenoCover2 = FALSE;
                $alleleGeno2 = "A";

        } elsif ($freqCGeno2 > $alleleFreq) {
        	$hetGeno2 = FALSE;
		$hetGenoCover2 = FALSE;
                $alleleGeno2 = "C";

        } elsif ($freqGGeno2 > $alleleFreq) {
        	$hetGeno2 = FALSE;
		$hetGenoCover2 = FALSE;
                $alleleGeno2 = "G";

        } elsif ($freqTGeno2 > $alleleFreq) {
        	$hetGeno2 = FALSE;
		$hetGenoCover2 = FALSE;
                $alleleGeno2 = "T";

        } else {
        	$hetGeno2 = TRUE;
		$alleleGeno2 = $snpData[13];

        	if($depthGeno2 > $minDepth ) {
		        $hetGenoCover2 = TRUE;
                } else {
                	$hetGenoCover2 = FALSE;
                } #end else
        } #end else

        #make an array of all consensus letters (A,C,G,T, R, w, N ...)
        @all_allele = ("$alleleGeno1", "$alleleGeno2");
        #print "@all_allele\n";
        $numberof_A = 0 ;
        $numberof_C = 0 ;
        $numberof_G = 0 ;
        $numberof_T = 0 ;
        $numberof_N = 0 ;
        $other_letters = 0;

        $curGeno = 0;

        #declare covered flags

         foreach $elem (@all_allele) {
                if ($elem eq "A") {
                    ++$numberof_A;
                    }
                 elsif($elem eq "C") {
                    ++$numberof_C;
                    }
                  elsif($elem eq "G") {
                    ++$numberof_G;
                    }
                   elsif($elem eq "T") {
                    ++$numberof_T;
                    }
                   elsif($elem eq "N") {
                    ++$numberof_N;
                    }
                  else {
                    ++$other_letters;
                    }

                  ++$curGeno;

        }

#        print join("\t", $snpData[0], $snpData[1], $snpData[2], $snpData[3], $consensusTrueCover, $polimorphicTrueCover);
#	print "\n";

#         print join ("\t", $numberof_A, $numberof_C, $numberof_G, $numberof_T, $numberof_N, $other_letters);
#         print "\n";
         #initiate the true SNP flag;
         $trueSNP = 0 ;

        if ($numberof_A  >= 1){
           if ($numberof_C !=0 || $numberof_G !=0 || $numberof_T != 0){

              $trueSNP = 1;
           }
         }
         elsif ($numberof_C  >= 1){
           if ($numberof_A !=0 || $numberof_G !=0 || $numberof_T != 0){

              $trueSNP = 1;
           }
         }
         elsif ($numberof_G >= 1){
           if ($numberof_A !=0 || $numberof_C !=0 || $numberof_T != 0){

              $trueSNP = 1;
           }
         }
         elsif ($numberof_T >= 1){
           if ($numberof_A !=0 || $numberof_C !=0 || $numberof_C != 0){

              $trueSNP = 1;
           }
         }
         else {
              $trueSNP = 0 ;

         }

        ###########################
        #                         #
        #  Count N's and Hets's   #
	#                         #
        ###########################

        $coverCount = 0;
        if($depthGeno1 >= $minDepth) { ++$coverCount;}
        if($depthGeno2 >= $minDepth) { ++$coverCount;}


        $nonNcount = 0;
        if($alleleGeno1 ne "N") { ++$nonNcount;}
        if($alleleGeno2 ne "N") { ++$nonNcount;}

	$hetCount = 0;
        if($hetGeno1 eq "TRUE") { ++$hetCount;}
        if($hetGeno2 eq "TRUE") { ++$hetCount;}

       	$hetCoverCount = 0;
        if($hetGenoCover1 eq "TRUE") { ++$hetCoverCount;}
        if($hetGenoCover2 eq "TRUE") { ++$hetCoverCount;}

        if($nonNcount == 2 && $coverCount == 2) {

               	if($hetCount == 0) {

                       	if ($trueSNP == 1) {

				if ($numProbes >= 3) {

                                        	print SNPWITHPROBES $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

                                       } else {

                                         	print SNPWITHOUTPROBES $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

                                       } #end else

        	               	} else {

              	                         print SAME $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

                                } #end else

                } else {

	                if ($hetCoverCount == 2) {

        	                print HETS $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

                	} else {

                        	print NONCOVER $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

	                } #end else

                } #end else

        } else {

        	print NONCOVER $snpData[0], "\t", $snpData[1], "\t", $numProbes, "\t", $refAllele, "\t", $allele_consensus, "\t", $depth,"\t", $alleleGeno1, "\t", $depthGeno1, "\t", $snpData[12], "\t", $alleleGeno2, "\t", $depthGeno2, "\t", $snpData[15], "\n";

        } #end else

} #end while

close(NONCOVER);
close(SAME);
close(SNPWITHPROBES);
close(SNPWITHOUTPROBES);
close(HETS);

exit;
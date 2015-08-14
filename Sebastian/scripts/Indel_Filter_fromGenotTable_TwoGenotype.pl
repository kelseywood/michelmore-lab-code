
if(@ARGV < 5) {

	print "Program usage: Four arguments needed, please input them.\n";
        print "First Genotype table file (from genotype table generator script).\n";
        print "Second and Third argument are the genotypes to compare (i.e 1 2).\n";
        print "Fourth argument is the minimmun depth per position (10).\n";
        print "Fifth argument is the output prefix name.\n";
        print "Sixth argument determine header presence (TRUE or FALSE).\n";
        exit 0;
}

my $minDepth = $ARGV[3];

###########################
#                         #
#   Calculate colums to   #
#        work with        #
#                         #
###########################

my $colFirstAlleleGeno1 = 8 + (($ARGV[1]-1)*4);
my $colDepthFirstAlleleGeno1 = 9 + (($ARGV[1]-1)*4);
my $colSecondAlleleGeno1 = 10 + (($ARGV[1]-1)*4);
my $colDepthSecondAlleleGeno1 = 11 + (($ARGV[1]-1)*4);

my $colFirstAlleleGeno2 = 8 + (($ARGV[2]-1)*4);
my $colDepthFirstAlleleGeno2 = 9 + (($ARGV[2]-1)*4);
my $colSecondAlleleGeno2 = 10 + (($ARGV[2]-1)*4);
my $colDepthSecondAlleleGeno2 = 11 + (($ARGV[2]-1)*4);



##################
#                #
#   Open files   #
#                #
##################

open(GENOTABLE,  $ARGV[0]);

open(VALID, ">$ARGV[4]_depth_$minDepth.validIndels.txt");
open(NONVALID, ">$ARGV[4]_depth_$minDepth.nonValidIndels.txt");

if($ARGV[5] eq "TRUE") {
	$fileNamesline = <GENOTABLE>;
	$headerline = <GENOTABLE>;
}

while ($indelLine = <GENOTABLE>) {

        ###########################
        #                         #
        #   Get Data from line    #
	#                         #
        ###########################

	chop $indelLine;
        @snpData = split("\t", $indelLine);


	$firstAlleleGeno1 = $snpData[$colFirstAlleleGeno1];
	$depthFirstAlleleGeno1 = $snpData[$colDepthFirstAlleleGeno1];
	$secondAlleleGeno1 = $snpData[$colSecondAlleleGeno1];
	$depthSecondAlleleGeno1 = $snpData[$colDepthSecondAlleleGeno1];

	$firstAlleleGeno2 = $snpData[$colFirstAlleleGeno2];
	$depthFirstAlleleGeno2 = $snpData[$colDepthFirstAlleleGeno2];
	$secondAlleleGeno2 = $snpData[$colSecondAlleleGeno2];
	$depthSecondAlleleGeno2 = $snpData[$colDepthSecondAlleleGeno2];

        ###########################
        #                         #
        #    Filter Indel Line    #
	#                         #
        ###########################

        if (($firstAlleleGeno1 ne "*" && $depthFirstAlleleGeno1 > $minDepth) || ($secondAlleleGeno1 ne "*" && $depthSecondAlleleGeno1 > $minDepth) || ($firstAlleleGeno2 ne "*" && $depthFirstAlleleGeno2 > $minDepth) || ($secondAlleleGeno2 ne "*" && $depthSecondAlleleGeno2 > $minDepth)) {

        	print VALID  $snpData[0], "\t", $snpData[1], "\t", $firstAlleleGeno1, "\t", $depthFirstAlleleGeno1, "\t", $secondAlleleGeno1, "\t", $depthSecondAlleleGeno1, "\t", $firstAlleleGeno2, "\t", $depthFirstAlleleGeno2, "\t", $secondAlleleGeno2, "\t", $depthSecondAlleleGeno2, "\n";

        } else {

        	print NONVALID  $snpData[0], "\t", $snpData[1], "\t", $firstAlleleGeno1, "\t", $depthFirstAlleleGeno1, "\t", $secondAlleleGeno1, "\t", $depthSecondAlleleGeno1, "\t", $firstAlleleGeno2, "\t", $depthFirstAlleleGeno2, "\t", $secondAlleleGeno2, "\t", $depthSecondAlleleGeno2, "\n";

        } #end else

} #end while
























exit;
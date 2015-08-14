

#######################################################################
#                                                                     #
#    This script is designed to filter the Indels output it           #
#    from BWA_AlleleComposition_Calculator base on a specific         #
#    depth.                                                           #
#                                                                     #
#######################################################################

if(@ARGV < 1) {

	print "Program usage: Two arguments needed, please input them.\n";
        print "First indel file (from BWA_AlleleComposition_Calculator script).\n";
        print "Second argument is the minimmun depth per position (10).\n";
        exit 0;
}

my $minDepth = $ARGV[1];

##################
#                #
#   Open files   #
#                #
##################

open(GENOTABLE,  $ARGV[0]);

open(VALID, ">$ARGV[0]_depth_$minDepth.validIndels.txt");
open(NONVALID, ">$ARGV[0]_depth_$minDepth.nonValidIndels.txt");

$headerline = <GENOTABLE>;

print VALID $headerline;
print NONVALID $headerline;

while ($indelLine = <GENOTABLE>) {

        ###########################
        #                         #
        #   Get Data from line    #
	#                         #
        ###########################

	chop $indelLine;
        @snpData = split("\t", $indelLine);


	$firstAlleleGeno1 = $snpData[4];
	$depthFirstAlleleGeno1 = $snpData[5];
	$secondAlleleGeno1 = $snpData[6];
	$depthSecondAlleleGeno1 = $snpData[7];

        ###########################
        #                         #
        #    Filter Indel Line    #
	#                         #
        ###########################

        if (($firstAlleleGeno1 ne "*" && $depthFirstAlleleGeno1 > $minDepth) || ($secondAlleleGeno1 ne "*" && $depthSecondAlleleGeno1 > $minDepth)) {

        	print VALID  $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $snpData[3], "\t", $firstAlleleGeno1, "\t", $depthFirstAlleleGeno1, "\t", $secondAlleleGeno1, "\t", $depthSecondAlleleGeno1, "\n";


        } else {

        	print NONVALID  $snpData[0], "\t", $snpData[1], "\t", $snpData[2], "\t", $snpData[3], "\t", $firstAlleleGeno1, "\t", $depthFirstAlleleGeno1, "\t", $secondAlleleGeno1, "\t", $depthSecondAlleleGeno1, "\n";

        } #end else

} #end while

exit;
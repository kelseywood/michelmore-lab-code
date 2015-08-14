######################################################################
######################################################################
#                                                                    #
#                                                                    #
#             Script design to generate GFF3 files for               #
#                          Probes PIP's                              #
#                                                                    #
#                    Sebastian Reyes Chin-Wo                         #
#                   Seed Biotechnology Center                        #
#                University of California at Davis                   #
#                      sreyesch@ucdavis.edu                          #
#                                                                    #
#                                                                    #
######################################################################
######################################################################

#

if (@ARGV < 1) {

	print "Program usage:\n";
	print "Two input arguments needed, please enter them.\n";
	print "- First argument: Coordinates file with the positions and LG of all the contigs on the sequence.\n";
	print "- Second argument: PIP positions.\n";
	die "Input arguments missing";

} #end if


################################
#                              #
# Global and storage variables #
#                              #
################################

my %blocks = ();


#############################
#############################
#                           #
#   Open Coordinates File   #
#                           #
#############################
#############################


# Extract the sequences position, LG and
# loaded into the sequences hash

open (COORDINATES, $ARGV[0]);

while ($line = <COORDINATES>) {

	chop($line);

	@coordinate = split("\t", $line);

	$blocks{$coordinate[0]} = [$coordinate[2], $coordinate[3]];

} #end while

close (COORDINATES);

####################################
####################################
#                                  #
#          Open PIP File           #
#                                  #
####################################
####################################

open (PROBEPIP, $ARGV[1]);
open (GFF3FILE, ">$ARGV[1].gff3");

while($PIPLine = <PROBEPIP>) {

	chop($PIPLine);

	@PIPInfo = split("\t", $PIPLine);

	if($blocks{$PIPInfo[1]}[0] ne "") {

		$startPos = ($PIPInfo[2] + $blocks{$PIPInfo[1]}[1] - 13);

		$endPos = ($PIPInfo[2] + $blocks{$PIPInfo[1]}[1] + 11);

		print GFF3FILE "LG", $blocks{$PIPInfo[1]}[0],"\t.\tAffyProbe\t", $startPos, "\t", $endPos,"\t.\t.\t.\tName=", $PIPInfo[0], ";Contig=", $PIPInfo[1], "\t\t";
		print GFF3FILE "\n";

	} #end if

} #end while

close (PROBEPIP);
close (GFF3FILE);

end;




















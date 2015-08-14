######################################################################
######################################################################
#                                                                    #
#                                                                    #
#             Script design to generate complete Map                 #
#                                                                    #
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
	print "- First argument: File with the genetic bins (Result from MadMapper), it should have two columns,\n";
	print " first column with the MadMapper group ID and second with the marker code.\n";
	print "- Second argument: File with the map information.\n";
	die "Input arguments missing";

} #end if

#Create global storage variables

my %bins = ();
my %classes = ();


#################
#               #
#   Open Bin    #
#     file      #
#               #
#################

#Load bin information into storage hashes

open (BINFILE, $ARGV[0]);

while ($binline = <BINFILE>) {

	chop($binline);

	@markerInfo = split("\t", $binline);

	push (@{$bins{$markerInfo[0]}}, $markerInfo[1]);

} #end while

close (BINFILE);

print "Finish loading the bin data\n";


#################
#               #
#   Open Map    #
#     file      #
#               #
#################

#Loop through the map file and print new map with bins

open (MAPFILE, $ARGV[1]);
open (RESULTS, ">$ARGV[1].withBins.txt");

while ($mapline = <MAPFILE>) {

	chop($mapline);

	@posInfo = split("\t", $mapline);

	foreach $contig (@{$bins{$posInfo[0]}}) {

		print RESULTS $posInfo[0], "\t", $contig, "\t", $posInfo[1], "\t", $posInfo[2], "\n";

	} #end foreach

	@bins{$curClass} = ();

} #end while

print "Done printing the new map\n\n";

close (MAPFILE);
close (RESULTS);
######################################################################
######################################################################
#                                                                    #
#                                                                    #
#             Script design to generate GFF3 files for               #
#                      Blast Results PIP's                           #
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

if (@ARGV < 5) {

	print "Program usage:\n";
	print "Six input arguments needed, please enter them.\n";
	print "- First argument: Info2 file result from tcl_blast_parser.\n";
	print "- Second argument: Minimun lenght of the alignment (100).\n";
	print "- Third argument: Minimun percent of identity (90).\n";
	print "- Forth argument: Number of hits to export (1).\n";
	print "- Fifth argument: Class of the Features.\n";
	die "Input arguments missing";

} #end if


################################
#                              #
# Global and storage variables #
#                              #
################################

my %blocks = ();
$minLength = $ARGV[1];
$minIdentity = $ARGV[2];
$hits = $ARGV[3];

print "Input Arguments\n";
print "Blast result file: ", $ARGV[0], "\n";
print "Minimun length to be exported: ", $minLength, "\n";
print "Minimun percent of identity to be exported: ", $minIdentity, "\n";
print "Number of hits to export: ", $hits, "\n";


####################################
####################################
#                                  #
#      Open BLAST RESULT File      #
#                                  #
####################################
####################################

open (INFO2, $ARGV[0]);
open (GFF3FILE, ">$ARGV[0].gff3");

while($blastLine = <INFO2>) {

	chomp($blastLine);

	@blastInfo = split("\t", $blastLine);

	if($blastInfo[2] ne "no_hits_found") {

		if($blastInfo[6] > $minLength && $blastInfo[4] gt $minIdentity && $blastInfo[7] le $hits) {

			@strand = split("/", $blastInfo[8]);

			if($strand[0] eq "+") {

			print GFF3FILE $blastInfo[1],"\t.\t", $ARGV[4], "\t", $blastInfo[9], "\t", $blastInfo[10],"\t", $blastInfo[3],"\t", $strand[1], "\t.\tID=", $blastInfo[0], ";Name=", $blastInfo[0], ";Length Alignment=", $blastInfo[6], ";Percent Identity=", $blastInfo[4], ";Hit Number=", $blastInfo[7];

			} else {

			print GFF3FILE $blastInfo[1],"\t.\t", $ARGV[4], "\t", $blastInfo[10],"\t", $blastInfo[9], "\t", $blastInfo[3],"\t", $strand[1], "\t.\tID=", $blastInfo[0], ";Name=", $blastInfo[0], ";Length Alignment=", $blastInfo[6], ";Percent Identity=", $blastInfo[4], ";Hit Number=", $blastInfo[7];

			} #end else

			print GFF3FILE "\n";

		} #end if

	} #end if

} #end while

close (PROBEPIP);
close (GFF3FILE);

end;























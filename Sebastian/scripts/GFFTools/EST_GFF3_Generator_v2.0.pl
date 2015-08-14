######################################################################
######################################################################
#                                                                    #
#                                                                    #
#             Script design to generate GFF3 files for               #
#                         polymorphism's                             #
#                                                                    #
#                    Sebastian Reyes Chin-Wo                         #
#                   Seed Biotechnology Center                        #
#                University of California at Davis                   #
#                      sreyesch@ucdavis.edu                          #
#                                                                    #
#                                                                    #
######################################################################
######################################################################

################## UPDATES FROM v1.0 ###################################
#Added funcionality to work with multiple positions per contig (duplicated contigs.
#
######################################################################

if (@ARGV < 2) {

	print "Program usage:\n";
	print "Three input arguments needed, please enter them.\n";
	print "- First argument: Coordinates file with the positions and LG of all the contigs on the sequence.\n";
	print "- Second argument: EST information .\n";
	print "- Third argument: Source of the EST.\n";
	die "Input arguments missing";

} #end if


################################
#                              #
# Global and storage variables #
#                              #
################################

my $source = $ARGV[2];
my %blocks = ();


#############################
#############################
#                           #
#   Open Coordinates dile   #
#                           #
#############################
#############################


# Extract the sequences position, LG and
# loaded into the sequences hash

open (COORDINATES, $ARGV[0]);

while ($line = <COORDINATES>) {

	chop($line);

	@coordinate = split("\t", $line);

	push (@{$blocks{$coordinate[0]}}, "$coordinate[2]-$coordinate[3]");

} #end while

close (COORDINATES);

####################################
####################################
#                                  #
#          Open EST File           #
#                                  #
####################################
####################################

open (CONTIGESTFILE, $ARGV[1]);
open (GFF3FILE, ">$ARGV[1].gff3");

while($contigLine = <CONTIGESTFILE>) {

	chop($contigLine);

	@contigInfo = split("\t", $contigLine);

	if($blocks{$contigInfo[0]}[0] ne "") {

                foreach $coord (@{$blocks{$contigInfo[0]}}) {

                	@cordinateSet = split("-", $coord);

	                $contigInfo[2] =~ s/\(//g;
	                $contigInfo[2] =~ s/\)//g;

	                $contigInfo[2] =~ s/\|/\t/g;
	                @ESTsInfo = split("\t", $contigInfo[2]);

	                foreach $EST (@ESTsInfo) {

	                        $EST =~ s/\+/\t/g;
	                        @ESTdata = split("\t", $EST);

	                        $ESTdata[1] =~ s/\,/\t/g;
	                        @ESTpos = split("\t", $ESTdata[1]);

	                        $startPos = $ESTpos[0] + $cordinateSet[1];

	                        --$startPos;

	                        $endPos = $ESTpos[1] + $cordinateSet[1];

	                        --$endPos;

	                        print GFF3FILE "LG", $cordinateSet[0],"\t", $source,"\tEST\t", $startPos, "\t", $endPos,"\t.\t.\t.\tName=", $ESTdata[0],";Contig=", $contigInfo[0];
	                        print GFF3FILE "\n";

	                } #end foreach

        	} #end foreach

	} #end if

} #end while

close (POLIMORPHISMFILE);
close (GFF3FILE);

end;
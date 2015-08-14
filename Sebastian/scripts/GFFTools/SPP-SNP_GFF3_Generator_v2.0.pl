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

################## UPDATES FROM v1 ###################################
#Added funcionality to work with multiple positions per contig (duplicated contigs.
#
######################################################################


if (@ARGV < 4) {

	print "Program usage:\n";
	print "Five input arguments needed, please enter them.\n";
	print "- First argument: Coordinates file with the positions and LG of all the contigs on the sequence.\n";
	print "- Second argument: File with the SNP positions.\n";
	print "- Third argument: Kind of polimorphism, type SNP or SPP.\n";
	print "- Forth argument: Type of identifier of the SNP, type 1 for Contig name (CLS_S3_Contig1) or 2 for Code name (AAAA).\n";
	print "- Fifth argument: Source of the polymorphism.\n";
	die "Input arguments missing";

} #end if


################################
#                              #
# Global and storage variables #
#                              #
################################

my $type = $ARGV[2];
my $identifier = $ARGV[3];
my $source = $ARGV[4];
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

	if($identifier eq 1) {

		push (@{$blocks{$coordinate[0]}}, "$coordinate[2]-$coordinate[3]");

	} elsif ($identifier eq 2) {

		push (@{$blocks{$coordinate[1]}}, "$coordinate[2]-$coordinate[3]");

	} else {

		print "Couldn't identify the identifier, please input 1 for Contig name (CLS_S3_Contig1) or 2 for Code name (AAAA).\n";
		die "Wrong input identifier";

	}

} #end while

close (COORDINATES);

####################################
####################################
#                                  #
#          Open SNP File           #
#                                  #
####################################
####################################

open (POLIMORPHISMFILE, $ARGV[1]);
open (GFF3FILE, ">$ARGV[1].$ARGV[2].gff3");

while($poliLine = <POLIMORPHISMFILE>) {

	chop($poliLine);

	@poliInfo = split("\t", $poliLine);

	if($blocks{$poliInfo[0]}[0] ne "") {

        	foreach $coord (@{$blocks{$poliInfo[0]}}) {

                	@cordinateSet = split("-", $coord);

	                if($type eq "SNP") {

	                        $pos = ($poliInfo[1] + $cordinateSet[1]);

	                        --$pos;

	                        print GFF3FILE "LG", $cordinateSet[0],"\t", $source,"\t", $type, "\t", $pos, "\t", $pos,"\t.\t.\t.\tName=", $poliInfo[0],"-", $poliInfo[1], ";Alleles=", $poliInfo[2], ";Alleles=", $poliInfo[3];
	                        print GFF3FILE "\n";

	                } elsif ($type eq "SPP") {

	                        $startPos = $poliInfo[1] + $cordinateSet[1];

	                        --$startPos;

	                        $endPos = $poliInfo[2] + $cordinateSet[1];

	                        --$endPos;

	                        print GFF3FILE "LG", $cordinateSet[0],"\t", $source,"\t", $type, "\t", $startPos, "\t", $endPos,"\t.\t.\t.\tName=", $poliInfo[0],"-", $poliInfo[1],"_", $poliInfo[2];
	                        print GFF3FILE "\n";

	                } else {

	                        print "Couldn't identify the type of polimorphism, please input SNP or SPP.\n";
	                        die "Wrong input type";

	                } #end else

		} #end foreach

	} #end if

} #end while

close (POLIMORPHISMFILE);
close (GFF3FILE);

end;
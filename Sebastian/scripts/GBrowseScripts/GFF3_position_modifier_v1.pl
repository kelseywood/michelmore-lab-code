#!/usr/bin/perl

######################################################################
######################################################################
#                                                                    #
#                                                                    #
#               Script design to rePosition tracks                   #
#                                                                    #
#                    Sebastian Reyes Chin-Wo                         #
#                   Seed Biotechnology Center                        #
#                University of California at Davis                   #
#                      sreyesch@ucdavis.edu                          #
#                                                                    #
#                                                                    #
######################################################################
######################################################################

if (@ARGV < 1) {

	print "Program usage:\n";
	print "Two input arguments needed, please enter them.\n";
	print "- First argument: Coordinates file with the positions and LG of all the contigs on the sequence.\n";
	print "- Next arguments: GFF3 files to modify.\n";
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

my $coordinates = shift(@ARGV);

open (COORDINATES, $coordinates);

while (my $line = <COORDINATES>) {

	chomp($line);

	my @coordinate = split("\t", $line);

	push (@{$blocks{$coordinate[0]}}, "$coordinate[1]\t$coordinate[2]\t$coordinate[3]\t$coordinate[4]");

} #end while

close (COORDINATES);

####################################
####################################
#                                  #
#          Open GFF3 File          #
#                                  #
####################################
####################################

foreach my $gff (@ARGV) {

	open (GFF3INPUT, $gff);

	my $base = $gff;

	$base =~ s{.*/}{};
	$base =~ s{\.[^.]+$}{};

	open (GFF3OUTPUT, ">$base.rePositioned.gff3");

	while(my $GFF3LINE = <GFF3INPUT>) {

		chomp($GFF3LINE);

		my @GFF3Info = split("\t", $GFF3LINE);

		if($blocks{$GFF3Info[0]}[0] ne "") {

		        foreach $coord (@{$blocks{$GFF3Info[0]}}) {

		        	my @cordinateSet = split("\t", $coord);

				my $startPos;
				my $endPos;

		                if($cordinateSet[1] eq "+") {

					$startPos = ($GFF3Info[3] + $cordinateSet[2]);

					$endPos = ($GFF3Info[4] + $cordinateSet[2]);

					--$startPos;
					--$endPos;

		                } else {

					$startPos = ($cordinateSet[3] - $GFF3Info[4]);

					$endPos = ($cordinateSet[3] - $GFF3Info[3]);

					++$startPos;
					++$endPos;

		                } #end else

		                if ($GFF3Info[6] eq "-" && $cordinateSet[1] eq "-") {

		                	$strand = "+";

		                } elsif ($GFF3Info[6] eq "+" && $cordinateSet[1] eq "+") {

		                	$strand = "+";

		                } elsif ($GFF3Info[6] eq "+" && $cordinateSet[1] eq "-") {

		                	$strand = "-";

		                } elsif ($GFF3Info[6] eq "-" && $cordinateSet[1] eq "+") {

		                	$strand = "-";

		                } elsif ($GFF3Info[6] eq "+" && $cordinateSet[1] eq ".") {

		                	$strand = "+";

		                } elsif ($GFF3Info[6] eq "-" && $cordinateSet[1] eq ".") {

		                	$strand = "-";

		                } else {

		                	$strand = ".";

		                } #end else

				print GFF3OUTPUT $cordinateSet[0],"\t", $GFF3Info[1], "\t", $GFF3Info[2], "\t", $startPos, "\t", $endPos,"\t", $GFF3Info[5], "\t", $strand, "\t", $GFF3Info[7], "\t", $GFF3Info[8], "Scaffold=", $GFF3Info[0],";";
				print GFF3OUTPUT "\n";

		        } #end foreach

		} #end if

	} #end while

	close (GFF3INPUT);
	close (GFF3OUTPUT);

} #end looping through gff's

exit;

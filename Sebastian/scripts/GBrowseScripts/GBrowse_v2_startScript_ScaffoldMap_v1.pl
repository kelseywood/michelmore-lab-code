######################################################################
######################################################################
#                                                                    #
#                                                                    #
#         Script design to generate "Whole Genome Sequence"          #
#          for display on Generic Genome Browser (GBrowse)           #
#         base on individual sequences and genetic distances         #
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

# The script is desgined to paste individual sequences (EST or contig) to generate
# an appropiate sequence to display on the GBrowse, it used genetic distance between
# them to calculculate the number of N's that separete each sequence. Also it generate
# the partner GFF file with the cordinates for each of the sequence and a complementary
# table with the start potition of each sequence for potitioning of other
# information linked to the sequences (SNP, genes, etc.)

if (@ARGV < 4) {

	print "Program usage:\n";
	print "Seven input arguments needed, please enter them.\n";
	print "- First argument: fasta file with the sequences to assemble.\n";
	print "- Second argument: Map file with the genetic distance between the sequences.\n";
	print "- Third argument: Assembly version\n";
	print "- Forth argument: Number of bases per cM\n";
	die "Input arguments missing";

} #end if


################################
#                              #
# Global and storage variables #
#                              #
################################

my %sequences = ();
my $numberSequences = 0;
my @wholeGenomicSequence = ();
my $curPosition = 1;
my $numberNs = 0;
my $featureNumber = 1;
my $prevMapPos = 0;
my $endPos = 0;
my $curLG = ();
my %LGpositions = ();
my $curBin = ();
my $prevEndPos = 0;

my $cM2bases = $ARGV[3];

############################
############################
#                          #
#   Open Input files and   #
#       load hashes        #
#                          #
############################
############################

#######################
#                     #
#   Open fasta file   #
#                     #
#######################

# Extract the sequences from the fasta file and
# loaded into the sequences hash

open (FASTA, $ARGV[0]);

while ($fastaLine0 = <FASTA>) {

	$fastaLine = reverse($fastaLine0);

	chop($fastaLine);

	$fastaLine1 = reverse($fastaLine);

	@header = split(" ", $fastaLine1);

	$curSequence = <FASTA>;

	chomp($curSequence);

	$sequences {$header[0]} = $curSequence;

	$numberSequences++;

} #end while

close (FASTA);

print "$numberSequences were loaded from $ARGV[0]\n";

####################################
####################################
#                                  #
#   Open Map file and construct    #
#   the "Whole Genomic Sequence"   #
#                                  #
####################################
####################################

$version = sprintf("%02d", $ARGV[2]);

open (GFF3FILE, ">$ARGV[1].markers.gff3");

open (MAP, $ARGV[1]);
open (SEQUENCE, ">WholeGenomicSequence.v$version.fa");
open (COORDINATESFILE, ">$ARGV[1]_coordinates.txt");

while($mapLine = <MAP>) {

	chop($mapLine);

	@lociInfo = split("\t", $mapLine);

	if($curLG ne $lociInfo[0]) {

		print SEQUENCE "\n";
		print SEQUENCE ">LG", $lociInfo[0], "\n";

		$LGpositions {$curLG} = $endPos;

		$featureNumber = 1;
		$prevMapPos = 0;
		$curPosition = 1;

	} #enf if

	$curScaffold = $lociInfo[1];

	$curSequence = $sequences {$curScaffold};

	my $lengthSequence = length ($curSequence);

	$mapDistance = $lociInfo[2] - $prevMapPos;

	if ($mapDistance eq 0) {

		$numberNs = 1000;

	} else {

		$numberNs0 = $cM2bases * $mapDistance;

		$numberNs = sprintf("%.0f", $numberNs0);

	} #end else

        $seqReps = $numberNs/8;

	my $Nsequence = "NANCNGNT" x $seqReps;

	print SEQUENCE $Nsequence;

        if($lociInfo[3] eq "+") {

		print SEQUENCE $curSequence;

        } else {

       		print SEQUENCE reverse($curSequence);

        } #end else

	$posPlusN = $curPosition + $numberNs;

	$curPosition = $posPlusN;

	$endPos = $curPosition + $lengthSequence;

	$endPos--;

	$featureID = sprintf("%07d", $featureNumber);

	@GFFID = ("LS_LG", $lociInfo[0], "_V", $version, "__", $featureID);

	print GFF3FILE "LG", $lociInfo[0],"\t.\tscaffold\t", $curPosition, "\t", $endPos,"\t.\t", $lociInfo[3], "\t.\t;ID=", @GFFID, ";Name=",  @GFFID, ";Alias=", $curScaffold, "\n";

	print COORDINATESFILE $curScaffold, "\t", $lociInfo[0], "\t", $curPosition, "\n";

	$endPos++;

	$curPosition = $endPos;
	$prevEndPos = $endPos;

	$featureNumber++;

	$prevMapPos = $lociInfo[2];

	$curLG = $lociInfo[0];

} #end while

$LGpositions {$curLG} = $endPos;
delete $LGpositions {""};

open(CHROMOSOMES, ">position_linkage_groups.gff3");

foreach $k (keys (%LGpositions)) {

	print CHROMOSOMES "LG", $k,"\t.\tdna\t1\t", $LGpositions{$k},"\t.\t.\t.\tName=LG", $k, ";ID=LG", $k, "\n";

} #end foreach

close(CHROMOSOMES);
close (MAP);
close (GFF3FILE);
close (COORDINATESFILE);
close (SEQUENCE);

end;





















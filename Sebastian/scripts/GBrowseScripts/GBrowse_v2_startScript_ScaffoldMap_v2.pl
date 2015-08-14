#!/usr/bin/perl
use strict;
use warnings;

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
	print "		With the following fields.\n";
	print "		LinkegeGroup	ScaffoldName Bin	Position	Strand\n";
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
my $curLG = "";
my %LGpositions = ();
my $curBin = "";
my $prevEndPos = 0;

my $cM2bases = ($ARGV[3]/8);


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

while (my $fastaLine = <FASTA>) {

	chomp($fastaLine);

	$fastaLine =~ s/>//;

	my @header = split(" ", $fastaLine);

	my $curSequence = <FASTA>;

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

my $version = $ARGV[2];


open (COORDINATESFILE, ">Lsat.1.v$version.lg.coordinates.txt");

open (GFF3FILE, ">Lsat.1.v$version.scaffolds.gff3");
open (BINFILE, ">Lsat.1.v$version.bins.gff3");

open (MAP, $ARGV[1]);
open (SEQUENCE, ">Lsat.1.v$version.lg.fasta");

while(my $mapLine = <MAP>) {

	chop($mapLine);

	my @lociInfo = split("\t", $mapLine);

	if($curLG ne $lociInfo[0]) {

		print SEQUENCE "\n";
		print SEQUENCE ">Lsat_1_v", $version, "_lg_", $lociInfo[0], "\n";

		$LGpositions {$curLG} = $endPos;

		$featureNumber = 1;
		$prevMapPos = 0;
		$curPosition = 1;

	} #enf if

	my $curScaffold = $lociInfo[1];

	my $curSequence = $sequences {$curScaffold};

	if ($sequences {$curScaffold} eq "" ) {

		die "Sequence for $curScaffold not found";

	} #end if

	my $lengthSequence = length ($curSequence);

	my $mapDistance = $lociInfo[3] - $prevMapPos;

	if ($mapDistance eq 0) {

		$numberNs = 125;

	} else {

		$numberNs = sprintf("%.0f", ($cM2bases * $mapDistance));

	} #end else

	my $Nsequence = "XAXCXGXT" x $numberNs;
#	my $Nsequence = "NANCNGNT" x $numberNs;

	print SEQUENCE $Nsequence;

        if($lociInfo[3] eq "+") {

		print SEQUENCE $curSequence;

        } else {

		my $revSequence = reverse($curSequence);

		$revSequence =~ y/ACGTacgt/TGCAtgca/;

       		print SEQUENCE $revSequence;

        } #end else

        my $unsequenceGap = length ($Nsequence);

	my $posPlusN = $curPosition + $unsequenceGap;

	$curPosition = $posPlusN;

	$endPos = $curPosition + $lengthSequence;

	$endPos--;

	my $featureID = sprintf("%07d", $featureNumber);

	print GFF3FILE "Lsat_1_v", $version, "_lg_", $lociInfo[0],"\t.\tscaffold\t", $curPosition, "\t", $endPos,"\t.\t", $lociInfo[4], "\t.\tID=", $curScaffold, ";Name=",  $curScaffold,";LG=", $lociInfo[0], ";GeneticPosition=", $lociInfo[3], ";Bin=", $lociInfo[2],";\n";

	print COORDINATESFILE $curScaffold, "\t", "Lsat_1_v", $version, "_lg_", $lociInfo[0], "\t", $lociInfo[3], "\t", $curPosition, "\t", $endPos, "\n";

        #Output data to the bin's GFF file
	if($curBin eq "") {

		print BINFILE "Lsat_1_v", $version, "_lg_", $lociInfo[0],"\t.\tbin\t", $curPosition, "\t";

	} elsif ($curBin ne $lociInfo[2]) {

		print BINFILE $prevEndPos,"\t.\t.\t.\tID=", $curBin, ";Name=",  $curBin,";Bin=", $curBin, ";\n";
		print BINFILE "Lsat_1_v", $version, "_lg_", $lociInfo[0],"\t.\tbin\t", $curPosition, "\t";

	} #end elseif


	$endPos++;

	$curPosition = $endPos;
	$prevEndPos = $endPos;

	$featureNumber++;

	$prevMapPos = $lociInfo[3];

	$curLG = $lociInfo[0];

	$curBin = $lociInfo[2];

} #end while

print BINFILE $prevEndPos,"\t.\t.\t.\tID=", $curBin, ";Name=",  $curBin,";bin=", $curBin, ";\n";

$LGpositions {$curLG} = $endPos;
delete $LGpositions {""};

open(CHROMOSOMES, ">Lsat.1.v$version.lg.gff3");

foreach my $k (keys (%LGpositions)) {

	print CHROMOSOMES "Lsat_1_v", $version, "_lg_", $k,"\t.\tdna\t1\t", $LGpositions{$k},"\t.\t.\t.\tName=Lsat_1_v", $version, "_lg_", $k, ";ID=Lsat_1_v", $version, "_lg_", $k, "\n";

} #end foreach

close (CHROMOSOMES);
close (MAP);
close (GFF3FILE);
close (COORDINATESFILE);
close (SEQUENCE);

exit;

















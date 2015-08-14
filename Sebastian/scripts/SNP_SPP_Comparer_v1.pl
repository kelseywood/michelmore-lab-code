#!/usr/bin/perl

if(@ARGV < 2) {

	print "Program Usage: Three arguments are needed, please input them.\n";
        print "First argument: file with SNP positions.\n";
        print "Second argument: file with SFPDev positions.\n";
        print "Third argument: output files prefix.\n";
        exit 0;

} #end if

my %SNPPos = ();
my %SNPData = ();
my %SNPcontigs = ();
my %SPPcontigs = ();

my $totalSPPPos = 0;
my $totalSNPPos = 0;

my $missSNP = 0;
my $missSPP = 0;
my $corSNP = 0;
my $corSPP = 0;

my %missSNPContigs = ();
my %missSPPContigs = ();

open(SNPFILE, $ARGV[0]);

while ($snpLine = <SNPFILE>) {

	chop($snpLine);

        @snpData = split("\t", $snpLine);

        $posId = "$snpData[0]-$snpData[1]";

#        print $posId, "\t";

        $SNPPos{$posId} = 1;

        $SNPData{$posId} = [@snpData];

        $SNPcontigs{$snpData[0]} = 1;

        ++$totalSNPPos;

} #end while

$numColumns = scalar (@snpData);

close (SNPFILE);

open(SPPFILE, $ARGV[1]);

open(CORRECTSPP, ">Correct_SFPCallings_$ARGV[2].txt");
open(FALSESPP, ">False_SFPCallings_$ARGV[2].txt");
open(MISSCONTIGS, ">MissContigs_$ARGV[2].txt");

my $detected = "";

while($sppLine = <SPPFILE>) {

	chop($sppLine);

        @sppData = split("\t", $sppLine);

        if ($SNPcontigs{$sppData[0]} ne "") {

	        $detected = "FALSE";

        	for ($pos = $sppData[1] - 8; $pos <= $sppData[2] + 8; ++$pos) {

	        	$posId = "$sppData[0]-$pos";

		        if ($SNPPos{$posId} ne "") {

	                        if ($detected eq "FALSE") {

        	           		print CORRECTSPP join ("\t", @sppData);
                                        print CORRECTSPP "\n";
                                        ++$corSPP;

                	        } #end if

	                   	$SNPPos{$posId} = 2;
        	                $detected = "TRUE";

                	} #end if

	        } #end for

        	if ($detected eq "FALSE") {

	       		print FALSESPP join ("\t", @sppData);
                        print FALSESPP "\n";
	                ++$missSPP;

		} #end if

        } else {

        	$missSPPContigs{$sppData[0]} = 1;

                print MISSCONTIGS $sppData[0], "\n";

       		print FALSESPP join ("\t", @sppData);
		print FALSESPP "\tNoSNPContig";
                print FALSESPP "\n";
                ++$missSPP;

        } #end else

        $SPPcontigs{$sppData[0]} = 1;

        ++$totalSPPPos;

} #end while

close(SPPFILE);
close(CORRECTSFP);
close(FALSESFP);

open(DETECTEDSNP, ">DectectedSNP_$ARGV[2].txt");
open(NONDETECTEDSNP, ">NonDectectedSNP_$ARGV[2].txt");

foreach $key (sort keys %SNPPos) {

#	print $key, "\t";
       	@posInfo = split("-", $key);

        if ($SPPcontigs{$posInfo[0]} ne "") {

		if ($SNPPos{$key} eq 1) {

			for ($i = 0; $i < $numColumns; ++$i) {
				print NONDETECTEDSNP $SNPData{$key}[$i], "\t";
			} #end foreach
                	print NONDETECTEDSNP "\n";
	                ++$missSNP;

        	} else {

			for ($i = 0; $i < $numColumns; ++$i) {
				print DETECTEDSNP $SNPData{$key}[$i], "\t";
			} #end foreach
        	        print DETECTEDSNP "\n";

                        ++$corSNP;

                } #end else

        } else {

                $missSNPContigs{$posInfo[0]} = 1;

        } #end else

} #end while

close (DETECTEDSNP);
close (NONDETECTEDSNP);

$numMissSNPContigs = scalar (keys %missSNPContigs);
$numMissSPPContigs = scalar (keys %missSPPContigs);
$numSNPContigs = scalar (keys %SNPcontigs);
$numSPPContigs = scalar (keys %SPPcontigs);


open(STATS, ">Stats_$ARGV[2].txt");

print STATS "Comparison ID: $ARGV[2]\n\n";

print STATS "Total number of SNP positions on ",  $ARGV[0], ": ", $totalSNPPos, "\n";
print STATS "Total number of contigs with SNPs: ", $numSNPContigs, "\n";
print STATS "Number of missed contigs from the SNPs: ", $numMissSNPContigs, "\n";
print STATS "Number of dectected SNPs: ", $corSNP, "\n";
print STATS "Number of not detected SNPs: ", $missSNP, "\n\n";

print STATS "Total number of SPP callings on ", $ARGV[1], ": ", $totalSPPPos, "\n";
print STATS "Total number of contigs with SPPs: ", $numSPPContigs, "\n";
print STATS "Number of missed contigs from the SPPs: ", $numMissSPPContigs, "\n";
print STATS "Number of false SPP callings: ", $missSPP, "\n";
print STATS "Number of correct SPP callings: ", $corSPP, "\n";

close(STATS);

exit;
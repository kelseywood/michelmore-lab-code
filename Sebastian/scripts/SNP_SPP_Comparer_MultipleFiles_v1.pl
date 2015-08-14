#!/usr/bin/perl

if(@ARGV < 2) {

	print "Program Usage: Three arguments are needed, please input them.\n";
        print "First argument: file with SNP positions.\n";
        print "Second argument: Directory with the SPP files.\n";
        print "Third argument: output files prefix.\n";
        exit 0;

} #end if

my %SNPPos = ();
my %SNPData = ();
my %SNPcontigs = ();
my $totalSNPPos = 0;

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

print "Finish loading SNP data\n\n";

$numColumns = scalar (@snpData);

close (SNPFILE);

my $Dir = $ARGV[1];
my @SPPFiles = ();

# open and read DIR
$numFiles = 0;
opendir(DIR, $Dir);

while ( defined( $file = readdir(DIR) ) ) {

	# skip '.' and '..'
	if ( $file ne "." && $file ne ".." ) {

       		push(@SPPFiles, $file);

        }    # end if

      	++$numFiles;

}    # end while

closedir(DIR);

print "$numFiles will be analyze against $ARGV[0]\n\n";

foreach $SPPfile (@SPPFiles) {

	print "Working on $SPPfile\n\n";

	my %SPPcontigs = ();
	my $totalSPPPos = 0;
	my $missSNP = 0;
	my $missSPP = 0;
	my $corSNP = 0;
	my $corSPP = 0;
	my %missSNPContigs = ();
	my %missSPPContigs = ();

	open(SPPFILE, $SPPfile);

	open(CORRECTSPP, ">Correct_SFPCallings_$ARGV[2]_$SPPfile.txt");
	open(FALSESPP, ">False_SFPCallings_$ARGV[2]_$SPPfile.txt");
	open(MISSCONTIGS, ">MissContigs_$ARGV[2]_$SPPfile.txt");

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
			print FALSESPP "NoSNPContig\t";
                	print FALSESPP "\n";
	                ++$missSPP;

        	} #end else

	        $SPPcontigs{$sppData[0]} = 1;

        	++$totalSPPPos;

	} #end while

	close(SPPFILE);
	close(CORRECTSFP);
	close(FALSESFP);

	$numMissSPPContigs = scalar (keys %missSPPContigs);
	$numSPPContigs = scalar (keys %SPPcontigs);

	open(STATS, ">Stats_$ARGV[2]_$SPPfile.txt");

	print STATS "Comparison ID: $ARGV[2]_$SPPfile\n\n";

	print STATS "Total number of SPP callings on ", $SPPfile, ": ", $totalSPPPos, "\n";
	print STATS "Total number of contigs with SPPs: ", $numSPPContigs, "\n";
	print STATS "Number of missed contigs from the SPPs: ", $numMissSPPContigs, "\n";
	print STATS "Number of false SPP callings: ", $missSPP, "\n";
	print STATS "Number of correct SPP callings: ", $corSPP, "\n";

	close(STATS);

} #end while

exit;
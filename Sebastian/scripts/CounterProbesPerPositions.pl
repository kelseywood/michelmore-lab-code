
open(INPUT1, $ARGV[0]) || die "Unable to open the file $ARGV[0]"; # First input, file with the contigs length

open(RESULTS, ">$ARGV[2]_ext$ARGV[3].txt") || die "Unable to open the results file";

my %pips = ();

while ($line = <INPUT1>) {

	chop($line);
	chop($line);

	@PIPInfo = split ("\t", $line);

        $posID = "$PIPInfo[0].$PIPInfo[1]";

#        print $posID, "\n";

	$pips{$posID} = "TRUE";

#        print $pips{$posID};

}# end while

close(INPUT1);

print "\n\n1", $pips{$posID};
print "\n\n2", $pips{QGJ9P24.yg.ab1.247}, "\n";
print "\n3", $posID, "\n";

my $coverage = 0;

my $extension = 12; # third input with the extension of the probe coverage from the PIP

print RESULTS "Contig\tPos\tDepth\n";

open (INPUT2, $ARGV[1]) || die "Unable to open the file $ARGV[1]"; # Second input, file with the Probe PIP's

while ($line2 = <INPUT2>) {

#	print "\n\n4", $pips{$posID};

	chop($line2);

	my @contigInfo = split ("\t", $line2);

 	my $curLength = $contigInfo[1];

	for ( my $pos = 1; $pos <= $curLength; ++$pos) {

		for ($i = $pos - $extension; $i <= $pos + $extension; ++$i) {

                	$posId = "$contigInfo[0].$i";

#                        print "_", $posId, "_", $posID, "_", "\n";
#                        print $pips{CGP_A0000000001_12}, "\n";

                        ++$coverage if defined $pips{$posId};

		}# end foreach

	print RESULTS $contigInfo[0], "\t", $pos, "\t", $coverage, "\n";

	$coverage = 0;

	}# end for

}# end while

close(INPUT2);
close(RESULTS);
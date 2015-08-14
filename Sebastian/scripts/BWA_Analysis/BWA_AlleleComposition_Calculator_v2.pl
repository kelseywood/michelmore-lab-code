if (@ARGV < 1) {

	print "Program usage: pileup file is require, please input it.\n";
	exit 0;

}

open(INPUT1, $ARGV[0]) || die "Cannot open input file $ARGV[0]"; # File with maq Pileup Results

$fileName = $ARGV[0];

open(RESULTS, ">$ARGV[0]_Allele_Comp.txt");

print RESULTS "Contig\tPos\tRef\tDepth\tA\tC\tG\tT\n";

while($line = <INPUT1>) {

	chop($line);

	@pos = split("\t", $line);

		###################
		#                 #
		#   Calculate     #
		#  allele comp    #
		#                 #
		###################

	 	$readprofile0 = $pos[4];

        	$readprofile = uc($readprofile0);

		@readsBases = split("", $readprofile);

		$countA = 0;
		$countC = 0;
		$countG = 0;
		$countT = 0;

        	my $indel = 0;

		foreach (@readsBases) {

        		if ($indel eq 0) {

				if ($_ eq "." || $_ eq ",") {
					$_ = $pos[2];
				}

				if($_ eq "A") {
					++$countA;
				} elsif ($_ eq "C") {
					++$countC;
				} elsif ($_ eq "G") {
					++$countG;
				} elsif ($_ eq "T") {
					++$countT;
	 			} elsif (/^[+-]?\d+$/) {
					$indel = $_;
				} #end elsif

			} else {

        	        	--$indel;

                	} #end elsif

		} #end foreach

        	if ($pos[2] ne "*") {

			print RESULTS $pos[0], "\t", $pos[1], "\t", $pos[2], "\t", $pos[3], "\t", $countA,  "\t", $countC,  "\t", $countG,  "\t", $countT, "\n";

	        } #end else

} #end while

close(RESULTS);
close(INPUT1);

end;
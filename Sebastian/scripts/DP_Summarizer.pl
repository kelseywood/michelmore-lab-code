#Open DP results file
open(INPUT1, $ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # the file containing the DP results
open(RESULTS,">$ARGV[0]_summarize.txt") || die "Cannot open the Results file"; # output file to write spp's information

#Create Storage Variables
my $prevContig = "";
my $prevEnd = 0;

print RESULST "ContigID\tStarts\tEnd\n";

#Start loop read file
while($line = <INPUT1>) {

	@SPP = split("\t",$line);

	#Check for working Contig

	if(@SPP[0] eq $prevContig) {

		#Same Contig

		$starts = @SPP[1];
		$starts_8 = $starts - 8;

		#Check End of past row, with end of current row

			if($starts_8<$prevEnd && $prevEnd<$starts) {
				#If they are in range of +50
				#Change End of last for End of current

				$prevEnd = @SPP[2];

				$prevContig = @SPP[0];

			} else {
				#If not print the end pos

				print RESULTS $prevEnd;

				#Start new line and storage variables
				print RESULTS @SPP[0],"\t",@SPP[1],"\t";

				$prevEnd = @SPP[2];

				$prevContig = @SPP[0];

			}

	} else {
		#Different Contig

		if($prevContig ne "") {

			#Print previous information
			print RESULTS $prevEnd;

		}

		#Start new line and storage variables
		print RESULTS @SPP[0],"\t",@SPP[1],"\t";

		$prevEnd = @SPP[2];

		$prevContig = @SPP[0];

	}

}

print RESULTS $prevEnd;

close(INPUT1);
close(RESULTS);
#End
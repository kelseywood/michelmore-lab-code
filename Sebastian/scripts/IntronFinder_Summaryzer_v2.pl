open(INPUT1,$ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # the file containing the results from intron finder

open(RESULTSINTRONS,">IntronPositions_PerContig.txt")|| die "Cannot open the Results file"; # output file to write intron positions per contig

my $contig = "";

while(<INPUT1>) {

	chomp(my $line = $_);
	$forPrint = $line;
	my @line = split('\s', $line);

	#Check if the current line is the line with the contig name and print
	if(@line[0] eq "Results") {

		#Eliminate additional characters from the line
		$contig = substr($forPrint, 12, -1);

	#Check if the line has the intron positions and print
	} elsif (@line[0] eq "positions:") {

		#Eliminate additional characters from the line
		foreach $value (@line) {

                	if ($value =~ /^[+-]?\d+$/) {

				print RESULTSINTRONS $contig, "\t";
				print RESULTSINTRONS $value, "\n";

                        } #end if
          	} #end foreach

	}# end elseif

}#end while

close (INPUT1);
close (RESULTSINTRONS);
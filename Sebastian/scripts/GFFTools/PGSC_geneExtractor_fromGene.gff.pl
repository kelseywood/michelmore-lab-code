
open(FASTA, $ARGV[0]);

print "ChrID\tGeneID\tStart\tEnd\tStrand\n";

while ($line = <FASTA>) {

	chomp($line);

	@fields = split("\t", $line);

        my @attributes = split(";", $fields[8]);

        my $ID = substr($attributes[0], 3, 100);

	if($fields[2] eq "mRNA") {

#		my $midPosition = ($fields[4]-$fields[3])/2;

		print $fields[0], "\t", $ID, "\t", $fields[3], "\t", $fields[4], "\t", $fields[6], "\n";

	} #end if

} #end while

close(FASTA);

exit;
		 

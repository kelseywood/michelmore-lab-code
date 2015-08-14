# Expected fasta header format
# ID Description [Specie Genus]


open(FASTA, $ARGV[0]);

my %fastaDescription;

while($fastaLine = <FASTA>) {

	chomp($fastaLine);

	$firstCharacter = substr($fastaLine,0,1);

#	print "$firstCharacter\n";

	if($firstCharacter eq ">") {

		$header = substr($fastaLine,1,200);

		chop($header);

		@data = split(" ", $header);

		$ID = $data[0];

#		print $ID, "\n";

		$genus = pop(@data);

		$specie = pop(@data);

		$fastaDescription{$ID} = ($genus, substr($specie,2,20), join(" ", @data));

	} #end if

} #end while

close(FASTA);

open(GFF, $ARGV[1]);

while($GFFLine = <GFF>) {

		chomp($GFFLine);

		@fields = split("\t", $GFFLine);

		if ($fields[1] eq "blastx" && $fields[2] eq "protein_match") {

			@atrributes = split(";", $fields[8]);

			$id = substr($atrributes[1], 5, 100);

#			print $id, "\n";

			print $GFFLine,"Specie=", $fastaDescription[$id][0], " ", $fastaDescription[$id][1], ";Description=", $fastaDescription[$id][2], ";\n";

		} else {

			print $GFFLine, "\n";

		} #end else

} #end while

close(GFF);

exit;
		


# Expected fasta header format
# ID Description [Specie Genus]


$fasta = $ARGV[0];
$fasta = $opt_f if $opt_f;

# TRY TO MAKE A REGULAR EXPRESSION WITH THE ENTIRE SPECIE NAME

system("cp $fasta temp.fasta");
system("perl -p -i -e 's/\\| /\t/g' temp.fasta");

#Open the fasta file and load header data into hash

open(FASTA, "temp.fasta");

my %tairDescription;

while($fastaLine = <FASTA>) {

	chomp($fastaLine);

	$firstCharacter = substr($fastaLine,0,1);

	if($firstCharacter eq ">") {

		$header = substr($fastaLine,1,200);

		@data = split("\t", $header);

		$data[0] =~ s/ //g;
		$data[2] =~ s/;/ \| /g;
		$data[3] =~ s/;/ \| /g;

#		print $data[0], "\n";

		$tairDescription{$data[0]} = [$data[2], $data[3]];

	} #end if

} #end while

close(FASTA);
system("chmod 666 temp.fasta");
system("rm temp.fasta");


open(GFF, $ARGV[1]);

while($GFFLine = <GFF>) {

		chomp($GFFLine);

		@fields = split("\t", $GFFLine);

		if ($fields[2] eq "mRNA") {

			@atrributes = split(";", $fields[8]);

			@id = split("-", substr($atrributes[1],3,100));

#			print "$id[0]", "\n";

			if(defined($tairDescription{$id[0]}[0])) {

				print $GFFLine,"Description=", $tairDescription{$id[0]}[0],";Location=", $tairDescription{$id[0]}[1], ";\n";

			} else {

				print $GFFLine, "\n";

			} #end else

		} else {

			print $GFFLine, "\n";

		} #end else

} #end while

close(GFF);

exit;
		


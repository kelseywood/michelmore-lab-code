# Expected fasta header format
# ID Description [Specie Genus]


open(TAIR, $ARGV[0]) or die "can't open file $ARGV[0]\n";

my %BGIDescription;

while($line = <TAIR>) {

	chomp($line);

	@data = split("\t", $line);

	$data[15] =~ s/;/ \| /g;

	$BGIDescription{$data[0]} = [$data[4],$data[8],$data[13],$data[15]];

#	print $tairDescription{$data[0]}[0], "\n";

} #end while

close(TAIR);

open(GFF, $ARGV[1]);

while($GFFLine = <GFF>) {

		chomp($GFFLine);

		@fields = split("\t", $GFFLine);

		if ($fields[2] eq "mRNA") {

			@atrributes = split(";", $fields[8]);

			$id = substr($atrributes[0],3,100);

#			print "$id[0]", "\n";

			if(defined($BGIDescription{$id}[0])) {

				print $GFFLine,"$ARGV[2]-Organism=", $BGIDescription{$id}[0],";$ARGV[2]-Identity=", $BGIDescription{$id}[1],";$ARGV[2]-Expected value=", $BGIDescription{$id}[2],";$ARGV[2]-Description=", $BGIDescription{$id}[3], ";\n";

			} else {

				print $GFFLine, "\n";

			} #end else

		} else {

			print $GFFLine, "\n";

		} #end else

} #end while

close(GFF);

exit;
		


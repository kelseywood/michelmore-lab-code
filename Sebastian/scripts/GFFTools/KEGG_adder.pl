#!/usr/bin/perl

if (@ARGV < 1) {

	print "Two are arguments are needed, please input them.\n";
	print "GFF File were to add the KEGG annotations.\n";
        print "KEGG File.\n\n";
	print "perl Fasta_renamer.pl <KEGGFile> <GFFFile>\n";
        exit 0;

} #end if


my %KEGG_data;

system("cp $ARGV[0] kegg.temp");
system("perl -p -i -e 's/;/\\|/g' kegg.temp");

open(KEGGFILE, "kegg.temp");
while ($line = <KEGGFILE>) {

	chomp($line);

	@KEGGs = split("\t", $line);

	$KEGG_data{$KEGGs[0]} = [$KEGGs[4], $KEGGs[9], $KEGGs[13], $KEGGs[15]];

} #end while

close(GOFILE);
system("chmod 666 kegg.temp");
system("rm kegg.temp");

open(GFF, $ARGV[1]);

while($GFFLine = <GFF>) {

		chomp($GFFLine);

		@fields = split("\t", $GFFLine);

		if ($fields[2] eq "mRNA") {

			@atrributes = split(";", $fields[8]);

			my $id="";

			foreach $attribute (@atrributes) {
	
				@value = split("=", $attribute);

			        if($value[0] eq "ID") {
					$id = $value[1];
				} #end if

			} #end foreach

			if(defined($KEGG_data{$id}[0])) {

				print $GFFLine,"KEGG_organism=", $KEGG_data{$id}[0], ";","KEGG_identity=", $KEGG_data{$id}[1], ";", "KEGG_pValue=", $KEGG_data{$id}[2], ";", "KEGG_description=", $KEGG_data{$id}[3], ";\n";

			} else {

				print $GFFLine, "\n";

			} #end else

		} else {

			print $GFFLine, "\n";

		} #end else

} #end while

close(GFF);

exit;

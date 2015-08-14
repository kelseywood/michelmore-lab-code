#!/usr/bin/perl

if (@ARGV < 1) {

	print "Two are arguments are needed, please input them.\n";
	print "GFF File were to add the GO annotations.\n";
        print "Go filee.\n\n";
	print "perl Fasta_renamer.pl <GFFFile> <GOFile>\n";
        exit 0;

} #end if

my %GO_genes;

open(GOFILE, $ARGV[0]);
while ($line = <GOFILE>) {

	chomp($line);

	@GOs = split("\t", $line);

#	print $GOs[0], "\n";

	for($i = 2; $i < (1 + $GOs[1]); ++$i) {

		$GOs[$i] =~ s/;/:/g;

		push(@{$GO_genes{$GOs[0]}}, $GOs[$i]);
	
	} #end for

} #end while

close(GOFILE);

open(GFF, $ARGV[1]);

while($GFFLine = <GFF>) {

		chomp($GFFLine);

		@fields = split("\t", $GFFLine);

		if ($fields[2] eq "mRNA") {

			@atrributes = split(";", $fields[8]);

			$id = substr($atrributes[0],3,100);

#			print "$id[0]", "\n";

			if(defined($GO_genes{$id}[0])) {

				print $GFFLine,"GO_annotations=", join(" \| ",@{$GO_genes{$id}}), ";\n";

			} else {

				print $GFFLine, "\n";

			} #end else

		} else {

			print $GFFLine, "\n";

		} #end else

} #end while

close(GFF);

exit;

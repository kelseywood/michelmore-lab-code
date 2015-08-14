#!/usr/bin/perl
use warnings;
use Getopt::Std;
use vars qw($opt_d $opt_f);
getopts('d:f');

#################################
#                               #
#     Get Header Line Data      #
#                               #
#################################

#Reformat the headers line to made them readable by split
#removing special characters [ and |


$fasta = $ARGV[0];
$fasta = $opt_f if $opt_f;

# TRY TO MAKE A REGULAR EXPRESSION WITH THE ENTIRE SPECIE NAME

system("cp $fasta temp.fasta");
system("perl -p -i -e 's/\\| /\|\t/g' temp.fasta");
system("perl -p -i -e 's/\\[/\t/g' temp.fasta");
system("perl -p -i -e 's/\\]//g' temp.fasta");

#Open the fasta file and load header data into hash

open(FASTA, "temp.fasta");

my %fastaDescription;

while($fastaLine = <FASTA>) {

	chomp($fastaLine);

	$firstCharacter = substr($fastaLine,0,1);

	if($firstCharacter eq ">") {

		$header = substr($fastaLine,1,200);

		@data = split("\t", $header);

		$data[1] =~ s/;/ \| /g;

		$fastaDescription{$data[0]} = [$data[2], $data[1]];

	} #end if

} #end while

close(FASTA);
system("chmod 666 temp.fasta");
system("rm temp.fasta");

#Open current GFF file and temporal output GFF				
open(GFF, $ARGV[1]) or print "Can't open file $ARGV[1]";

while($GFFLine = <GFF>) {

	chomp($GFFLine);

	@fields = split("\t", $GFFLine);

	if (($fields[1] eq "blastx" && $fields[2] eq "protein_match") || ($fields[1] eq "protein2genome" && $fields[2] eq "protein_match")) {

		@atrributes = split(";", $fields[8]);

		$id = substr($atrributes[1], 5, 100);

		if($fastaDescription{$id}[0] ne "") {

			print $GFFLine,"Specie=", $fastaDescription{$id}[0], ";Description=", $fastaDescription{$id}[1], "\n";

		} else {

			print STDERR "No specie found for $id, printing original GFF line\n";
			print $GFFLine, "\n";

		} #end else

	} else {

		print $GFFLine, "\n";

	} #end else

} #end while

close(GFF);

exit;

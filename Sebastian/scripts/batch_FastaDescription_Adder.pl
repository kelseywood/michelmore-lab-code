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


$fasta = "proteins.fa";
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

		$fastaDescription{$data[0]} = [$data[2], $data[1]];

	} #end if

} #end while

close(FASTA);

#Remove temp.fasta (was just junk)
system("chmod 666 temp.fasta");
system("rm temp.fasta");


##############################
#                            #
#      Open directories      #
#    and add description     #
#                            #
##############################

my @chr = qw(1 2 3 4 5 6 7 8 9 A); #
my $Dir = "GFF/$chr/$i/";
my $Dir = $opt_d if $opt_d

foreach my $chr (@chr) {
	print STDERR "processing chr $chr\n";

	for (my $i = 0; $i < 10; $i++) {

		print STDERR "\tprocessing chunk $i\n";

		system("mkdir GFF/$chr/$i") unless -d "GFF/$chr/$i";

		# open and read DIR
		my @myFiles = "";

		while ( defined( my $file = readdir(DIR) ) ) {

			# skip '.' and '..'
			if ( $file ne "." && $file ne ".." ) {

		       		push( @myFiles, "$Dir$file");

	        	}    # end if

		}    # end while

		closedir(DIR);

		#Do printing and renaming loop per file

		foreach my $file (@myFiles) {

			#if the file has a name
			if ($file ne "") {

				#Open current GFF file and temporal output GFF				
				open(GFF, $file) or print "Can't open file $file";
				$tempFile = "$file.temp";
				open(TEMPGFF, ">$tempFile");

				while($GFFLine = <GFF>) {

						chomp($GFFLine);

						@fields = split("\t", $GFFLine);

						if (substr($GFFLine,0,1) eq "L" && $fields[1] eq "blastx" && $fields[2] eq "protein_match") {

							@atrributes = split(";", $fields[8]);

							$id = substr($atrributes[1], 5, 100);

							if($fastaDescription{$id}[0] ne "") {

								print TEMPGFF $GFFLine,"Specie=", $fastaDescription{$id}[0], ";Description=", $fastaDescription{$id}[1], "\n";

							} else {

								print STDERR "No specie found for $id, printing original GFF line\n";
								print TEMPGFF $GFFLine, "\n";

						} else {

							print TEMPGFF $GFFLine, "\n";

						} #end else

				} #end while

				close(GFF);

				system("mv $tempFile $file");

			} #end if
		} #end foreach
	} #end for
} #end foreach

exit;

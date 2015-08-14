#!/usr/bin/perl
use warnings;

################################
#                              #
# Global and storage variables #
#                              #
################################

my %sequences = ();
my %GFF;

my @file = (
'All_All_postMaker_exonerate_plusLink.gff3',
'Lsat.1.v3.AffyContigs.gmap.onePath.gff3',
'Lsat.1.v3.geneModels.gff3',
'Lsat.1.v3.genePredict.AbInitio.gff3',
'Lsat.1.v3.genePredict.geneWise_plusLink.gff3',
'Lsat.1.v3.MapPositions.gff3',
'Lsat.1.v3.repeatElements.gff3',
'Lsat.1.v3.ncRNA.gff3',
#'Lsat.1.v3.scaffolds.gff3',
'unMaker_scaffolds.gff'
);


############################
############################
#                          #
#   Open Input files and   #
#       load hashes        #
#                          #
############################
############################

#######################
#                     #
#   Open fasta file   #
#                     #
#######################

# Extract the sequences from the fasta file and
# loaded into the sequences hash

open (FASTA, "/Data/Genome_Data_Final/Lsat.1.v3.ScaffoldsSequence.bigger_1kb.fasta");

while ($fastaLine0 = <FASTA>) {

	$fastaLine = reverse($fastaLine0);

	chop($fastaLine);

	$fastaLine1 = reverse($fastaLine);

	@header = split(" ", $fastaLine1);

	$curSequence = <FASTA>;

	chomp($curSequence);

	$sequences {$header[0]} = $curSequence;

} #end while

close (FASTA);
print STDERR "Fasta file load it\n";


#######################
#                     #
#    Open GFF files   #
#                     #
#######################

# Extract the features from the GFF files and
# loaded it into the hash

foreach my $file (@file) {
	open(GFF, $file) or die "Can't open $file.\n";
	while (my $line = <GFF>) {
		next if $line =~ /#/;
		@fields = split("\t", $line);
		$line =~ s/Transposon/match/;
		$line =~ s/TEprotein/match/;
		$line =~ s/TandemRepeat/match/;

		$line =~ s/UTR_3/exon/;
		$line =~ s/UTR_5/exon/;

		$line =~ s/ESTs\tmRNA/EST\tmatch/;
		$line =~ s/PredictedExon/match_part/;

		$line =~ s/AUGUSTUS\tmRNA/AUGUSTUS\tmatch/;
		$line =~ s/AUGUSTUS\tCDS/AUGUSTUS\tmatch_part/;

		$line =~ s/GlimmerHMM\tmRNA/GlimmerHMM\tmatch/;
		$line =~ s/GlimmerHMM\tCDS/GlimmerHMM\tmatch_part/;

		$line =~ s/GeneWise_At\tmRNA/GeneWise_At\tmatch/;
		$line =~ s/GeneWise_At\tCDS/GeneWise_At\tmatch_part/;
		$line =~ s/GeneWise_Gm\tmRNA/GeneWise_Gm\tmatch/;
		$line =~ s/GeneWise_Gm\tCDS/GeneWise_Gm\tmatch_part/;
		$line =~ s/GeneWise_St\tmRNA/GeneWise_St\tmatch/;
		$line =~ s/GeneWise_St\tCDS/GeneWise_St\tmatch_part/;
		$line =~ s/GeneWise_Tc\tmRNA/GeneWise_Tc\tmatch/;
		$line =~ s/GeneWise_Tc\tCDS/GeneWise_Tc\tmatch_part/;
		$line =~ s/GeneWise_Vv\tmRNA/GeneWise_Vv\tmatch/;
		$line =~ s/GeneWise_Vv\tCDS/GeneWise_Vv\tmatch_part/;

		$line =~ s/blastn\trRNA/rRNA\tmatch/;
		$line =~ s/cmsearch\tmiRNA/miRNA\tmatch/;
		$line =~ s/cmsearch\tsnRNA/snRNA\tmatch/;
		$line =~ s/tRNAscan-SE\ttRNA/tRNA\tmatch/;

		$line =~ s/Map_Locations/match/;
		push @{$GFF{$fields[0]}}, $line;

		if ($fields[1] eq "GLEAN" && $fields[2] eq "CDS") {

			$line =~ s/CDS/exon/;
			push @{$GFF{$fields[0]}}, $line;

		} #end if

	}
	close GFF;

	print "$file load it into memory.\n";
	
}
print STDERR "GFF files load it\n";

$file = 'Lsat.1.v3.checkAnnotations.gff3';

my %annotatedScafolds;

open(GFF, $file) or die "Can't open $file.\n";

while (my $line = <GFF>) {

	next if $line =~ /#/;
	@fields = split("\t", $line);
	push @{$GFF{$fields[0]}}, $line;

	$annotatedScafolds{$fields[0]} = 1;

}
close GFF;

print "$file load it into memory.\n";


#######################
#                     #
#   Write GFF files   #
#                     #
#######################

# write GFF
my $i = 0;
my $n = keys %GFF;
my $m = int $n / 100;
print STDERR "writing $n files\n";

system("mkdir Database_hasApollo/") unless -d "Database_hasApollo/";
system("mkdir Database_hasApollo/0/") unless -d "Database_hasApollo/0/";
system("mkdir Database_hasApollo/0/Annotated/") unless -d "Database_hasApollo/0/Annotated/";
system("mkdir Database_hasApollo/0/UnAnnotated/") unless -d "Database_hasApollo/0/UnAnnotated/";
system("mkdir Database_hasApollo/1/") unless -d "Database_hasApollo/1/";
system("mkdir Database_hasApollo/1/Annotated/") unless -d "Database_hasApollo/1/Annotated/";
system("mkdir Database_hasApollo/1/UnAnnotated/") unless -d "Database_hasApollo/1/UnAnnotated/";
system("mkdir Database_hasApollo/2/") unless -d "Database_hasApollo/2/";
system("mkdir Database_hasApollo/2/Annotated/") unless -d "Database_hasApollo/2/Annotated/";
system("mkdir Database_hasApollo/2/UnAnnotated/") unless -d "Database_hasApollo/2/UnAnnotated/";
system("mkdir Database_hasApollo/3/") unless -d "Database_hasApollo/3/";
system("mkdir Database_hasApollo/3/Annotated/") unless -d "Database_hasApollo/3/Annotated/";
system("mkdir Database_hasApollo/3/UnAnnotated/") unless -d "Database_hasApollo/3/UnAnnotated/";
system("mkdir Database_hasApollo/4/") unless -d "Database_hasApollo/4/";
system("mkdir Database_hasApollo/4/Annotated/") unless -d "Database_hasApollo/4/Annotated/";
system("mkdir Database_hasApollo/4/UnAnnotated/") unless -d "Database_hasApollo/4/UnAnnotated/";
system("mkdir Database_hasApollo/5/") unless -d "Database_hasApollo/5/";
system("mkdir Database_hasApollo/5/Annotated/") unless -d "Database_hasApollo/5/Annotated/";
system("mkdir Database_hasApollo/5/UnAnnotated/") unless -d "Database_hasApollo/5/UnAnnotated/";
system("mkdir Database_hasApollo/6/") unless -d "Database_hasApollo/6/";
system("mkdir Database_hasApollo/6/Annotated/") unless -d "Database_hasApollo/6/Annotated/";
system("mkdir Database_hasApollo/6/UnAnnotated/") unless -d "Database_hasApollo/6/UnAnnotated/";
system("mkdir Database_hasApollo/7/") unless -d "Database_hasApollo/7/";
system("mkdir Database_hasApollo/7/Annotated/") unless -d "Database_hasApollo/7/Annotated/";
system("mkdir Database_hasApollo/7/UnAnnotated/") unless -d "Database_hasApollo/7/UnAnnotated/";
system("mkdir Database_hasApollo/8/") unless -d "Database_hasApollo/8/";
system("mkdir Database_hasApollo/8/Annotated/") unless -d "Database_hasApollo/8/Annotated/";
system("mkdir Database_hasApollo/8/UnAnnotated/") unless -d "Database_hasApollo/8/UnAnnotated/";
system("mkdir Database_hasApollo/9/") unless -d "Database_hasApollo/9/";
system("mkdir Database_hasApollo/9/Annotated/") unless -d "Database_hasApollo/9/Annotated/";
system("mkdir Database_hasApollo/9/UnAnnotated/") unless -d "Database_hasApollo/9/UnAnnotated/";
system("mkdir Database_hasApollo/A/") unless -d "Database_hasApollo/A/";
system("mkdir Database_hasApollo/A/Annotated/") unless -d "Database_hasApollo/A/Annotated/";
system("mkdir Database_hasApollo/A/UnAnnotated/") unless -d "Database_hasApollo/A/UnAnnotated/";

foreach my $scaf (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;

	if (exists($sequences{$scaf})) {

		if (exists($annotatedScafolds{$scaf})) {

			@scafName = split("_", $scaf);
			open(OUT, ">Database_hasApollo/$scafName[4]/Annotated/$scaf.gff") or die "Couldn't open file for Database_hasApollo/$scafName[4]/Annotated/$scaf.gff.\n";
			foreach my $line (@{$GFF{$scaf}}) {
				print OUT $line;
			}
			print OUT ">", $scaf, "\n";
			print OUT $sequences{$scaf}, "\n";
			close (OUT);

		} else {

			@scafName = split("_", $scaf);
			open(OUT, ">Database_hasApollo/$scafName[4]/UnAnnotated/$scaf.gff") or die "Couldn't open file for Database_hasApollo/$scafName[4]/UnAnnotated/$scaf.gff.\n";
			foreach my $line (@{$GFF{$scaf}}) {
				print OUT $line;
			}
			print OUT ">", $scaf, "\n";
			print OUT $sequences{$scaf}, "\n";
			close (OUT);

		} #end else

	} else {

		print "\nSequence not found for $scaf.\n";

	} #end else

#	system("tar -czf Database/$scafName[4]/$scaf.gff.tar.gz Database/$scafName[4]/$scaf.gff");
}
print STDERR "done\n";

































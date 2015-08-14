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
#'All_All_postMaker_exonerate_plusLink.gff3',
#'Lsat.1.v3.AffyContigs.gmap.onePath.gff3',
#'Lsat.1.v3-4.geneModels.gff',
#'Lsat.1.v3.genePredict.AbInitio.gff3',
#'Lsat.1.v3.genePredict.geneWise_plusLink.gff3',
#'Lsat.1.v3.MapPositions.gff3',
#'Lsat.1.v3.repeatElements.gff3',
#'Lsat.1.v3.ncRNA.gff3',
#'Lsat.1.v3.scaffolds.gff3',
#'unMaker_scaffolds.gff',
'All_All_postMaker_exonerate_plusLink.gff3.reLocated.gff.masked.gff.withoutOrphans.gff',
'Lsat.1.v3-4.geneModels.gff.withoutRedundant.gff.withoutRedundant.gff.reLocated.gff.masked.gff.clean.gff.withoutOrphans.gff',
'Lsat.1.v3.AffyContigs.gmap.onePath.gff3.reLocated.gff.masked.gff.withoutOrphans.gff',
'Lsat.1.v3.genePredict.AbInitio.gff3.reLocated.gff.masked.gff.withoutOrphans.gff',
'Lsat.1.v3.genePredict.geneWise_plusLink.gff3.reLocated.gff.masked.gff.withoutOrphans.gff',
#'Lsat.1.v3.MapPositions.gff3.reLocated.gff.masked.gff',
#'Lsat.1.v3.ncRNA.gff3.reLocated.gff.masked.gff',
'Lsat.1.v3.repeatElements.gff3.reLocated.gff.masked.gff',
'Lsat.1.v4.repeatElements.gff3.reLocated.gff.masked.gff',
#'Lsat.1.v4.scaffolds.gff3.reLocated.gff.rename.gff',
#'unMaker_scaffolds.gff.reLocated.gff'
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

open (FASTA, "Lsat.1.v4.ScaffoldsSequence.All.fasta_biggerthan_1000.fasta");

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
	open(GFF, $file) or die;
	while (my $line = <GFF>) {
		next if $line =~ /#/;
		@fields = split("\t", $line);
		push @{$GFF{$fields[0]}}, $line;

	}
	close GFF;

	print "$file load it into memory.\n";
	
}
print STDERR "GFF files load it\n";


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

system("mkdir Database_hasGBrowse") unless -d "Database_hasGBrowse";

system("mkdir Database_hasGBrowse/0/") unless -d "Database_hasGBrowse/0/";
system("mkdir Database_hasGBrowse/1/") unless -d "Database_hasGBrowse/1/";
system("mkdir Database_hasGBrowse/2/") unless -d "Database_hasGBrowse/2/";
system("mkdir Database_hasGBrowse/3/") unless -d "Database_hasGBrowse/3/";
system("mkdir Database_hasGBrowse/4/") unless -d "Database_hasGBrowse/4/";
system("mkdir Database_hasGBrowse/5/") unless -d "Database_hasGBrowse/5/";
system("mkdir Database_hasGBrowse/6/") unless -d "Database_hasGBrowse/6/";
system("mkdir Database_hasGBrowse/7/") unless -d "Database_hasGBrowse/7/";
system("mkdir Database_hasGBrowse/8/") unless -d "Database_hasGBrowse/8/";
system("mkdir Database_hasGBrowse/9/") unless -d "Database_hasGBrowse/9/";
system("mkdir Database_hasGBrowse/A/") unless -d "Database_hasGBrowse/A/";

foreach my $scaf (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;

#	if (exists($sequences{$scaf})) {

		@scafName = split("_", $scaf);
		open(OUT, ">Database_hasGBrowse/$scafName[4]/$scaf.gff") or die "Couldn't open file for $scaf.\n";
		foreach my $line (@{$GFF{$scaf}}) {
			print OUT $line;
		}
#		print OUT ">", $scaf, "\n";
#		print OUT $sequences{$scaf}, "\n";
		close (OUT);

#	} else {

#		print "Sequence not found for $scaf.\n";

#	} #end else

#	system("tar -czf Database/$scafName[4]/$scaf.gff.tar.gz Database/$scafName[4]/$scaf.gff");
}
print STDERR "done\n";

































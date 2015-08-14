#!/usr/bin/perl
use warnings;

################################
#                              #
# Global and storage variables #
#                              #
################################

my %GFF;

my @file = (
#'All_All_postMaker_exonerate_plusLink.gff3',
#'Lsat.1.v3.AffyContigs.gmap.onePath.gff3',
#'Lsat.1.v3.geneModels.gff3',
#'Lsat.1.v3.genePredict.AbInitio.gff3',
#'Lsat.1.v3.genePredict.geneWise_plusLink.gff3',
#'Lsat.1.v3.MapPositions.gff3',
#'Lsat.1.v3.repeatElements.gff3',
#'Lsat.1.v3.ncRNA.gff3',
#'Lsat.1.v3.scaffolds.gff3',
#'unMaker_scaffolds.gff'
'salinas.filtered_subreads_without_CCS.gmapouput.id80.cov90.good.gff',
'/home/sreyesch/GBrowse_v2_2/GBrowse_ready/Lsat.1.v4.geneModels.gff'
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

system("mkdir GFF_Database/") unless -d "GFF_Database/";

foreach my $scaf (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;

	open(OUT, ">GFF_Database/$scaf.gff") or die "GFF_Database/$scaf.gff\n";
	foreach my $line (@{$GFF{$scaf}}) {
		print OUT $line;
	}
	close (OUT);

}
print STDERR "done\n";

exit;


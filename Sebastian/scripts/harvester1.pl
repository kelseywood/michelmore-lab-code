#!/usr/bin/perl
use strict; use warnings;
use DataBrowser;
use FileHandle;

#die "this has already been run, and probably does not need running again";

my $wd = `pwd`; chomp $wd;
die "run only from jamboree" unless $wd eq "/home/ikorf/maker/jamboree";

##########################
# part 1: modify BGI GFF #
##########################

# determining the legal sequences
print STDERR "reading genome.fa...";
my %legal;
open(IN, "genome.fa") or die;
while (<IN>) {
	if (/^>(\S+)\s+(\S+)/) {
		my ($lsat, $len) = ($1, $2);
		$lsat =~ s/\./_/g;
		$legal{$lsat} = $len;
	}
}
close IN;
print STDERR "done\n";

# mapping from BGI names to Lsat names
print STDERR "mapping names...";
my %bgi_to_lsat;
open(IN, "BGI_to_Lsat.txt") or die;
while (<IN>) {
	chomp;
	my ($bgi, $lsat) = split;
	$lsat =~ s/\./_/g;
	($bgi) = $bgi =~ /^BGI_V3_(.+)/;
	$bgi_to_lsat{$bgi} = $lsat;
}
close IN;
print STDERR "done\n";

# read all GFF, edit as appropriate
print STDERR "munging GFFs...";
my @file = (
#	'BGI_GFF_Files/gene_prediction/lettuce.v3.scaffold.Glean.filter.final.gff',
	'BGI_GFF_Files/ncRNA/lettuce_S_miRNA.gff',
	'BGI_GFF_Files/ncRNA/lettuce_S_rRNA.gff',
	'BGI_GFF_Files/ncRNA/lettuce_S_snRNA.gff',
	'BGI_GFF_Files/ncRNA/lettuce_S_tRNA.gff',
	'BGI_GFF_Files/repeats/denovo.gff',
	'BGI_GFF_Files/repeats/lettuce.v3.scaffold.Proteinmask.annot.gff',
	'BGI_GFF_Files/repeats/lettuce.v3.scaffold.RepeatMasker.out.gff',
	'BGI_GFF_Files/repeats/lettuce.v3.scaffold.trf.dat.gff',
	'BGI_GFF_Files/RNASeqv3_annotation_1_1/longest/Lactuca_sativa.after_RNAseq.gff.noisoforms.gff_v3',
	'BGI_GFF_Files/gene_prediction/At.final.cds.filter.pep.solar.genewise.gff',
#	'BGI_GFF_Files/RNASeqv3_annotation_1_1/splicing/Lactuca_sativa.after_RNAseq.gff_v3',
);

my %GFF;
foreach my $file (@file) {
	my $limit;
	open(GFF, $file) or die;
	while (my $line = <GFF>) {
		next if $line =~ /#/;
		my ($bgi) = $line =~ /^(\S+)/;
		next if not defined $legal{$bgi_to_lsat{$bgi}};
		my $lsat = $bgi_to_lsat{$bgi};
		$line =~ s/$bgi/$lsat/g;
		$line =~ s/Transposon/transposable_element/;
		$line =~ s/TEprotein/transposable_element/;
		$line =~ s/TandemRepeat/region/;
		$line =~ s/UTR_3/region/;
		$line =~ s/UTR_5/region/;
		push @{$GFF{$lsat}}, $line;
#		last if $limit++ == 100; # testing
	}
	close GFF;
}
print STDERR "done\n";

# write GFF
my $i = 0;
my $n = keys %GFF;
my $m = int $n / 100;
print STDERR "writing $n files\n";
foreach my $scaf (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;
	open(OUT, ">BGI/$scaf.gff") or die;
	foreach my $line (@{$GFF{$scaf}}) {
		print OUT $line;
	}
	close OUT;
}
print STDERR "done\n";


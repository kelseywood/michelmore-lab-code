#!/opt/rocks/lib/perl5
############################################
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
#############################################
if( @ARGV < 4){

      print "usage: [Provide fasta file of your assembly] space [Provide the SNP file format from Bowtie] space [Provide output file name] space [Provide SNP Prefix name]\n";
      exit 0;
}


my $idx = Bio::Index::Fasta->new(
                                 -filename => "$ARGV[0].idx",
                                 -write_flag => 1
                                );
$idx->make_index("$ARGV[0]");


my $ambig = {
             'R' => { 'A' => 'G', 'G' => 'A'},
             'Y' => { 'C' => 'T', 'T' => 'C'},
             'S' =>	{ 'G' => 'C', 'C' => 'G'},
             'W' =>	{ 'A' => 'T', 'T' => 'A'},
             'K' =>	{ 'G' => 'T', 'T' => 'G'},
             'M' =>	{ 'A' => 'C', 'C' => 'A'},
            };

open(my $fh_ill, ">$ARGV[2]" ) or die;
print $fh_ill join(",", ('SNP_Name', 'Contig', 'Pos', 'Sequence') )."\n";

my $id = 0;

open(my $fh, $ARGV[1]) or die;
SNP:
while (<$fh>) {
    #handle output here
    my ($SNPid, $contig, $pos, $ref_base, $snp_base, undef ) = split(/\t/);
    print STDERR "$SNPid\t$contig\t$pos\t$ref_base\t$snp_base\n";

    $snp_base = uc($snp_base);
    $ref_base = uc($ref_base);

    #resolve ambig SNP
    my $snp_allele;
    if ( $snp_base =~ /[^ATCG]/) {
        $snp_allele = $ambig->{$snp_base}->{$ref_base};
        next SNP if ! $snp_allele; #means it is ambig without either base being a ref base
    }else {
        $snp_allele =  $snp_base;
    }

    my $snp_seq = create_snp_sequence($contig, $pos, $ref_base, $snp_allele);
    next SNP if ! $snp_seq;
    #my $snp_name = "G_hirs_snp_$id";
    print_illumina_line($SNPid, $contig, $pos, $snp_seq);
}
close($fh);
close($fh_ill);


sub create_snp_sequence {
    my ( $ref_acc, $ref_pos, $ref_allele, $snp_allele) = @_;
    my $flanking_span = 100;
    #get the seq obj
    my $seq = $idx->fetch($ref_acc);
    #some sanity checking - this should have been farther upstream
    if ( ($ref_pos-100 < 1) ||  ($ref_pos+100 > $seq->length) ) {
        return;
    }
    my $left_flank = $seq->subseq(($ref_pos-100), ($ref_pos-1));
    my $right_flank = $seq->subseq(($ref_pos+1), ($ref_pos+100));

    $left_flank = uc($left_flank);
    $right_flank = uc($right_flank);

    my $snp = "[$ref_allele/$snp_allele]";
    return $left_flank.$snp.$right_flank;
}

sub print_illumina_line {
    my ($SNPid, $contig, $pos, $sequence) = @_;
    #SNP_name sequence genome_build_ver chr coordinate source dbsnp_ver ploidy, species, Customer strand
    print $fh_ill join(",", ($SNPid, $contig, $pos, $sequence) )."\n";
    return;
}
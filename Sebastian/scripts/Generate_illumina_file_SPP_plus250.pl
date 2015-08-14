#!/opt/rocks/lib/perl5
############################################
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::Index::Fasta;
#############################################
if( @ARGV < 3){

      print "usage: [Provide fasta file of your assembly] space [Provide the SPP file] space [Provide output file name]\n";
      exit 0;
}


my $idx = Bio::Index::Fasta->new(
                                 -filename => "$ARGV[0].idx",
                                 -write_flag => 1
                                );
$idx->make_index("$ARGV[0]");


open(my $fh_ill, ">$ARGV[2]" ) or die;

my $id = 0;

open(my $fh, $ARGV[1]) or die;
SNP:
while (<$fh>) {
    #handle output here
    my ($contig, $startPos, $endPos, undef ) = split(/\t/);
    chomp($endPos);
    print STDERR "$contig\t$startPos\t$endPos\n";

    my $snp_seq = create_snp_sequence($contig, $startPos, $endPos);
    next SNP if ! $snp_seq;
    #my $snp_name = "G_hirs_snp_$id";
    print_illumina_line($contig, $startPos,$endPos, $snp_seq);
}
close($fh);
close($fh_ill);


sub create_snp_sequence {
    my ( $ref_acc, $ref_startPos, $ref_endPos) = @_;
    my $flanking_span = 16;
    #get the seq obj
    my $seq = $idx->fetch($ref_acc);
    #some sanity checking - this should have been farther upstream

    my $left_flank = "";
    my $right_flank = "";

    if ($ref_startPos-250 < 1) {
	$left_flank = $seq->subseq(1, ($ref_startPos));
    } else {
	$left_flank = $seq->subseq(($ref_startPos-250), ($ref_startPos));
    }

    if ($ref_endPos+250 > $seq->length) {
	$right_flank = $seq->subseq(($ref_startPos+1), ($seq->length));
    } else {
    	$right_flank = $seq->subseq(($ref_startPos+1), ($ref_endPos+250));
    }


    $left_flank = uc($left_flank);
    $right_flank = uc($right_flank);

    return $left_flank.$right_flank;
}

sub print_illumina_line {
    my ($contig, $startPos, $endPos, $sequence) = @_;
    #SNP_name sequence genome_build_ver chr coordinate source dbsnp_ver ploidy, species, Customer strand
    print $fh_ill ">", $contig, "-", $startPos, "_", $endPos, "\n", $sequence, "\n";
    return;
}
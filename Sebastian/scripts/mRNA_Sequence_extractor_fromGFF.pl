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

      print "usage: [Provide fasta file of your assembly] space [Provide the GFF file] space [Provide output file name]\n";
      exit 0;
}


my $idx = Bio::Index::Fasta->new(
                                 -filename => "$ARGV[0].idx",
                                 -write_flag => 1
                                );
$idx->make_index("$ARGV[0]");



open(TABLE, ">$ARGV[2].txt" ) or die;
print TABLE join("\t", ('ID', 'Scaffold', 'Start', 'End') )."\n";

open(FASTA, ">$ARGV[2].fasta" ) or die;


open(my $fh, $ARGV[1]) or die;

while (<$fh>) {
    #handle output here
    my ($scaffold, $source, $type, $start, $end, $score, $strand, $phase, $attribute, undef ) = split(/\t/);

    if($type eq "mRNA") {

        my @attributes = split(";", $attribute);

        my $ID = substr($attributes[0], 3, 100);

	print TABLE join("\t", $ID, $scaffold, $start, $end). "\n";

#####################################
#                                   #
#  If you want to get the sequence  #
#        activate this code         #
#                                  #
####################################
#        my $mRNA_seq = create_snp_sequence($scaffold, $start, $end);
#	print FASTA ">";
#	print FASTA join(" ", $ID, $scaffold, $start, $end). "\n";
#	print FASTA $mRNA_seq, "\n"

    } #end if

}
close($fh);
close(TABLE);
close(FASTA);


sub create_snp_sequence {
    my ($scaffold, $start, $end) = @_;

    my $seq = $idx->fetch($scaffold);

    my $mRNA = $seq->subseq($start, $end);

    return $mRNA;
}


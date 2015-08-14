#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_h $opt_f $opt_r $opt_t $opt_M $opt_G);
getopts('hf:r:t:M:G');

my $files = $opt_f if $opt_f;
my $reference = $opt_r if $opt_r;
my $cores = $opt_t if $opt_t;
my $mismatch = $opt_M if $opt_M;
my $gaps = $opt_G if $opt_G;

if( (defined($opt_h)) ) {

	print "Missing Arguments.\n";
	print "Please input all the arguments need it\n";
	print "BWA_bath_alginer.pl -h -f fileList -r reference.fasta -M 0.5 -G 1\n";
	print "-h output help\n";
	print "-f file with the list of fastq files to align\n";
	print "-r Name of the formatted reference fasta\n";
	print "-t number of cores to use\n";
	print "-M Mismmatch allowances\n";
	print "-G Number of gaps allowed\n";
	print "Help menu trigger\n\n";
	exit;

} #end if

if( !(defined($files)) || !(defined($reference)) || !(defined($cores))|| !(defined($mismatch)) || !(defined($gaps)) ) {

	print "Missing Arguments.\n";
	print "Please input all the arguments need it\n";
	print "BWA_bath_alginer.pl -h -f fileList -r reference.fasta -M 0.5 -G 1\n";
	print "-h output help\n";
	print "-f file with the list of fastq files to align\n";
	print "-r Name of the formatted reference fasta\n";
	print "-t number of cores to use\n";
	print "-M Mismmatch allowances\n";
	print "-G Number of gaps allowed\n";
	die "Missing Arguments";

} #end if

open(INPUT, $files);

open(LOG, ">$files.aligner.log.txt");
open(BAMS, ">$files.BAMlist.txt");

while(my $forward = <INPUT>) {

	chomp($forward);

	print "Processing $forward files.\n";
	print  LOG "Processing $forward files.\n";

	print LOG "bwa aln -t $cores -n $mismatch -o $gaps $reference $forward > $forward.sai\n";

	if( system("bwa aln -t $cores -n $mismatch -o $gaps $reference $forward > $forward.sai") != 0 ) {

		print "ERROR: Couldn't perform alignment for $forward.\n";

	} #end if

	print LOG "bwa samse $reference $forward.sai $forward > $forward.sam\n";

	if( system("bwa samse $reference $forward.sai $forward > $forward.sam") != 0 ) {

		print "ERROR: Couldn't create SAM file for $forward.\n";

	} #end if

	print LOG "samtools view -b -S -o $forward.bam $forward.sam\n\n";

	print BAMS "$forward.bam\n";

	if( system("samtools view -b -S -o $forward.bam $forward.sam") != 0 ) {

		print "ERROR: Couldn't transform SAM to BAM $forward.\n";

	} #end if

#	system("rm $forward.sam");

} #end while for input files

exit;


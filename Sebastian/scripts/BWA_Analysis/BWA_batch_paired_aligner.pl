#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_h $opt_f $opt_d $opt_r $opt_t $opt_M $opt_G);
getopts('hf:d:r:t:M:G');

my $files = $opt_f if $opt_f;
my $outDir = $opt_d if $opt_d;
my $reference = $opt_r if $opt_r;
my $cores = $opt_t if $opt_t;
my $gaps = $opt_G if $opt_G;
my $mismatch = $opt_M if $opt_M;

if( (defined($opt_h)) ) {

	print "Missing Arguments.\n";
	print "Please input all the arguments need it\n";
	print "BWA_bath_alginer.pl -h -f fileList -r reference.fasta -M 0.05 -G 1\n";
	print "-h output help\n";
	print "-f tab-delimited file with the list of fastq files to align (each line with two columns, forward and reverse sense)\n";
	print "-d output directory for alignment files\n";
	print "-r Name of the indexed reference (prefix)\n";
	print "-t number of cores to use\n";
	print "-M Mismmatch allowances\n";
	print "-G Number of gaps allowed\n";
	print "Help menu trigger\n\n";
	exit;

} #end if


if( !(defined($files)) || !(defined($reference)) || !(defined($cores))|| !(defined($mismatch)) || !(defined($gaps)) || !(defined($outDir)) ) {

	print "Missing Arguments.\n";
	print "Please input all the arguments need it\n";
	print "BWA_bath_alginer.pl -h -f fileList -r reference.fasta -M 0.05 -G 1\n";
	print "-f tab-delimited file with the list of fastq files to align (each line with two columns, forward and reverse sense)\n";
	print "-d output directory for alignment files\n";
	print "-r Name of the indexed reference (prefix)\n";
	print "-t number of cores to use\n";
	print "-M Mismmatch allowances\n";
	print "-G Number of gaps allowed\n";
	die "Missing Arguments";

} #end if

open(INPUT, $files);

system("mkdir $outDir") unless -d $outDir;

open(LOG, ">$outDir/$files.aligner.log.txt");

print LOG "Input fastq files: ", $files, "\n";
print LOG "Reference index: ", $reference, "\n";
print LOG "Output dir: ", $outDir, "\n";
print LOG "Number of threads to use: ", $cores, "\n";
print LOG "Mismatch percentage: ", $mismatch, "\n";
print LOG "Max gaps allow per read: ", $gaps, "\n";

open(BAMS, ">$outDir/$files.BAMlist.txt");
open(SAMS, ">$outDir/$files.SAMlist.txt");

while(my $line = <INPUT>) {

	chomp($line);

	my @files = split ("\t", $line);

	print "Processing $files[0] files.\n";
	print  LOG "Processing $files[0] files.\n";

	my @arrayName1 = split("/", $files[0]);

	my $filename1 = pop(@arrayName1);

	print LOG "bwa aln -t $cores -n $mismatch -o $gaps $reference $files[0] > $outDir/$files[0].sai\n";

	if( system("bwa aln -t $cores -n $mismatch -o $gaps $reference $files[0] > $outDir/$filename1.sai") != 0 ) {

		print "ERROR: Couldn't perform alignment for $files[0].\n";

	} #end if

	print "Processing $files[1] files.\n";
	print  LOG "Processing $files[1] files.\n";

	my @arrayName2 = split("/", $files[1]);

	my $filename2 = pop(@arrayName2);


	print LOG "bwa aln -t $cores -n $mismatch -o $gaps $reference $files[1] > $outDir/$files[1].sai\n";

	if( system("bwa aln -t $cores -n $mismatch -o $gaps $reference $files[1] > $outDir/$filename2.sai") != 0 ) {

		print "ERROR: Couldn't perform alignment for $files[1].\n";

	} #end if

	print LOG "bwa sampe $reference $outDir/$filename1.sai $outDir/$filename2.sai $files[0] $files[1] -f $outDir/$filename1.paired.sam 2> $outDir/$filename1.paired.sam.log\n";

	if( system("bwa sampe $reference $outDir/$filename1.sai $outDir/$filename2.sai $files[0] $files[1] -f $outDir/$filename1.paired.sam 2> $outDir/$filename1.paired.sam.log") != 0 ) {

		print "ERROR: Couldn't create SAM file for $files[1] and $files[1].\n";

	} #end if

	print SAMS "$files[1].paired.sam\n";

	print BAMS "$files[1].paired.bam\n";

	print LOG "samtools view -b -S -o $outDir/$filename1.paired.bam $outDir/$filename1.paired.sam\n\n";

	if( system("samtools view -b -S -o $outDir/$filename1.paired.bam $outDir/$filename1.paired.sam") != 0 ) {

		print "ERROR: Couldn't transform SAM to BAM for $files[1] and $files[1].\n";

	} #end if

} #end while for input files

exit;


#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_f $opt_o $opt_s);
getopts('hi:f:o:s:');

############################################
#                                          #
#           GFF_harvester_v12.pl            #
#                                          #
#         Sebastian Reyes Chin-Wo          #
#          sreyesch@ucdavis.edu            #
#            Genome Center                 #
#               UC Davis                   #
#                                          #
############################################

# Based on Ian Korf (ikorf@ucdavis.edu) harvester.pl script


###############################################
###############################################
#
# Diff v1 add split functionality for gff files
#
###############################################
###############################################


if(defined($opt_h)) {

	print "Detail help\n\n";

        print "Usage: GFF_harvester.pl -i gff_list.txt -o outpuDir  <-f fasta.file> <-s int>\n\n";
        print "Require arguments:\n";
        print "  -i <text.file>   text file with list of gff files to analyze (one per line) (can be generated from 'ls -1 > text.file'\n";
        print "  -o <directory>   output directory wher to print the data\n\n";
        print "Optional:\n";
        print "  -f <file.fasta>  fasta sequence to add into gff files for gff/fasta files\n";
        print "  -s <integer>     split gff randomly in directory with <int> files per directory\n\n";

        exit;

} #end if detail option was ask

if( !defined($opt_i) || !defined($opt_o) ) {

	print "Error: Missing required arguments.\n\n";
        print "Usage: GFF_harvester.pl -options -i gff_list.txt -o outputDir\n";
        print "For detail help please use GFF_harvester.pl -h\n";
        die "Missing arguments\n";

} #end if missing arguments

################################
#                              #
# Global and storage variables #
#                              #
################################

my $gff_list = $opt_i;

my $fastaFile = $opt_f;

my $outputDir = $opt_o;

my $filesperdir = "";

print "File containing the gff's to use : $gff_list\n";
print "Output directory : $outputDir\n";

if(defined($opt_s)) {$filesperdir = $opt_s; print "Files will be split in directory, files per directory: $filesperdir\n"}

open(GFFLIST, $gff_list);

my @file;

while(my $gffFile = <GFFLIST>) {

	if(substr($gffFile,0,1) ne "#") {

		chomp($gffFile);

		push (@file, $gffFile);

	} #end if commented line

} #end loading gff list

print STDERR scalar(@file), " files will load it as sources of information\n";

my %sequences = ();
my %GFF;


#######################
#                     #
#   Open fasta file   #
#                     #
#######################

if(defined($fastaFile)) {

	print "Fasta file use for retrieve sequence: $fastaFile\n";

	# Extract the sequences from the fasta file and
	# loaded into the sequences hash

	open (FASTA, $fastaFile);

	while (my $fastaLine = <FASTA>) {

	        my @header = split(" ", substr($fastaLine,1,200));

	        my $curSequence = <FASTA>;

	        chomp($curSequence);

	        $sequences {$header[0]} = $curSequence;

	} #end while

	close (FASTA);
	print STDERR "Fasta file load it\n";

} #end if loading fasta file

#######################
#                     #
#    Open GFF files   #
#                     #
#######################

# Extract the features from the GFF files and
# loaded it into the hash
# also create source-type list

my %sources;

foreach my $file (@file) {
	open(GFF, $file) or die;
	while (my $line = <GFF>) {
		next if $line =~ /^#/;
		my @fields = split("\t", $line);
		push @{$GFF{$fields[0]}}, $line;
		
		my $sourceType = $fields[1]."|".$fields[2];
		$sources{$sourceType} = "";

	}# End loading data
	close GFF;

	print "$file load it into memory.\n";

} #end of foreach gff file
print STDERR "GFF files load it\n";

#######################
#                     #
# Write Sources/types #
#                     #
#######################

print "Will generate sources/types file\n";

open(SOURCES, ">$gff_list.sources-types.txt");

foreach my $source (sort keys %sources) { $source =~ s/\|/\t/; print SOURCES $source, "\n" }

#######################
#                     #
#   Write GFF files   #
#                     #
#######################

# write GFF
my $i = 0;
my $n = keys %GFF;
my $m = int $n / 100;
if($m < 1) { $m = 1 }

print STDERR "writing $n files\n";

my $subDir = "";

if(defined($opt_s)) {

	$subDir = 1;

} #end if defined sub dir number

system("mkdir $outputDir") unless -d $outputDir;

system("mkdir $outputDir/$subDir") unless -d "$outputDir/$subDir";

if($n == 0) { die "Couldn't load any data into features array\n\n"}

foreach my $reference (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;

	if(defined($opt_s)) {
	        if($i > ($filesperdir*$subDir)) {
	                ++$subDir;
	                system("mkdir $outputDir/$subDir") unless -d "$outputDir/$subDir";
	        } #end if
        } #end if defined sub dir number

	if(defined($fastaFile)) {

		if (exists($sequences{$reference})) {

			open(OUT, ">$outputDir/$subDir/$reference.gff") or die "Couldn't open file for $reference.\n";
			foreach my $line (@{$GFF{$reference}}) {
				print OUT $line;
			}

			print OUT ">", $reference, "\n";
			print OUT $sequences{$reference}, "\n";

			close (OUT);

		} else {

			print "Sequence not found for $reference.\n";

		} #end else

        } else {

		open(OUT, ">$outputDir/$subDir/$reference.gff") or die "Couldn't open file for $reference.\n";
		foreach my $line (@{$GFF{$reference}}) {
			print OUT $line;
		}

        } #end else sequence or not sequence

} #end foreach references

print STDERR "done\n";

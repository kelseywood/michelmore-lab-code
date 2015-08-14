#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_i $opt_f $opt_o $opt_s);
getopts('hi:f:o:s');

############################################
#                                          #
#           GFF_harvester_v1.pl            #
#                                          #
#         Sebastian Reyes Chin-Wo          #
#          sreyesch@ucdavis.edu            #
#            Genome Center                 #
#               UC Davis                   #
#                                          #
############################################

# Based on Ian Korf (ikorf@ucdavis.edu) harvester.pl script

if(defined($opt_h)) {

	print "Detail help\n\n";

        print "Usage: GFF_harvester.pl -i gff_list.txt -o outpuDir  <-f fasta.file>\n\n";
        print "Require arguments:\n";
        print "  -i text file with list of gff files to analyze (one per line)\n";
        print "  -o output directory wher to print the data\n\n";
        print "Optional:\n";
        print "  -f fasta sequence to add into gff files for gff/fasta files\n";

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

open(GFFLIST, $gff_list);

my @file;

while(my $gffFile = <GFFLIST>) {

	if(substr($gffFile,0,1) ne "#") {

		chomp($gffFile);

		push (@file, $gffFile);

	} #end if commented line

} #end loading gff list

my %sequences = ();
my %GFF;


#######################
#                     #
#   Open fasta file   #
#                     #
#######################

if(defined($fastaFile)) {

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

foreach my $file (@file) {
	open(GFF, $file) or die;
	while (my $line = <GFF>) {
		next if $line =~ /#/;
		my @fields = split("\t", $line);
		push @{$GFF{$fields[0]}}, $line;

	}
	close GFF;

	print "$file load it into memory.\n";

} #end of foreach gff file
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

system("mkdir $outputDir") unless -d $outputDir;

foreach my $reference (keys %GFF) {
	print STDERR int(100 * $i / $n), " " if $i++ % $m == 0;

	if(defined($fastaFile)) {

		if (exists($sequences{$reference})) {

			open(OUT, ">$outputDir/$reference.gff") or die "Couldn't open file for $reference.\n";
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

		open(OUT, ">$outputDir/$reference.gff") or die "Couldn't open file for $reference.\n";
		foreach my $line (@{$GFF{$reference}}) {
			print OUT $line;
		}

        } #end else sequence or not sequence

} #end foreach references

print STDERR "done\n";

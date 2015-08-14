#!/usr/bin/perl
use strict; use warnings;

if ( !(defined($ARGV[0])) ) { die "Missing arguments: Please give a genbank file\n"}

open(GENBANK, $ARGV[0]);

my $root = $ARGV[0];

$root =~ s/\.[^.]+$//;

my $organism;
my @definition;
my $gi;
my %sequences;

my $numSequences = 0;

while(my $line = <GENBANK>) {

	chomp($line);

	if( $line =~ /DEFINITION/) {

		push(@definition, substr($line, 12, length($line)) );

		while(my $definitionLine = <GENBANK>) {

			chomp($definitionLine);

			if( $definitionLine =~ /ACCESSION/) {

				last

			} else {

				push(@definition, substr($definitionLine, 12, length($definitionLine)) )

			} #end else

		} #end definition line while

	} elsif ($line =~ /ORGANISM/) {

		$organism = substr($line, 12, length($line));

		$organism =~ s/ /_/g;

	} elsif ($line =~ /VERSION/) {

		my @giLine = split(" ", $line);

		$gi = pop(@giLine);

	} elsif ($line =~ /ORIGIN/) {

		my @sequence;

		while(my $sequenceline = <GENBANK>) {


			if( $sequenceline =~ /\/\//) {

				last;

			} else {

				$sequenceline =~ s/[0-9]|\s//g;

				push(@sequence, $sequenceline);

			} #end if is sequence

		} #end while sequence

#		print "Found sequence\n";
#		print "Description: ", join(" ", @definition), "\n";
#		print "Organism: $organism\n";
#		print "Organism: $gi\n";
#		print "Sequence: ", join("", @sequence), "\n\n";

		$sequences{$organism}{$gi} = [join(" ", @definition), join("", @sequence)];

		@definition = "";
		$organism = "";
		++$numSequences;

	} # end of line

} #end while

print "$numSequences were found\n";

foreach my $organism (keys %sequences) {

	print "Printing fasta file for $organism\n";

	open(FASTA, ">".$root.".".$organism.".FASTA");

	foreach my $gi ( keys %{$sequences{$organism}} ) {

		print FASTA ">".$sequences{$organism}{$gi}[0]."\n";
		print FASTA $sequences{$organism}{$gi}[1]."\n";

	} #end foreach gi

	close(FASTA)

} #end foreach organism


exit;
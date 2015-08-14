#!/usr/bin/perl
use strict;
use warnings;

## Script design to find families of orthologs from a blast out put (specifically an info2 file from Alex Kozik tcl_blast_parser)

if(@ARGV < 2) {

	print "Usage: Info1_PrefixAdder_Query-Subject.pl <info1 file> <Query prefix> <Subject prefix>\n";
	die "Missing arguments\n";

} #end if arguments missing

my $info1 = $ARGV[0];

my $queryPrefix = $ARGV[1];

my $subjectPrefix = $ARGV[2];

open(INFO1, $info1) or die "Can't open $info1\n";

open(RENAMED, ">$info1.renamed.info1");

while(my $line = <INFO1>) {

	chomp($line);

	my @fields = split("\t", $line);

	my $query = shift(@fields);

	print RENAMED $queryPrefix, $query, "\t";

	if( defined($fields[7]) ) {

		my $subject = shift(@fields);

		print RENAMED $subjectPrefix, $subject;

	} #end if hit

	foreach my $field (@fields) {

		print RENAMED "\t", $field;

	} #end printing rest

	print RENAMED "\n";

} #end while

close(INFO1);
close(RENAMED);

exit;

	




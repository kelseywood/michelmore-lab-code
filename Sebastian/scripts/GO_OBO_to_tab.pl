#!/usr/bin/perl

use strict;
use warnings;

# This script was designed the file obtain from http://www.geneontology.org/GO.downloads.ontology.shtml
# Using wget to download the generated file from OBO-Edit 2.3, not using the actual OBO XML file

if (@ARGV < 1) { 

	die "missing arguments\n";

} #end if



my $foundTerms = 0;

open (OBO, $ARGV[0]);
open (TEXT, ">$ARGV[1]");

print TEXT join("\t", ("id", "name", "namespace", "definition", "synonym", "is_obsolete", "alt_id", "is_a") ), "\n";

while (my $line = <OBO>) {

	if ($line =~ /\[Term\]/) {

		++$foundTerms;

		my $id = "";
		my $name = "";
		my $namespace = "";
		my $definition = "";
		my $synonym = "";
		my @alt_id;
		my @is_a;
		my $is_obsolete = "false";
			

		while (my $termLine = <OBO>) {

			chomp($termLine);

			if($termLine eq "") {last}

			my @info = split(": ", $termLine);

			if($info[0] =~ /id/) { $id = $info[1] }
			if($info[0] =~ /name/) { $name = $info[1] }
			if($info[0] =~ /namespace/) { $namespace = $info[1] }
			if($info[0] =~ /def/) { $definition = $info[1]; $definition =~ s/&quot;//g }
			if($info[0] =~ /synonym/) { $synonym = $info[1]; $synonym =~ s/(&quot;|&quot| EXACT \[\])//g }
			if($info[0] =~ /is_obsolete/) { $is_obsolete = $info[1] }

			if($info[0] =~ /alt_id/) { push(@alt_id, $info[1]) }
			if($info[0] =~ /is_a/) { push(@is_a, $info[1]) }

		} #end while reading term info

		print TEXT join("\t", ($id, $name, $namespace, $definition, $synonym, $is_obsolete, join(",", @alt_id), join(",", @is_a) ) ), "\n";

	} #end if term located

} #end while reading file

print "Found terms $foundTerms\n\n";

close(OBO);
close(TEXT);

exit;


#!/usr/bin/perl
use strict; use warnings;

my $dir = $ARGV[0];

opendir(DIR, $dir);

my @myFiles = "";

while ( defined( my $file = readdir(DIR) ) ) {

#	my @filename = split ("\.", $file);

	my $ext = substr($file,-4,4);

	print $ext, "\n";

	# skip '.' and '..'
	if ($ext eq "gff3") {

       		push( @myFiles, "$file");

      	}    # end if

}    # end while

closedir(DIR);

my @chrs = qw(0 1 2 3 4 5 6 7 8 9 A);

foreach my $file (@myFiles) {

	open(GFFINPUT, $file);
	open(GFFOUTPUT, ">$file.newName.gff3");

	while(my $line = <GFFINPUT>) {

		foreach my $chr (@chrs) {

			my $orgName = "Lsat.1.v3.g.$chr\.";

			my $newName = "Lsat_1_v3_g_$chr\_";

			$line =~ s/$orgName/$newName/g;

		} #end foreach

		print GFFOUTPUT $line;

	} #end while

	close(GFFOUTPUT);
	close(GFFINPUT);

} #end foreach

exit;




















#!/usr/bin/perl
use strict; use warnings;

if (@ARGV < 1) {

	print "Two arguments are need it please input them.\n";
	print "Number of folders to sort the files.\n";
	print "Folder that contains the files to sort.\n";
	print "Ussage: Random_sorter.pl <Number> <Dir>.\n";
	die;

} #end if

my @myFiles;
my $folders = $ARGV[0];
my $Dir = $ARGV[1];

# open and read DIR
opendir(DIR, $Dir);

while ( defined( my $file = readdir(DIR) ) ) {

	# skip '.' and '..'
	if ( $file ne "." && $file ne ".." ) {

       		push( @myFiles, $file );

        }    # end if
}    # end while

closedir(DIR);

for (my $i=1; $i<=$folders;++$i) {

	system("mkdir $Dir/$i") unless -d "$Dir/$i/";

} #end for

foreach my $file (@myFiles) {

	my $rand = rand(($folders));

	my $roundRand = sprintf("%.0f", $rand);

	if($roundRand < $rand) {

		++$roundRand;

	} #end if

#	print "\nprint into random $roundRand\n";

	system("mv $Dir/$file $Dir/$roundRand/$file");

} #end foreach



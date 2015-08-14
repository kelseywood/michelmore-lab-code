#!/usr/bin/perl

use warnings;

my $Dir = $ARGV[0];
my @SPPFiles = ();
my $SNPFile = $ARGV[1];

# open and read DIR
opendir(DIR, $Dir);

while ( defined( $file = readdir(DIR) ) ) {

	# skip '.' and '..'
	if ( $file ne "." && $file ne ".." ) {

       		push(@SPPFiles, $file);

        }    # end if
}    # end while

closedir(DIR);

print "\n";
print join("\n", @SPPFiles);
print "\n";

foreach $SPPfile (@SPPFiles) {

	print "$Dir/$SPPfile.$SNPFile.Results";
        print "\n";

	if(system("perl ../SNP_SPP_Comparer_v1.pl $SNPFile $SPPfile Results_$ARGV[2]_$SPPfile")) {
		print "error executing system command for comparing", $SNPfile, " and ", $SPPfile, "\n";
	}  #end if

} #end foreach

exit;
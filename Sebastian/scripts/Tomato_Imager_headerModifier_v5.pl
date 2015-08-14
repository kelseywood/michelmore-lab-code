#!/usr/bin/perl -w
use warnings;

my $Dir = $ARGV[0];
my @myFiles = ();

# open and read DIR
opendir(DIR, $Dir);

while ( defined( my $file = readdir(DIR) ) ) {

	# skip '.' and '..'
	if ( $file ne "." && $file ne ".." ) {

       		push( @myFiles, $file );

        }    # end if
}    # end while

closedir(DIR);

@myFiles = sort(@myFiles);

@directory_path = split("/", $Dir);

$directoryName = pop(@directory_path);

chomp($directoryName);

open(TABLE, ">Correlation_Table_Filenames_$directoryName.txt");

print TABLE $Dir, "\n";

foreach $file (@myFiles) {

	$image = $file;

	chop ($image);
	chop ($image);
	chop ($image);
	chop ($image);

	$section = substr($image, 0, 1);

	if($section eq "L") {

		$section = "t";

	} #end if

        $plot = uc(substr($image, 1, 30));

	$plot =~ s/ /-/;

	$newFile = join("", "2010_ca_", $plot, "_fruit_", $section,".jpg");

	if (system ("mv", "$Dir/$file", "$Dir/$newFile") != 0) {
		print "error executing system command for move ", $file, "\n";
	} #end if

        print TABLE $file , "\t", $newFile, "\n";

} #end foreach

close(TABLE);
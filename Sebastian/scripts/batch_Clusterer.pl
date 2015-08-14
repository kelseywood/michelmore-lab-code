#!/usr/bin/perl
use strict; use warnings;

my @chr = qw(1 2 3 4 5 6 7 8 9 A);

system("mkdir GFF_clustered") unless -d "GFF_clustred";

foreach my $chr (@chr) {
	print STDERR "processing chr $chr\n";
	system("mkdir GFF_clustered/$chr") unless -d "GFF_clustred/$chr";
	for (my $i = 0; $i < 10; $i++) {
		print STDERR "\tprocessing chunk $i\n";
		system("mkdir GFF/$chr/$i") unless -d "GFF/$chr/$i";

		# open and read DIR

		my $Dir = "GFF/$chr/$i/";

		system("mkdir GFF_clustered/$chr/$i") unless -d "GFF_clustred/$chr/$i";
	
		my $clusterdDir = "GFF_clustered/$chr/$i/";

		opendir(DIR, $Dir);

		my @myFiles = "";

		while ( defined( my $file = readdir(DIR) ) ) {

			# skip '.' and '..'
			if ( $file ne "." && $file ne ".." ) {

		       		push( @myFiles, $file );

	        	}    # end if

		}    # end while

		closedir(DIR);

		#Run clustering script

		foreach my $file (@myFiles) {


			if ($file ne "") {

				print STDERR "\tprocessing file $file\n";

				my $newFile = "$file.clustered.gff";

				if (system ("perl procgff.pl -x 20 -y 10 $Dir/$file > $clusterdDir/$newFile") != 0) {
					print "error executing system command for move ", $file, "\n";
				} #end if

			} #end if
		
		} #end foreach

	} #end for

} #end foreach

exit;

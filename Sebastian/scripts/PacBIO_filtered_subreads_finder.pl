#!/usr/bin/perl

######################################################
#                                                    #
#         PacBio Filtered SubReads Finder            #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to search the secondary analysis folders and retrieve the filtered files
# appending the folder id to the fasta file name

# You will provide the folder where the main secondary analysis SMRT cell folders are and the script will go in and find the filtered files

use strict;
use warnings;

if (@ARGV < 2) {

	print "Usage\n";
	print "PacBIO_filtered_subreads_finder.pl <Folder with SMRTCells> <Out Folder>\n";
	die "Missing arguments\n";

} #end if arguments

my $inFolder = $ARGV[0];
my $outFolder = $ARGV[1];

system("mkdir $outFolder") unless -d $outFolder;

opendir(DIR, $inFolder);
while (my $directory = readdir(DIR)) { 

	if($directory =~ /_secondary_analysis_/) {

		my @dirName = split("_", $directory);

		my $root = $dirName[3];

		opendir(SECONDANA, "$inFolder/$directory");
		while (my $subDir = readdir(SECONDANA)) { 

			if($subDir =~ /^data$/) {

#				print "$directory/$subDir it's a data folder of secondary analysis folder\n";

				opendir(DATA, "$inFolder/$directory/$subDir");
				while (my $subSubDir = readdir(DATA)) { 

					if ($subSubDir =~ /subreads/) {

						system("cp $inFolder/$directory/$subDir/$subSubDir $outFolder/$root.$subSubDir");

					} #end if a CSS file

				} #end if are CSS
				close(DATA);

			} #end is data folder

		} #end while Second analysis folder
		close(SECONDANA);
	
	} #end is a secondary analysis folder

} #end while for folder
closedir(DIR);

exit;

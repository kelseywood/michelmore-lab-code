#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

if(!defined($ARGV[1])) { print "No enough arguments provided, please provide a directory to look for maker outputs and a prefix for the output files\n"; die "No input\n"}

print STDERR "Extracting gff out of datastore\n";

my $gff_dir = $ARGV[1];
my $gffprefix = $ARGV[1];
my $threads = 

system("mkdir", $gff_dir) unless -d $gff_dir;

my @resultGFFs;
my $makerOutputDir = $ARGV[0];

my $pm = new Parallel::ForkManager($threads);
##### Parallelizing Conditions #####
$pm->run_on_finish( sub { my ($pid, $exit_code, $ident) = @_; print MAKERJOBSLOG "$ident finish with PID $pid, at ", getTime(), "\n"; } );
$pm->run_on_start( sub { my ($pid,$ident)=@_; print MAKERJOBSLOG "$ident started, pid: $pid, at ", getTime(),"\n"; } );
$pm->run_on_wait( sub {print STDERR "."}, 10 );

opendir(DIR, $makerOutputDir);
while (my $makerRun = readdir(DIR)) { find_gff($makerRun, $makerOutputDir) } #end foreach directory
closedir(DIR);

print STDERR scalar(@resultGFFs), " number of gff's were generated\n\n";

##### Generate complete GFF #####

print STDERR "Consolidating annotation info by source\n\n";

if (-e "$gffprefix*.gff") {

	system("rm $gffprefix*.gff");

} #end if output file exist

foreach my $gff (@resultGFFs) {

        my %annotation;

	open(GFF, $gff);

        while(my $gffLine = <GFF>) {

        	if($gffLine =~ /^[[:alnum:]]/) {

                	my @fields = split("\t", $gffLine);

                        if($fields[1] =~ /^[[:alnum:]]/) {

                        	push(@{$annotation{$fields[1]}}, $gffLine);

                        } #end if source starts with a letter

                } elsif ($gffLine =~ /^>/) {

                	last;

                } #end elsif line

        } #end while

        close(GFF);

        foreach my $source (keys(%annotation)) {

        	open(SOURCE, ">>$gffprefix.maker.output.$source.gff");

                foreach my $line (@{$annotation{$source}}) { print SOURCE $line }

                close(SOURCE);

        } #end foreach

} #end foreach gff files

exit;

####################################
#                                  #
#         Find GFF file sub        #
#                                  #
####################################

sub find_gff {

	my ($makerRun, $makerOutputDir) = @_;

	#Find output directories
	if($makerRun =~ /.output$/ && -d "$makerOutputDir/$makerRun") {

		print "$makerOutputDir/$makerRun\n";

		#open subdir and make new path
		opendir(RUNSRESULTS, "$makerOutputDir/$makerRun");

		while(my $datastore = readdir(RUNSRESULTS) ) {

			#Find the datastore
			if($datastore =~ /.datastore$/) {

				opendir(DATASTORE, "$makerOutputDir/$makerRun/$datastore");

				while(my $subDatastore = readdir(DATASTORE) ) {

					if( -d "$makerOutputDir/$makerRun/$datastore/$subDatastore") {

						opendir(SUBDATASTORE, "$makerOutputDir/$makerRun/$datastore/$subDatastore");

						while(my $subSubDatastore = readdir(SUBDATASTORE) ) {

							if( -d "$makerOutputDir/$makerRun/$datastore/$subDatastore/$subSubDatastore") {

								opendir(SUBSUBDATASTORE, "$makerOutputDir/$makerRun/$datastore/$subDatastore/$subSubDatastore");

								while(my $sequence = readdir(SUBSUBDATASTORE) ) {

									if( -d "$makerOutputDir/$makerRun/$datastore/$subDatastore/$subSubDatastore/$sequence") {

										opendir(SEQUENCERESULT, "$makerOutputDir/$makerRun/$datastore/$subDatastore/$subSubDatastore/$sequence");

										while( my $result = readdir(SEQUENCERESULT) ) {

											my $extension = substr($result,-3);

											if( $extension eq "gff") {

												system("cp $makerOutputDir/$makerRun/$datastore/$subDatastore/$subSubDatastore/$sequence/$result $gff_dir/$result");
												push(@resultGFFs, "$gff_dir/$result");

											} #end if for gff files

										} #end while reading result dir

									} #end if it is a sequence directory

								} #end while subsubdatastore

							} #end if subsubdatastored directory found

						} #end while subdatastore files

					} #end uf subdatore directory found

				} #end reading datastore directory

			} #end if is datastore

		} #end while maker run directory

	} #end if is an output directory

} #end find_gff sub


















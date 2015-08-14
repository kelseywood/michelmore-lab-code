#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
use warnings;
use Getopt::Long;

###################################################
#
#
#	GFF3 Genomic Sequence Extractor
#
#
#
###################################################

#####################################
#
# Diff from v1
# add flaking reginos
#

if (!defined($ARGV[2]) ) {

	print "Three argument is needed, please input it.\n";
	print ".\n";
	print "GFF_mRNA_integrity_checker.pl <BackUpGFF> <GFFToCheck> <OutPrefix>\n";
        exit;

} #end if

my %ExonPool;

open(BACKUPGFF, $ARGV[0]) or die "Can't open file $ARGV[0]\n";

while (my $line = <BACKUPGFF>) {

	chop($line);
        my @GFF_line = split ("\t", $line);

	if($GFF_line[2] eq "exon") {
		$ExonPool{$GFF_line[0]}{$GFF_line[3]}{$GFF_line[4]} = join("\t", ($GFF_line[0],$GFF_line[1],$GFF_line[2],$GFF_line[3],$GFF_line[4],$GFF_line[5],$GFF_line[6],$GFF_line[7]));
		$ExonPool{$GFF_line[0]}{$GFF_line[4]}{$GFF_line[3]} = join("\t", ($GFF_line[0],$GFF_line[1],$GFF_line[2],$GFF_line[3],$GFF_line[4],$GFF_line[5],$GFF_line[6],$GFF_line[7]));
	} #end if is an exon

} #end while

close(BACKUPGFF);

print STDERR "Exons from backup file loaded\n";

my %mRNAs;

open(GFFTOCHECK, $ARGV[1]) or die "Can't open $ARGV[1]\n";

while (my $line2 = <GFFTOCHECK>) {

	chop($line2);

        my @GFF_line2 = split ("\t", $line2);

	if ( ($GFF_line2[2] eq "exon") || ($GFF_line2[2] eq "CDS") ) {

		my @features = split (";", $GFF_line2[8]);

		my $name;

		foreach my $feature (@features) {
			if ( $feature =~ m/Parent/) {
				my @value = split("=", $feature);
				$name = $value[1];
			} #end else
		} #end foreach

		$mRNAs{$GFF_line2[0]}{$name}{$GFF_line2[2]}{$GFF_line2[3]}{$GFF_line2[4]} = $GFF_line2[8];

	} #end if it's an exon or CDS

} #end while

close(GFFTOCHECK);

print STDERR "Exons and CDSs loaded from file to check\n";

open(MISSING, ">$ARGV[2].missingExons.gff3") or die "Can't open output file\n";
open(STILLMISSING, ">$ARGV[2].stillMissingExons.gff3") or die "Can't open output file\n";

my $analyzedmRNA = 0;
my $addedFirstExons = 0;
my $addedMiddleExons = 0;
my $addedLastExons = 0;

foreach my $scaffold (keys %mRNAs) {

	foreach my $mRNA (keys %{$mRNAs{$scaffold}}) {

		my @CDSs = sort {$a <=> $b} keys %{$mRNAs{$scaffold}{$mRNA}{'CDS'}};

		my $firstCDS = shift(@CDSs);
		my $lastCDS = pop(@CDSs);

		
		if( !defined($mRNAs{$scaffold}{$mRNA}{'exon'}{$firstCDS}) ) {

			my %exonEnds;

			foreach my $key ( keys %{$mRNAs{$scaffold}{$mRNA}{'exon'}} ) { foreach my $key2 (keys %{$mRNAs{$scaffold}{$mRNA}{'exon'}{$key}} )  { $exonEnds{$key2} = "1" } }

			my @firstCDSthing = keys %{$mRNAs{$scaffold}{$mRNA}{'CDS'}{$firstCDS}};

			if (!defined($exonEnds{$firstCDSthing[0]}) ) {
				
				if( !defined($ExonPool{$scaffold}{$firstCDSthing[0]}) ) {

					print STILLMISSING "Exon ", $mRNA, " ", $firstCDS, " wasn't found on backup\n";

		                } elsif ( (!keys %{$ExonPool{$scaffold}{$firstCDSthing[0]}}) ) {

					print STILLMISSING "Exon ", $mRNA, " ", $firstCDS, " wasn't found on backup, empty hash\n";

				} else {

					my @Exonthing = keys %{$ExonPool{$scaffold}{$firstCDSthing[0]}};

					$mRNAs{$scaffold}{$mRNA}{'CDS'}{$firstCDS}{$firstCDSthing[0]} =~ s/cds/exon/g;
					print MISSING $ExonPool{$scaffold}{$firstCDSthing[0]}{$Exonthing[0]}, "\t", $mRNAs{$scaffold}{$mRNA}{'CDS'}{$firstCDS}{$firstCDSthing[0]}, "\n";

		                       	delete($mRNAs{$scaffold}{$mRNA}{'CDS'}{$firstCDS});

					++$addedFirstExons;

				} #end if exon present in the pool

			} #end if it's indeed missing

		} #end if missing first exon

		foreach my $CDS (@CDSs) {

                	if(!defined($mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS})) {next}

			my @CDSthing = keys %{$mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS}};

			if (!defined($mRNAs{$scaffold}{$mRNA}{'exon'}{$CDS}{$CDSthing[0]}) ) {

				if(!defined($ExonPool{$scaffold}{$CDS}{$CDSthing[0]})) {

					print STILLMISSING "Exon ", $mRNA, " ", $CDS, " wasn't found on backup\n";

				} else {

#                                	print $mRNA, "\t", $CDS, "\t", $CDSthing[0], "\t", keys %{$mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS}}, "\n";
					$mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS}{$CDSthing[0]} =~ s/cds/exon/g;
					print MISSING $ExonPool{$scaffold}{$CDS}{$CDSthing[0]}, "\t", $mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS}{$CDSthing[0]}, "\n";

                                        delete($mRNAs{$scaffold}{$mRNA}{'CDS'}{$CDS});

                                        ++$addedMiddleExons;

				} #end if is in exon pool
			}

		} #end foreach CDS

		if(defined($lastCDS) ) {

			my @lastCDSthing = keys %{$mRNAs{$scaffold}{$mRNA}{'CDS'}{$lastCDS}};

			if( !defined($mRNAs{$scaffold}{$mRNA}{'exon'}{$lastCDS}) ) {

				if( !defined($ExonPool{$scaffold}{$lastCDS}) ) {

					print STILLMISSING "Exon ", $mRNA, " ", $lastCDS, " wasn't found on backup\n";

		                } elsif ( (!keys %{$ExonPool{$scaffold}{$lastCDS}}) ) {

					print STILLMISSING "Exon ", $mRNA, " ", $lastCDS, " wasn't found on backup, empty hash\n";

				} else {

					if(!defined($mRNAs{$scaffold}{$mRNA}{'CDS'}{$lastCDS})) {next}

					my @Exonthing = keys %{$ExonPool{$scaffold}{$lastCDS}};

	#				print @lastCDSthing, "\n";

					$mRNAs{$scaffold}{$mRNA}{'CDS'}{$lastCDS}{$lastCDSthing[0]} =~ s/cds/exon/g;
					print MISSING $ExonPool{$scaffold}{$lastCDS}{$Exonthing[0]}, "\t", $mRNAs{$scaffold}{$mRNA}{'CDS'}{$lastCDS}{$lastCDSthing[0]}, "\n";

		                       	delete($mRNAs{$scaffold}{$mRNA}{'CDS'}{$lastCDS});

					++$addedLastExons;

				} #end if exon present in the pool


			} #end if missing last exon exon

		} #end if defined last exon

	} #end foreach mRNA

} #end printing scaffold

close(MISSING);

print "Added first exons $addedFirstExons\n";
print "Added middle exons $addedMiddleExons\n";
print "Added last exons $addedLastExons\n";
print "Total added exons ", $addedFirstExons + $addedMiddleExons + $addedLastExons, "\n";
exit;

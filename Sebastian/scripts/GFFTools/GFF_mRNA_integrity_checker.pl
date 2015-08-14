#!/usr/bin/perl
use lib "/opt/rocks/lib/perl5/site_perl/5.10.1/";
use strict;
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

if (!defined($ARGV[0]) ) {

	print "One argument is needed, please input it.\n";
	print "GFF File to verify integrity of gene model structurs.\n";
	print "GFF_mRNA_integrity_checker.pl <GFFFile>\n";
        exit;

} #end if

my %mRNAs;

my $countmRNAs = 0;

open(GFF, $ARGV[0]) or die "Can't open file $ARGV[0]\n";

while (my $line = <GFF>) {

	chop($line);

        my @GFF_line = split ("\t", $line);

	my @features = split (";", $GFF_line[8]);

	foreach my $feature (@features) {

		my @value = split("=", $feature);

	        if ( ($value[0] eq "ID") && ($GFF_line[2] eq "mRNA") ) {

			$mRNAs{$value[1]}{'data'} = [$GFF_line[0],$GFF_line[1],$GFF_line[6],$GFF_line[8],$GFF_line[3],$GFF_line[4]];

			++$countmRNAs;

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "UTR_5") ) {

			$mRNAs{$value[1]}{'UTR_5'}{$GFF_line[3]} = $GFF_line[4];

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "UTR_3") ) {

			$mRNAs{$value[1]}{'UTR_3'}{$GFF_line[3]} = $GFF_line[4];

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "CDS") ) {

			$mRNAs{$value[1]}{'CDS'}{$GFF_line[3]} = $GFF_line[4];

	        } elsif ( ($value[0] eq "Parent") && ($GFF_line[2] eq "exon") ) {

			$mRNAs{$value[1]}{'exon'}{$GFF_line[3]} = $GFF_line[4];

	        } #end else

	} #end foreach

} #end while

print "$countmRNAs mRNA's were loaded and will be use for sequence extraction\n\n";

close(GFF);


open(ERRORS, ">$ARGV[0].errors.txt");

my $analyzedmRNA = 0;

foreach my $mRNA (sort {$a <=> $b} keys %mRNAs) {

	if( ($analyzedmRNA/1000) =~ m/^\d+$/) {

		print "Done with $analyzedmRNA mRNA's of $countmRNAs\n";

	} #end if integer

	++$analyzedmRNA;

	my @mRNA_info = @{$mRNAs{$mRNA}{'data'}};

	delete $mRNAs{$mRNA}{'data'};

        my @genomic;
        my @transcript;
        my @CDS;
	my $codingSequence;
        my $proteinSequence;

	###  sequence is on plus ###
        if ($mRNA_info[2] eq "+") {

#		print ERRORS $mRNA, "\n";

		if(!defined($mRNAs{$mRNA}{'CDS'})) { print "Error: $mRNA don't have CDS\n"}
		if(!defined($mRNAs{$mRNA}{'exon'})) { print "Error: $mRNA don't have Exons\n"}

		my @CDSStarts = sort {$a<=>$b} keys %{$mRNAs{$mRNA}{'CDS'}} ;
		my @ExonStarts = sort {$a<=>$b} keys %{$mRNAs{$mRNA}{'exon'}} ;

		if (scalar(@CDSStarts) > scalar(@ExonStarts) ) {
#			print ERRORS "$mRNA contain more CDS's than exons, please check GFF file\n";
		} #end if

		#Retrieve 5' UTR if present

		my $badCount = 0;

		#Get spliced UTR exons
		while($mRNAs{$mRNA}{'exon'}{$ExonStarts[0]} < $mRNAs{$mRNA}{'CDS'}{$CDSStarts[0]}) {

			if(!defined($ExonStarts[0])) {last}

			my $firstUTR = shift(@ExonStarts);

		} #end while

		

		# Get first CDS
		#Get UTR piece attach to first CDS
#		if( ($ExonStarts[0] < $CDSStarts[0]) && ($mRNAs{$mRNA}{'CDS'}{$CDSStarts[0]} == $mRNAs{$mRNA}{'exon'}{$ExonStarts[0]}) ) {
#
#			my $firstCDS = shift(@CDSStarts);
#		        my $firstExon = shift(@ExonStarts);
#
#		} #end of printing 5' UTR

                my $CDSIndex = 0;
                my $ExonIndex = 0;

		# Get information out of CDS's
                while ($CDSIndex < scalar(@CDSStarts) ) {

                       		my $start_start;
                                my $start_end;
                                my $end_start;
                                my $end_end;

                                #Start of the gene against start of the evidence
                                if( ($CDSStarts[$CDSIndex]) > $ExonStarts[$ExonIndex]) { $start_start = "before" }
                                elsif ( ($CDSStarts[$CDSIndex]) < $ExonStarts[$ExonIndex]) {    $start_start = "after" }
                                else { $start_start = "match" } #end else for start-start comparison

                                #Start of the gene against end of the evidence
                                if( ($CDSStarts[$CDSIndex]) > $mRNAs{$mRNA}{'exon'}{$ExonStarts[$ExonIndex]}) { $start_end = "before" }
                                elsif ( ($CDSStarts[$CDSIndex]) < $mRNAs{$mRNA}{'exon'}{$ExonStarts[$ExonIndex]}) {    $start_end = "after" }
                                else { $start_end = "match" } #end else for start-end comparison

                                #End of the gene against start of the evidence
                                if( ($mRNAs{$mRNA}{'CDS'}{$CDSStarts[$CDSIndex]}) > $ExonStarts[$ExonIndex]) { $end_start = "before" }
                                elsif ( ($mRNAs{$mRNA}{'CDS'}{$CDSStarts[$CDSIndex]}) < $ExonStarts[$ExonIndex]) {    $end_start = "after" }
                                else { $end_start = "match" } #end else for end-start comparison

                                #End of the gene against end of the evidence
                                if( ($mRNAs{$mRNA}{'CDS'}{$CDSStarts[$CDSIndex]}) > $mRNAs{$mRNA}{'exon'}{$ExonStarts[$ExonIndex]}) { $end_end = "before" }
                                elsif ( ($mRNAs{$mRNA}{'CDS'}{$CDSStarts[$CDSIndex]}) < $mRNAs{$mRNA}{'exon'}{$ExonStarts[$ExonIndex]}) {    $end_end = "after" }
                                else { $end_end = "match" } #end else for end-end comparison


					 if ( ($CDSStarts[$CDSIndex] - $mRNAs{$mRNA}{'CDS'}{$CDSStarts[$CDSIndex]}) == 0 ) {

						#Status
						if($CDSIndex eq $CDSIndex) {++$ExonIndex}

						#Walking
	                                        ++$CDSIndex;

	                                #Exon evidence before exon of gene
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "after") && ($end_end eq "after") ) {

						#Example
						#Gene model |---|
						#Evidence           |---|
                                                print ERRORS "Error: state 1 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#walking
	                                        ++$CDSIndex;

	                                #Exon evidence before exon of gene
	                                } elsif ( ($start_start eq "before") && ($start_end eq "before") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model        |---|
						#Evidence   |---|
                                                print ERRORS "Error: state 2 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$ExonIndex;

	                                #Overlaping exons (evidence first)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model   |---|
						#Evidence   |---|
                                                print ERRORS "Error: state 3 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$ExonIndex;

	                                #Overlaping exons (gene first)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#Gene model |---|
						#Evidence     |---|
                                                print ERRORS "Error: state 4 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
	                                        ++$CDSIndex;

	                                #Partial match (end matches) (evidence shorter)
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#Gene model |---|
						#Evidence    |--|
                                                print ERRORS "Error: state 5 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
	                                        ++$CDSIndex;
	                                        ++$ExonIndex;

	                                #Partial match (end matches) (gene shorter)
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#Gene model  |--|
						#Evidence   |---|
						if($CDSIndex != 0) {
	                                                print ERRORS "Error: state 6 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#							print ERRORS "        $start_start $start_end $end_start $end_end\n";
						} # if is not the first CDS
						#Walking
	                                        ++$CDSIndex;
	                                        ++$ExonIndex;

	                                #Partial match (start matches) (evidence shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model |---|
						#Evidence   |--| 
						if($CDSIndex < scalar(@CDSStarts)) {
	                                                print ERRORS "Error: state 7 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#							print ERRORS "        $start_start $start_end $end_start $end_end\n";
						} # if is not the last CDS
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$ExonIndex;

	                                #Partial match (start matches) (gene shorter)
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#Gene model |--|
						#Evidence   |---|

						if($CDSIndex+1 != scalar(@CDSStarts) ) {
		                                        print ERRORS "Error: state 8 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#							print ERRORS "        $start_start $start_end $end_start $end_end\n";
						} #end if is not last CDS

						#Walking
	                                        ++$CDSIndex;

	                                #Evidence exon within gene exon
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model |---|
						#Evidence    |-|
                                                print ERRORS "Error: state 9 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$ExonIndex;

	                                #Gene exon within evidence exon
	                                } elsif ( ($start_start eq "before") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "after") ) {

						#Example
						#Gene model  |-|
						#Evidence   |---|

						if(scalar(keys %{$mRNAs{$mRNA}{'CDS'}}) != 1) {
	                                                print ERRORS "Error: state 10 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#							print ERRORS "        $start_start $start_end $end_start $end_end\n";
						} #end if
						#Walking
	                                        ++$CDSIndex;

	                                #Overlapping end-start
	                                } elsif ( ($start_start eq "before") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model     |---|
						#Evidence   |---|
                                                print ERRORS "Error: state 11 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$ExonIndex;

	                                #Overlapping start-end
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "after") ) {

						#Example
						#Gene model |---|
						#Evidence       |---|
                                                print ERRORS "Error: state 12 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
#						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
	                                        ++$CDSIndex;

	                                #
	                                } elsif ( ($start_start eq "match") && ($start_end eq "match") && ($end_start eq "before") && ($end_end eq "before") ) {

						#Example
						#Gene model |---|
						#Evidence   |
                                                print ERRORS "Error: state 13 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
						#Walking
	                                        ++$ExonIndex;

	                                #
	                                } elsif ( ($start_start eq "after") && ($start_end eq "after") && ($end_start eq "match") && ($end_end eq "match") ) {

						#Example
						#Gene model     |
						#Evidence   |---|
                                                print ERRORS "Error: state 14 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]\n";
						#Walking
						if($ExonIndex >= scalar(@ExonStarts) ) { ++$CDSIndex }
	                                        ++$CDSIndex;
	                                        ++$ExonIndex;

	                                #Perfect match
	                                } elsif ( ($start_start eq "match") && ($start_end eq "after") && ($end_start eq "before") && ($end_end eq "match") ) {

						#Example
						#Gene model |---|
						#Evidence   |---|
#                                                print ERRORS "Error: state 13 on CDS $CDSIndex on $mRNA\n";
						#Walking
	                                        ++$CDSIndex;
	                                        ++$ExonIndex;

					#Anything else
	                                } else {

                                                print ERRORS "Error: state 15 on CDS $CDSIndex on $mRNA for CDS $CDSStarts[$CDSIndex] with exon $ExonStarts[$ExonIndex]";
						#Status
#	                                        print ERRORS "None of the comparisons return TRUE for $mRNA, CDS $CDSIndex, exon $ExonIndex.\n";
						print ERRORS "        $start_start $start_end $end_start $end_end\n";
						#Walking
	                                        ++$CDSIndex;
	                                        ++$ExonIndex;

	                                } # end else check status


                } #end print CDS

		my $lastExon = pop(@ExonStarts);
		my $lastCDS = pop(@CDSStarts);

        ###  sequence is on minus ###
        } else {

		next;

        } #end print plus and minus strands

} #end printing mRNAs

exit;

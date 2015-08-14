###########################################################
###########################################################
#                                                         #
#                                                         #
#                   High_Troughput_MAQ                    #
#                                                         #
#                 Sebastian Reyes Chin-Wo                 #
#                Seed Biotechnology Center                #
#            University of California at Davis            #
#                                                         #
#                  sreyesch@ucdavis.edu                   #
#                                                         #
###########################################################
###########################################################

#Script design to run the MAQ application using has input raw reads files and a fasta file with the reference sequence.
#It can be use with unpaired or paired reads.

if (@ARGV < 6) {

	print "Program usage: Six inputs parameters are needed. Please enter them\n";
	print "First and Second input the files for paired end reads in case to have it, if not gave a 0\n";
	print "Third input is the file for non-paired reads, if not input 0\n";
	print "Fourth input is the Fasta file with the reference sequence\n";
	print "Fifth input is the experiment ID for the results file\n";
	print "Sixth input o is the maximum number of mismatches that can be found\n";
	exit 0;
}


#Initialization of global counting and storage variables
my $readIndex = 1;
my $readIndex2 = 1;
my $fileRecord = 1;
my $fileNumber = 1;
my $recordsPerFile = 2000000;
my @mapFiles = ();
my $experID = $ARGV[4];


#Generate output filenames

$bfaFile = "$ARGV[3].bfa";
$fastQ1 = "Temporary_$experID.1.fastQ";
$fastQ2 = "Temporary_$experID.2.fastQ";
$fastQ = "Temporary_$experID.fastQ";

$matchPaired = "Temporary_$experID.PAIRED";
$matchUnPaired = "Temporary_$experID.UNPAIRED";

$mapfile = "Results_$experID.map";
$assemblefile = "Results_$experID.assemble";
$checkfile = "Results_$experID.check";
$winfile = "Results_$experID.win";
$viewfile = "Results_$experID.view";
$snpfile = "Results_$experID.snp";
$pileupfile = "Results_$experID.pileup";
$mapViewfile = "Results_$experID.mapView";


#Convertion of reference fasta file into maq .bfa input file
if(system("./maq", "fasta2bfa", $ARGV[3], $bfaFile) != 0) {
	print "error executing system command for fasta2bfa on line", $line, "\n";
}  #end if


if ($ARGV[0] ne 0 && $ARGV[1] ne 0) {

	######################################
	######################################
	#                                    #
	#      Procedure for Paired Reads    #
	#                                    #
	######################################
	######################################


	#open reads file for Paired end reads first sense
	open(INPUT1, $ARGV[0]) || die "Cannot open file \"$ARGV[0]\""; # The file Containing reads sequence

	#Open temporary fastQ file for reads data
	open(RESULTS, ">$fastQ1") || die "Cannot open temporary fastQ file";
	open(RESULTS2, ">$fastQ2") || die "Cannot open temporary fastQ file";

	while ($read = <INPUT1>) {

		#Generate Unique readID for specific pair
		$readID = sprintf("%010d", $readIndex);

		#Check number of reads in the temporary file
		if($fileRecord < $recordsPerFile) {

			#Less than the stablish number, continue printing

			#Generate quality line
	                $lengthSequence = length ($read);
	                @quality = ();
	       	        for($i=1; $i < $lengthSequence; ++$i) {
	                        push (@quality, "b");
			} #end for

	       	        $sequence = $read;
	                chop($sequence);

			#Print read in fastQ format
	                print RESULTS "\@", $experID,"_IGA_",$readID,"\\1\n";
	               	print RESULTS $sequence, "\n";
        	        print RESULTS "+", $experID,"_IGA_",$readID,"\\1\n";
	                print RESULTS @quality, "\n";

		} else {

			print "Finish file $fileNumber for first sense of paired reads\n";

			#If the file reach the number of reads allow per file, read the data for the other sense

			#Open reads file for Paired end reads seconde sense
			open(INPUT2, $ARGV[1]) || die "Cannot open file \"$ARGV[0]\""; # The file Containing reads sequence

			$fileRecord = 1;
			$readIndex2 = 0;

			while ($read2 = <INPUT2>) {

				#Generate Unique readID for specific pair
				$readID2 = sprintf("%010d", $readIndex2);

#				print $fileRecord2,"\n";

				if($fileRecord2 > (($fileNumber - 1) * $recordsPerFile)) {

#					print "Working on row $fileRecord2, with $fileRecord\n";

					if($fileRecord < $recordsPerFile) {

#						print "Got into printing line for seconde sense line $read2\n";

						#Less than the stablish number, continue printing

						#Generate quality line
						$lengthSequence = length ($read2);
						@quality2 = ();
			                	for($i=1; $i < $lengthSequence; ++$i) {
		        	                	push (@quality2, "b");
			                	} #end for

						$sequence2 = $read2;
						chop($sequence2);

						#Print read in fastQ format
				                print RESULTS2 "\@", $experID,"_IGA_",$readID2,"\\2\n";
	       		 		        print RESULTS2 $sequence2, "\n";
						print RESULTS2 "+", $experID,"_IGA_",$readID2,"\\2\n";
				                print RESULTS2 @quality2, "\n";

						++$fileRecord;

					} else {

						print "Finish file $fileNumber for second sense of paired reads\n";

						#If file for the second sense reach the number of reads start the MAQ commands

						close(RESULTS);
						close(RESULTS2);

						# Convertion of fastQ files into .bfq maq input files
						if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ1, "$fastQ1.bfq") != 0) {
							print "error executing system command for fastq2bfq on line", $line, "\n";
						}  #end if
						if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ2, "$fastQ2.bfq") != 0) {
							print "error executing system command for fastq2bfq on line", $line, "\n";
		                                }  #end if

						# Align .bfq files with read data to the .bfa file with the reference sequence
						if(system("./maq", "map", "-C", "512", "-n", $ARGV[5], "$matchPaired$fileNumber.match", $bfaFile, "$fastQ1\@1.bfq", "[$fastQ2\@1.bfq]") != 0) {
							print "error executing system command match on line", $line, "\n";
						}  #end if

						# Load array with the .match file names
						push(@mapFiles, "$matchPaired$fileNumber.match");

						# Open a new temporary file for fastQ files for both senses and print the firs read
						open(RESULTS, ">$fastQ1") || die "Cannot open file results file";
						print RESULTS "\@", $curVar,"_IGA_",$readID,"\\1\n";
						print RESULTS $sequence, "\n";
						print RESULTS "+", $curVar,"_IGA_",$readID,"\\1\n";
						print RESULTS @quality0, "\n";

						open (RESULTS2, ">$fastQ2") || die "Cannot open file results file";
	                	                print RESULTS2 "\@", $experID,"_IGA_",$readID2,"\\2\n";
	        	                        print RESULTS2 $sequence2, "\n";
		                                print RESULTS2 "+", $experID,"_IGA_",$readID2,"\\2\n";
		                                print RESULTS2 @quality2, "\n";

						# Restore variables and count the file process
						++$fileNumber;
						$fileRecord = 0;
						$fileRecord2 = 0;

						last;

					} #end else

				} #end if

                        # Counting variables
                        ++$readIndex2;
                        ++$fileRecord2;

			} #end while

		} #end else

	# Counting variables
	++$readIndex;
	++$fileRecord;

	} #end while

	close(RESULTS);
	close(RESULTS2);


        # Convertion of fastQ files into .bfq maq input files
        if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ1, "$fastQ1.bfq") != 0) {
        	print "error executing system command for fastq2bfq on line", $line, "\n";
        }  #end if
        if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ2, "$fastQ2.bfq") != 0) {
        	print "error executing system command for fastq2bfq on line", $line, "\n";
        }  #end if

        # Align .bfq files with read data to the .bfa file with the reference sequence
        if(system("./maq", "map", "-C", "512", "-n", $ARGV[5], "$matchPaired$fileNumber.match", $bfaFile, "$fastQ1\@1.bfq", "[$fastQ2\@1.bfq]") != 0) {
        	print "error executing system command for match on line", $line, "\n";
        }  #end if

        # Load array with the .match file names
        push(@mapFiles, "$matchPaired$fileNumber.match");

} else {

	print "Non Paired-End reads were input\n"

} #end else

$fileRecord = 1;

if ($ARGV[2] ne 0) {

        ######################################
        ######################################
	#                                    #
        #    Procedure for UnPaired Reads    #
	#                                    #
        ######################################
        ######################################

	#open Raw Reads file for Paired end reads first sense
	open(INPUT3, $ARGV[2]) || die "Cannot open file \"$ARGV[0]\""; # The file Containing reads sequence

	open(RESULTS, ">$fastQ") || die "Cannot open temporary fastQ file";

	while ($read = <INPUT3>) {


                #Generate Unique readID for read
		$readID = sprintf("%010d", $readIndex);

                #Check number of reads in the temporary file
		if($fileRecord <= $recordsPerFile) {

                        #Less than the stablish number, continue printing

                        #Generate quality line
	                $lengthSequence = length ($read);
	                @quality = ();
        	        for($i=1; $i < $lengthSequence; ++$i) {
	                        push (@quality, "I");
			} #end for

        	        $sequence = $read;
	                chop($sequence);

			# Print read in fastQ format
	                print RESULTS "\@", $experID,"_IGA_",$readID,"\n";
                	print RESULTS $sequence, "\n";
        	        print RESULTS "+", $experID,"_IGA_",$readID,"\n";
	                print RESULTS @quality, "\n";

		} else {

			#If reach the number of reads per file start the MAQ commands

			print "File ", $fileNumber, " ready\n";

			close(RESULTS);

			# Convert the fastQ file into .bfq maq input
			if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ, "$fastQ.bfq") != 0) {
				print "error executing system command for fastq2bfq on line", $line, "\n";
			}  #end if

			# Run the alignment trough match of the reads against the reference sequence
			if(system("./maq", "map", "-C", "512", "-n", $ARGV[5], "$matchUnPaired$fileNumber.match", $bfaFile, "$fastQ\@1.bfq") != 0) {
				print "error executing system command for match on line", $line, "\n";
			}  #end if

			# Load array with the .match file name
			push(@mapFiles, "$matchUnPaired$fileNumber.match");

			# Open new temporary fastQ file
			open(RESULTS, ">$fastQ") || die "Cannot open file results file";

			# Print read in fastQ format
			print RESULTS "\@", $curVar,"_IGA_",$readID,"\n";
			print RESULTS $sequence, "\n";
			print RESULTS "+", $curVar,"_IGA_",$readID,"\n";
			print RESULTS @quality0, "\n";

			++$fileNumber;

			$fileRecord = 1;

		} #end else

	++$readIndex;
	++$fileRecord;

	} #end while

        # Convert the fastQ file into .bfq maq input
        if(system("./maq", "fastq2bfq", "-n", $recordsPerFile, $fastQ, "$fastQ.bfq") != 0) {
        	print "error executing system command for fastq2bfq on line", $line, "\n";
        }  #end if

        # Run the alignment trough match of the reads against the reference sequence
        if(system("./maq", "map", "-C", "512", "-n", $ARGV[5], "$matchUnPaired$fileNumber.match", $bfaFile, "$fastQ\@1.bfq") != 0) {
        	print "error executing system command for match on line", $line, "\n";
        }  #end if

        # Load array with the .match file name
        push(@mapFiles, "$matchUnPaired$fileNumber.match");

} else {

	print "Non Non-Paired-End reads were input"

} #end else


$MapList = join(" ", @mapFiles);

print "map Files generated ", $MapList, "\n";

if(system("./maq mapmerge $mapfile $MapList") != 0) {
	print "error executing system command for MapMerge", $line, "\n";
}  #end if

if(system("./maq pileup -p $bfaFile $mapfile > $pileupfile") != 0) {
        print "error executing system command for cns2view on line", $line, "\n";
}  #end if

if(system("./maq mapview -p $mapfile > $mapViewfile") != 0) {
        print "error executing system command for cns2view on line", $line, "\n";
}  #end if

if(system("./maq", "assemble", "-r", "0.01", "-m", "5", "-Q", "200", "-q", "0", $assemblefile, $bfaFile, $mapfile) != 0) {
	print "error executing system command for assemble on line", $line, "\n";
}  #end if

if(system("./maq mapcheck -q 0 $bfaFile $mapfile > $checkfile") != 0) {
        print "error executing system command for mapcheck on line", $line, "\n";
}  #end if

if(system("./maq cns2win -w 1000 $assemblefile > $winfile") != 0) {
        print "error executing system command for cns2win on line", $line, "\n";
}  #end if

if(system("./maq cns2view $assemblefile > $viewfile") != 0) {
	print "error executing system command for cns2view on line", $line, "\n";
}  #end if

if(system("./maq cns2snp $assemblefile > $snpfile") != 0) {
        print "error executing system command for cns2view on line", $line, "\n";
}  #end if

end;

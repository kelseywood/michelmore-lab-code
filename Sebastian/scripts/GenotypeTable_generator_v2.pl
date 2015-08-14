#!/usr/bin/perl

if(@ARGV < 1) {

	print "Program usage: at least two arguments required, please input them\n";
        print "Firs argument is the file with the reference information.\n";
        print "Later arguments, genotype data to merge on the base file.\n";
        exit 0;
}

$numFiles = scalar @ARGV;

my %posData = ();

for($file = 1; $file < $numFiles; ++$file) {

	open(INPUT, $ARGV[$file]);

       %posData = ();

        while ($line = <INPUT>) {

        	chop($line);

                @posInfo = split("\t", $line);

                $posID = join("-", $posInfo[0], $posInfo[1]);

                $posData{$posID} = [$posInfo[3], $posInfo[4], $posInfo[5], $posInfo[6], $posInfo[7], $posInfo[8]];

        } #end while

        if ($file == 1) {

        	open(BASEFILE, "$ARGV[0]");

        } else {

        	open(BASEFILE, "Genotype_Table.txt");

        } #end elsif

        open(RESULTS, ">Temporary_Genotype_Table.txt");

        if ($file eq 1) {

        	print RESULTS "\t\t\t",$ARGV[0],"\t\t\t\t\t\t" ,$ARGV[$file],"\t\t\n";

        } else {

               $baseLine = <BASEFILE>;

               chop($baseLine);

               print RESULTS $baseLine, "\t", $ARGV[$file],"\t\t\n";

        } #end else

        while ($baseLine = <BASEFILE>) {

	        chop($baseLine);

	        @baseData = split("\t", $baseLine);

                $posID = join ("-", $baseData[0], $baseData[1]);

                if($posData{$posID}[0] ne "") {

	                print RESULTS $baseLine, "\t", $posData{$posID}[0], "\t", $posData{$posID}[1], "\t", $posData{$posID}[2], ":", $posData{$posID}[3], ":", $posData{$posID}[4], ":", $posData{$posID}[5], "\n";

                } else {

                        print RESULTS $baseLine, "\tN\t0\t-:-:-:-\n";

                } #end else

        } #end while

        close(RESULTS);
        close(BASEFILE);

        if(system("mv Temporary_Genotype_Table.txt Genotype_Table.txt") != 0) {
		print "error executing system command for rename output file file ", $ARGV[$file], $line, "\n";
	}  #end if

} #end for

exit;
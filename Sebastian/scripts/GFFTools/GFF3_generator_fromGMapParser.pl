
my $source = $ARGV[1];

open(INPUT, $ARGV[0]);
open(OUTPUT, ">$ARGV[0].gff3");


while($line = <INPUT>) {

	chomp($line);

        @info = split("\t", $line);

        $info[4] =~ s/,//g;

        @range = split("--", $info[4]);

        $scaffold = $info[3];

        $scaffold =~ s/ //;

        print OUTPUT $scaffold, "\t", $source, "\tmRNA\t", $range[0], "\t", $range[1], "\t", $info[7], "\t", $info[2], "\t.\tName=",$info[0], ";ID=",$info[0], ";Coverage=", $info[6], ";Exons=",$info[5], "\n";

        @exons = split("...:  ", $info[8]);

        my $prevIntron = 0;

        if($info[2] eq "+") {

	        my $prevEndPos = $range[0];

	        foreach $exon (@exons) {

	                $modified = $exon;
	                $modified =~ s/\(/_/;
	                $modified =~ s/\)/_/;

	                @split1 = split("_",$modified);

	                @range = split("-", $split1[1]);

	                $exonStartPos = $prevEndPos + $prevIntron;

	                $exonEndPos = $exonStartPos + (abs($range[0] - $range[1]));

	                $prevEndPos = $exonEndPos;

	                @intron = split("->", $exon);
                        if($intron[1] eq "") {
		                @intron = split("<-", $exon);
                        }

	                $prevIntron = substr($intron[1], 6, 9);

	                ++$prevIntron;

	                print OUTPUT $scaffold, "\t", $source, "\tPredictedExon\t", $exonStartPos, "\t", $exonEndPos, "\t.\t", $info[2], "\t.\tParent=",$info[0], ";PositionEST=", $split1[1], "\n";

	        } #end foreach

        } else {

	        my $prevStartPos = $range[1];

	        foreach $exon (@exons) {

	                $modified = $exon;
	                $modified =~ s/\(/_/;
	                $modified =~ s/\)/_/;

	                @split1 = split("_",$modified);

	                @range = split("-", $split1[1]);

	                $exonEndPos = $prevStartPos - $prevIntron;

	                $exonStartPos = $exonEndPos - (abs($range[0] - $range[1]));

	                $prevStartPos = $exonStartPos;

	                @intron = split("->", $exon);
                        if($intron[1] eq "") {
		                @intron = split("<-", $exon);
                        }

	                $prevIntron = substr($intron[1], 6, 9);

	                ++$prevIntron;

	                print OUTPUT $scaffold, "\t", $source, "\tPredictedExon\t", $exonStartPos, "\t", $exonEndPos, "\t.\t", $info[2], "\t.\tParent=",$info[0], ";PositionEST=", $split1[1], "\n";

	        } #end foreach

	} #end else
} #end while

close(INPUT);
close(OUTPUT);

exit;





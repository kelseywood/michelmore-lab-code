#!/usr/bin/perl
use strict;
use warnings;

if(!defined($ARGV[0])) {
	print "Please provide an tRNAscan .ss file\n";
	die "Missing arguments\n";
} #end if not defined argument

## Script design to convert an tRNAscan ss file into a gff and a fasta file

my $ssFile = $ARGV[0];
my $base = $ARGV[0];
$base =~ s{\.[^.]+$}{};

my %SeqStr;

open(SS, $ssFile);

my $header;
my $annotation;
my $intron;
my $seq;
my $str;

open(GFF, ">$base.gff3");
open(FASTA, ">$base.fasta");

while(my $line = <SS>) {

	chomp($line);

	if ($line eq "") { next
			  
	} elsif ($line =~ m/Length/) { $header = $line
		
	} elsif ($line =~ /Type:/) { $annotation = $line
		
	} elsif ($line =~ /Seq:/) { $seq = $line
				  
	} elsif ($line =~ /intron/) { $intron = $line
				  
	} elsif ($line =~ /Str:/) {
		
		$header =~ s/\t*.$//;
		$header =~ s/\(|\)//g;
		my @headerInfo = split(" ", $header);
		my @coordinates = split("-", $headerInfo[1]);
		
		my $strand;
		my $start;
		my $end;
		if ( ($coordinates[0]-$coordinates[1]) < 0) { $strand = "+" ; $start =  $coordinates[0]; $end = $coordinates[1] }
		else { $strand = "-"; $start =  $coordinates[1]; $end = $coordinates[0] }		
		
		my $reference = $headerInfo[0];
		$reference =~ s/\.trna\d*//;
		
		$annotation =~ s/(Type: )|(Anticodon: )|(Score: )//g;
		$annotation =~ s/ /\t/g;
		
		my @annotationInfo = split("\t", $annotation);
		
		my $score = pop(@annotationInfo);
			
		my @seqInfo = split(" ", $seq);
		my @strInfo = split(" ", $line);
		
		print GFF $reference, "\tRNAScan-SE\ttRNA\t", $start, "\t", $end,"\t", $score, "\t", $strand,"\t.\t";
		print GFF "ID=",$headerInfo[0], ";Name=", $headerInfo[0], ";Type=", $annotationInfo[0], ";Anti-codon=", $annotationInfo[1], ";Structure=", $strInfo[1], ";";

		if (defined($intron)) {
			
			$intron =~ s/\(|\)//g;
			my @intronInfo = split(" ", $intron);
			print GFF "PossibleIntron=", $intronInfo[3], ";";
			
		} #end if there's an intron
		
		print GFF "\n";
		
		print FASTA ">", $header, "\n", $seqInfo[1], "\n";
		
		$header = "";
		$annotation = "";
		undef($intron);
		$seq = "";
		$str = "";

	} #end if is seq

} #end reading ss file

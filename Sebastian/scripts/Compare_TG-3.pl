#!usr/bin/perl
#Kevin Stoffel 1/28/10

#This file is to compare two TG files

use strict;
use warnings;

die "Error! TG_Compare.pl <TG_FILE1> <TG_FILE2>" unless @ARGV ==2;

open(TG1, "<$ARGV[0]") or die "Cannot open $ARGV[0]!"; #OPEN Files#
open(TG2, "<$ARGV[1]") or die "Cannot open $ARGV[1]!";
open(SAME, ">Common_SPPs_$ARGV[0]_vs_$ARGV[1]") or die "Cannot Create Common_SPPs_$ARGV[0]_vs_$ARGV[1].txt";
open(UNIQUE1, ">Unique_SPPs_$ARGV[0]") or die "Cannot Create Unique_SPPs_$ARGV[0].txt";
open(UNIQUE2, ">Unique_SPPs_$ARGV[1]") or die "Cannot Create Unique_SPPs_$ARGV[1].txt";
open(STATS, ">Stats_for_$ARGV[0]_vs_$ARGV[1]") or die "Cannot Create Stats_for_$ARGV[0]_vs_$ARGV[1].txt";

my @TG1_array;
my @TG2_array;
my $TG1_lines = 0;
my $TG2_lines = 0;
my (@No_match1, @No_match2);
my $done;
my ($counter1, $counter2);
my $counter3 = 0;
my $counter4 = 0;
my $counter5 = 0;
my $counter6 = 0;

	while(<TG1>){				#Get Data out of files and put into Arrays TG1/2_array#
		chomp;
		my @tmp = split "\t", $_;
		push @TG1_array, [@tmp];
		$TG1_lines++;
	}
print "Done getting data from $ARGV[0].\n";
	while(<TG2>){
		chomp;
		my @tmp = split "\t", $_;
		push @TG2_array, [@tmp];
		$TG2_lines++;
	}
print "Done getting data from $ARGV[1].\n";

close(TG1);					#Close Input Files#
close(TG2);

for (my $i = 0; $i < $TG1_lines; $i++){			#Compare TG1 to TG2 and report common SPPs to Common file and 
	$done = 1;					#Unique TG1 SPPs to Unique_TG1 file.
	$counter1 = 0;
	for (my $j = 0; $j < $TG2_lines and $done != 0; $j++){
		$counter1++;
		if ($TG1_array[$i][0] eq $TG2_array[$j][0]){
			if	(($TG1_array[$i][1] >= $TG2_array[$j][1] and $TG1_array[$i][1] <= $TG2_array[$j][2]+12) or
				($TG1_array[$i][2] <= $TG2_array[$j][2] and $TG1_array[$i][2] >= $TG2_array[$j][1]-12)){
					print SAME "$TG1_array[$i][0]\t$TG1_array[$i][1]\t$TG1_array[$i][2]\t$TG2_array[$j][1]\t$TG2_array[$j][2]\n";
					$done = 0;
			}
			else	{@No_match1 = ($TG1_array[$i][0], $TG1_array[$i][1], $TG1_array[$i][2]);
			}
		}
		else 	{@No_match1 = ($TG1_array[$i][0], $TG1_array[$i][1], $TG1_array[$i][2]);
		}
	}
	if ($counter1 == $TG2_lines){
		print UNIQUE1 "$No_match1[0]\t$No_match1[1]\t$No_match1[2]\t\n";  #Make counter to only print once.
		$counter5++;
	}
}

print "Done with $ARGV[0] by $ARGV[1]\n";

for (my $i = 0; $i < $TG2_lines; $i++){                 #Compare TG1 to TG2 and report common SPPs to Common file and
	$done = 1;                                      #Unique TG1 SPPs to Unique_TG1 file.
	$counter2 = 0;
	for (my $j = 0; $j < $TG1_lines and $done != 0; $j++){
		$counter2++;
		if ($TG2_array[$i][0] eq $TG1_array[$j][0]){
			if	(($TG2_array[$i][1] >= $TG1_array[$j][1] and $TG2_array[$i][1] <= $TG1_array[$j][2]+12) or
				($TG2_array[$i][2] <= $TG1_array[$j][2] and $TG2_array[$i][2] >= $TG1_array[$j][1]-12)){
#					print SAME "$TG1_array[$i][0]\t$TG1_array[$i][1]\t$TG1_array[$i][2]\t$TG2_array[$j][1]\t$TG2_array[$j][2]\n";
					$done = 0;
			}
			else    {@No_match2 = ($TG2_array[$i][0], $TG2_array[$i][1], $TG2_array[$i][2]);
			}
		}
		else    {@No_match2 = ($TG2_array[$i][0], $TG2_array[$i][1], $TG2_array[$i][2]);
		}
	}
	if ($counter2 == $TG1_lines){
		print UNIQUE2 "$No_match2[0]\t$No_match2[1]\t$No_match2[2]\t\n";  #Make counter to only print once.
		$counter6++;
	}
}

print "Done with $ARGV[1] by $ARGV[0]\n";
print "Getting Stats.\n";				#Getting stats


for (my $i = 0; $i < $TG1_lines; $i++){
	if 	($TG1_array[$i][0] ne $TG1_array[$i+1][0]){
		$counter3++;
	}
}
for (my $i = 0; $i < $TG2_lines; $i++){
	if	($TG2_array[$i][0] ne $TG2_array[$i+1][0]){
		$counter4++;
	}
}

print STATS "$ARGV[0] vs $ARGV[1]\n";
print STATS "Number of SPPs for $ARGV[0]: $TG1_lines\n";
print STATS "Number of SPPs for $ARGV[1]: $TG2_lines\n";
print STATS "Number of Contigs for $ARGV[0]: $counter3\n";
print STATS "Number of Contigs for $ARGV[1]: $counter4\n";
print STATS "Number of Unique SPPs for $ARGV[0]: $counter5\n";
print STATS "Number of Unique SPPs for $ARGV[1]: $counter6\n";



close (UNIQUE1);
close (UNIQUE2);


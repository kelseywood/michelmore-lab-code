#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

if(!defined($ARGV[3])) {
	print "Missing inputs\n";
	print "Usage: <TextFileWithLTR> <3LTR.fasta> <5LTR.fasta> <outDir>\n";
	die
} #end else

my $cluster = $ARGV[0];
my $ltr3 = $ARGV[1];
my $ltr5 = $ARGV[2];
my $outDir = $ARGV[3];

my $curDir = cwd();

my $clusterName = $cluster;
$clusterName =~ s{.*/}{};
$clusterName =~ s{\.[^.]+$}{};


if( !(-d $outDir) ) { system("mkdir $outDir") }
if( !(-d "$outDir/$clusterName") ) { system("mkdir $outDir/$clusterName") }
else { print "Analysis have already been done for $clusterName\n Skipping ...\n\n"; exit}


#########################
#
# Getting LTR sequences
#
#########################

print "Extracting LTR sequences out of fasta files\n";

# Extracting 3ltr
if( system("/home/sreyesch/scripts/FastaTools/Fasta_retriever.pl $ltr3 $cluster $outDir/$clusterName/$clusterName.3ltr") ) {
	die "FATAL ERROR: Couldn't extract 3ltr for $$clusterName\n";
} #end

# Extracting 5ltr
if( system("/home/sreyesch/scripts/FastaTools/Fasta_retriever.pl $ltr5 $cluster $outDir/$clusterName/$clusterName.5ltr") ) {
	die "FATAL ERROR: Couldn't extract 5ltr for $$clusterName\n";
} #end

# Adding prefix to fasta files
system("perl -p -i -e 's/>/>LTR3_/g' $outDir/$clusterName/$clusterName.3ltr.fasta\n");
system("perl -p -i -e 's/>/>LTR5_/g' $outDir/$clusterName/$clusterName.5ltr.fasta\n");

# Concatenating the files
system("cat $outDir/$clusterName/$clusterName.3ltr.fasta $outDir/$clusterName/$clusterName.5ltr.fasta > $outDir/$clusterName/$clusterName.ltr.fasta");

print "Sucessfully generated file with ltr sequences\n\n";


###################################
#
# Performing the clustar aligment
#
###################################

print "Performnig Clustal aligment of ltr's\n";

# Clustal aligment
if( system("/home/sreyesch/MichelmoreBin/clustalo-1.2.0-Ubuntu-x86_64 -i $outDir/$clusterName/$clusterName.ltr.fasta -o $outDir/$clusterName/$clusterName.ltr.clustal --force > $outDir/$clusterName/$clusterName.ltr.clustal.log") ) {
	die "FATAL ERROR: Couldn't perform aligment of ltr sequences\n";
} #end

open(CLUSTAL, "$outDir/$clusterName/$clusterName.ltr.clustal");

my %seqs;
my $seq;

while(my $aligmentLine = <CLUSTAL>) {

	chomp($aligmentLine);

	if($aligmentLine =~ /^>/) { $seq = $aligmentLine; $seq =~ s/>// }

	else { $seqs{$seq} .= $aligmentLine}

} #end reading aligment

my @sequences = sort keys(%seqs);

my $numSeq = scalar(@sequences);
my $seqSize = length($seqs{$sequences[0]});

open(PHYLIP, ">$outDir/$clusterName/$clusterName.ltr.phylip");
print PHYLIP " ", $numSeq, " ", $seqSize, "\n";

foreach my $sequence (@sequences) { print PHYLIP $sequence, "  ", $seqs{$sequence}, "\n" }

undef(%seqs);
undef($seq);
undef(@sequences);
undef($numSeq);
undef($seqSize);

print "Aligment and generation of PHYLIP file completed\n\n";


###################################
#
# Generate the tree with RAxML
#
###################################

print "Will generate tree with RAxML\n";

system("/home/sreyesch/MichelmoreBin/standard-RAxML-master/raxmlHPC-PTHREADS -T 2 -s $outDir/$clusterName/$clusterName.ltr.phylip -m GTRCAT -w $curDir/$outDir/$clusterName/ -n $clusterName.RAxML -p 1 > $outDir/$clusterName/RAxML.$clusterName.log");

print "Done generating tree\n\n";


###################################
#
# Calcualte divergence with baseML
#
###################################

print "Time to calculate the sequence divergence with baseML\n";

system("cp /home/sreyesch/MichelmoreBin/paml4.7a/baseml.ctl $outDir/$clusterName/baseml.ctl");

system("perl -p -i -e 's/brown.nuc/$clusterName.ltr.phylip/' $outDir/$clusterName/baseml.ctl");
system("perl -p -i -e 's/brown.trees/RAxML_bestTree.$clusterName.RAxML/' $outDir/$clusterName/baseml.ctl");
system("perl -p -i -e 's/cleandata = 1/cleandata = 0/' $outDir/$clusterName/baseml.ctl");

chdir("$outDir/$clusterName");

if(system("/home/sreyesch/MichelmoreBin/paml4.7a/bin/baseml baseml.ctl > baseml.log") ) {

	die "Couldn't run baseml for $clusterName\n";

} #end if

print "Divergence calculation completed\n\n";

###################################
#
# Estimation of LTR age
#
###################################

print "Finally, time to calculate the ages\n";

open(BASET, "2base.t");

open(AGES, ">$clusterName.LTR_Ages.txt");

$numSeq = <BASET>;

$numSeq =~ s/ //g;

my $seqCounter = 1;
my $LTRindex = 1;

while( my $sequence = <BASET>) {

	if($seqCounter <= $numSeq/2) {
		++$seqCounter;
	} else {
		chomp($sequence);

		$sequence =~ s/  / /g;
		$sequence =~ s/  / /g;
		$sequence =~ s/  / /g;

		my @divergence = split(" ", $sequence);

		print AGES $clusterName, "\t", $divergence[0], "\t", $divergence[$LTRindex], "\t", $divergence[$LTRindex]/(2*(10**-8)), "\n";

		++$seqCounter;
		++$LTRindex

	} #end calculating age

} #end

close(BASET);
close(AGES);

print "Age calculation completed\n\n";

if(-f "../$clusterName.LTR_Ages.txt") {system("rm ../$clusterName.LTR_Ages.txt") }

system("ln -s $clusterName/$clusterName.LTR_Ages.txt ../");

print "ALL DONE!!!!\n\n";

exit;

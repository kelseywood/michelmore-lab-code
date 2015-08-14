#!/usr/bin/perl

######################################################
#                                                    #
#                Exonerate parser                    #
#                                                    #
#                   Developed by                     #
#              Sebastian Reyes-Chin-Wo               #
#         University of California at Davis          #
#                   Genome Center                    #
#               sreyesch@ucdavis.edu                 #
#                                                    #
######################################################

# This script was designed to pares the output from an exonerate analysis with the following options
# --showalignment no --showtargetgff yes --ryo "# --- START OF GENOMIC SEQUENCE DUMP ---\n>%ti  Genomic (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tas\n# --- STOP OF GENOMIC SEQUENCE DUMP ---\n# --- START OF CODING SEQUENCE DUMP ---\n>%ti  Coding (%tab - %tae) strand= %g rank= %r ID= %pi qID= %qi\n%tcs\n# --- STOP OF CONDING SEQUENCE DUMP ---\n"
# It will parse the file and create three fasta files and a gff file

use strict;
use warnings;
use Bio::Seq;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

my $CodonTable = Bio::Tools::CodonTable->new();

if(!defined($ARGV[0])) {
	print "Please provide an exonerate output\n";
	die "Missing file to parse\n";
} #end if not file provided

open(INPUT, $ARGV[0]);

open(GENOMIC, ">$ARGV[0].filtered.genomic.fasta");
open(CODING, ">$ARGV[0].filtered.coding.fasta");
open(PROTEIN, ">$ARGV[0].filtered.protein.fasta");
open(GFF, ">$ARGV[0].filtered.gff");

open(STATS, ">$ARGV[0].filtered.stats.txt");

print STATS "#Target	Query	Start	End	Strand	Rank	Identity	LengthAlignment\n";

while(my $startLine = <INPUT>) {

	if($startLine =~ /^Command line:/) {

		if($startLine =~ /(Bremia_UCOS\.protein\.fasta\.46.fasta)|(Bremia_UCOS\.protein\.fasta\.53.fasta)|(Bremia_UCOS\.protein\.fasta\.55.fasta)|(Bremia_UCOS\.protein\.fasta\.68.fasta)|(Bremia_UCOS\.protein\.fasta\.98.fasta)|(Bremia_UCOS\.protein\.fasta\.133.fasta)|(Bremia_UCOS\.protein\.fasta\.153.fasta)|(Bremia_UCOS\.protein\.fasta\.154.fasta)|(Bremia_UCOS\.protein\.fasta\.160.fasta)|(Bremia_UCOS\.protein\.fasta\.162.fasta)|(Bremia_UCOS\.protein\.fasta\.204.fasta)|(Bremia_UCOS\.protein\.fasta\.212.fasta)|(Bremia_UCOS\.protein\.fasta\.216.fasta)|(Bremia_UCOS\.protein\.fasta\.225.fasta)|(Bremia_UCOS\.protein\.fasta\.237.fasta)|(Bremia_UCOS\.protein\.fasta\.241.fasta)|(Bremia_UCOS\.protein\.fasta\.246.fasta)/) {next}

		my $name;

if($startLine =~ /Bremia_UCOS\.protein\.fasta\.1\./) {$name = "ctg7180000009127_3566-4192" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.2\./) {$name = "ctg7180000010514_9671-10537" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.3\./) {$name = "ctg7180000010344_14279-10214" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.4\./) {$name = "ctg7180000010710_24906-26184" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.5\./) {$name = "ctg7180000015173_3030-3732" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.6\./) {$name = "ctg7180000015173_2870-1784" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.7\./) {$name = "ctg7180000008809_40209-39278" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.8\./) {$name = "ctg7180000011488_15940-14356" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.9\./) {$name = "ctg7180000009101_30930-31675" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.10\./) {$name = "ctg7180000014182_974-6992" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.11\./) {$name = "ctg7180000008743_19282-18231" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.12\./) {$name = "ctg7180000010310_21928-21220" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.13\./) {$name = "ctg7180000008208_46017-44848" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.14\./) {$name = "ctg7180000009063_21041-20291" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.15\./) {$name = "ctg7180000008583_39610-38208" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.16\./) {$name = "ctg7180000009316_13824-12180" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.17\./) {$name = "ctg7180000010126_4876-5701" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.18\./) {$name = "ctg7180000009609_7721-8914" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.19\./) {$name = "ctg7180000008414_16824-17818" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.20\./) {$name = "ctg7180000010166_31228-12643" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.21\./) {$name = "ctg7180000014203_2368-4855" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.22\./) {$name = "ctg7180000009892_1583-3026" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.23\./) {$name = "ctg7180000008181_62115-63352" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.24\./) {$name = "ctg7180000013762_2266-1573" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.25\./) {$name = "ctg7180000011081_1452-972" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.26\./) {$name = "ctg7180000009466_37831-36207" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.27\./) {$name = "ctg7180000009394_8928-7414" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.28\./) {$name = "ctg7180000008298_58781-59804" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.29\./) {$name = "ctg7180000008298_53085-54531" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.30\./) {$name = "ctg7180000012464_11141-9833" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.31\./) {$name = "ctg7180000009617_6852-5130" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.32\./) {$name = "ctg7180000013480_2908-3502" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.33\./) {$name = "ctg7180000010668_317-2156" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.34\./) {$name = "ctg7180000011198_10607-10960" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.35\./) {$name = "ctg7180000010964_5584-5044" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.36\./) {$name = "ctg7180000010964_8773-7646" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.37\./) {$name = "ctg7180000011482_8300-13477" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.38\./) {$name = "ctg7180000011386_2123-3328" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.39\./) {$name = "ctg7180000012162_14538-15525" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.40\./) {$name = "ctg7180000014184_11877-12967" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.41\./) {$name = "ctg7180000008555_17233-16267" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.42\./) {$name = "ctg7180000015156_1514-674" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.43\./) {$name = "ctg7180000012577_5350-8573" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.44\./) {$name = "ctg7180000008204_15156-16816" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.45\./) {$name = "ctg7180000009726_9863-10818" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.47\./) {$name = "ctg7180000009150_27765-25304" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.48\./) {$name = "ctg7180000008271_7820-6967" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.49\./) {$name = "ctg7180000013654_3284-4037" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.50\./) {$name = "ctg7180000015228_1016-50" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.51\./) {$name = "ctg7180000008649_1893-530" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.52\./) {$name = "ctg7180000008172_95559-96453" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.54\./) {$name = "ctg7180000013874_6436-5734" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.56\./) {$name = "ctg7180000008208_58010-54257" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.57\./) {$name = "ctg7180000011523_18510-19395" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.58\./) {$name = "ctg7180000009482_3530-2727" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.59\./) {$name = "ctg7180000010803_609-1542" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.60\./) {$name = "ctg7180000008460_49425-50427" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.61\./) {$name = "ctg7180000008934_54110-53420" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.62\./) {$name = "ctg7180000008934_18980-20993" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.63\./) {$name = "ctg7180000008949_36013-33442" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.64\./) {$name = "ctg7180000010596_2645-1475" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.65\./) {$name = "ctg7180000008454_5742-8019" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.66\./) {$name = "ctg7180000008454_11611-9694" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.67\./) {$name = "ctg7180000009221_7293-8343" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.69\./) {$name = "ctg7180000009963_9456-11325" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.70\./) {$name = "ctg7180000008690_776-3356" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.71\./) {$name = "ctg7180000008777_32039-33340" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.72\./) {$name = "ctg7180000010960_253-1846" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.73\./) {$name = "ctg7180000014992_532-2902" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.74\./) {$name = "ctg7180000013218_1901-4790" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.75\./) {$name = "ctg7180000009276_26045-24674" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.76\./) {$name = "ctg7180000009276_23004-24135" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.77\./) {$name = "ctg7180000013661_5447-4973" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.78\./) {$name = "ctg7180000015365_5599-9037" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.79\./) {$name = "ctg7180000010007_12311-11441" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.80\./) {$name = "ctg7180000008573_7606-11701" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.81\./) {$name = "ctg7180000012701_3849-2943" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.82\./) {$name = "ctg7180000008437_12178-13619" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.83\./) {$name = "ctg7180000011761_8379-10197" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.84\./) {$name = "ctg7180000013845_5050-4258" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.85\./) {$name = "ctg7180000011538_9710-11917" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.86\./) {$name = "ctg7180000015566_8371-7747" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.87\./) {$name = "ctg7180000012224_9678-11018" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.88\./) {$name = "ctg7180000008454_44816-46307" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.89\./) {$name = "ctg7180000010333_16763-13241" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.90\./) {$name = "ctg7180000012699_3398-2131" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.91\./) {$name = "ctg7180000009262_33709-34809" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.92\./) {$name = "ctg7180000015266_1021-4410" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.93\./) {$name = "ctg7180000012438_5421-3566" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.94\./) {$name = "ctg7180000012438_7758-5928" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.95\./) {$name = "ctg7180000012041_9566-8987" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.96\./) {$name = "ctg7180000014235_2925-3678" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.97\./) {$name = "ctg7180000011507_7108-6247" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.99\./) {$name = "ctg7180000008257_14704-11528" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.100\./) {$name = "ctg7180000011528_8851-10108" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.101\./) {$name = "ctg7180000011528_8738-7424" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.102\./) {$name = "ctg7180000010915_9001-7413" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.103\./) {$name = "ctg7180000009103_2138-562" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.104\./) {$name = "ctg7180000012530_4441-1354" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.105\./) {$name = "ctg7180000013049_13211-13606" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.106\./) {$name = "ctg7180000012592_3143-3629" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.107\./) {$name = "ctg7180000010932_9005-10293" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.108\./) {$name = "ctg7180000010991_1676-3768" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.109\./) {$name = "ctg7180000014131_19821-20571" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.110\./) {$name = "ctg7180000009677_3707-4450" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.111\./) {$name = "ctg7180000008280_21716-18809" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.112\./) {$name = "ctg7180000009497_9052-9895" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.113\./) {$name = "ctg7180000009497_13184-12494" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.114\./) {$name = "ctg7180000010501_9873-13197" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.115\./) {$name = "ctg7180000008464_1616-2447" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.116\./) {$name = "ctg7180000010146_31919-30553" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.117\./) {$name = "ctg7180000009906_14796-18183" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.118\./) {$name = "ctg7180000009130_16424-15016" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.119\./) {$name = "ctg7180000008957_32344-31676" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.120\./) {$name = "ctg7180000014211_8315-10157" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.121\./) {$name = "ctg7180000008871_5829-6872" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.122\./) {$name = "ctg7180000010149_22044-17745" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.123\./) {$name = "ctg7180000008146_6479-5621" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.124\./) {$name = "ctg7180000011087_3937-9218" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.125\./) {$name = "ctg7180000012661_632-57" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.126\./) {$name = "ctg7180000010736_16401-14757" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.127\./) {$name = "ctg7180000009473_22861-21905" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.128\./) {$name = "ctg7180000009954_7917-9184" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.129\./) {$name = "ctg7180000010961_20500-19550" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.130\./) {$name = "ctg7180000010823_7630-9576" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.131\./) {$name = "ctg7180000015757_3457-1494" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.132\./) {$name = "ctg7180000010904_16038-15575" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.134\./) {$name = "ctg7180000011896_7875-5328" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.135\./) {$name = "ctg7180000012302_8817-12268" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.136\./) {$name = "ctg7180000009623_17572-19132" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.137\./) {$name = "ctg7180000012162_12145-14322" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.138\./) {$name = "ctg7180000010057_19544-20462" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.139\./) {$name = "ctg7180000008516_25540-23371" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.140\./) {$name = "ctg7180000008516_25792-29636" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.141\./) {$name = "ctg7180000008488_40778-38494" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.142\./) {$name = "ctg7180000009006_34718-32972" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.143\./) {$name = "ctg7180000011763_7046-8330" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.144\./) {$name = "ctg7180000008669_24074-26687" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.145\./) {$name = "ctg7180000008227_39487-42529" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.146\./) {$name = "ctg7180000009529_22274-21713" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.147\./) {$name = "ctg7180000009962_8464-7037" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.148\./) {$name = "ctg7180000008431_3119-1721" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.149\./) {$name = "ctg7180000015577_10616-9204" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.150\./) {$name = "ctg7180000013935_8846-7398" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.151\./) {$name = "ctg7180000009992_14931-12903" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.152\./) {$name = "ctg7180000010274_13757-14444" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.155\./) {$name = "ctg7180000010768_13280-14277" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.156\./) {$name = "ctg7180000008328_23292-22656" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.157\./) {$name = "ctg7180000012311_4780-2901" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.158\./) {$name = "ctg7180000008939_5691-11244" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.159\./) {$name = "ctg7180000008504_14014-12596" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.161\./) {$name = "ctg7180000008829_29458-27920" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.163\./) {$name = "ctg7180000009054_12746-11006" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.164\./) {$name = "ctg7180000010254_1751-2999" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.165\./) {$name = "ctg7180000010941_8276-6819" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.166\./) {$name = "ctg7180000009479_18027-20553" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.167\./) {$name = "ctg7180000008399_12153-10951" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.168\./) {$name = "ctg7180000008695_24491-22890" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.169\./) {$name = "ctg7180000009414_20315-21485" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.170\./) {$name = "ctg7180000008672_23879-14878" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.171\./) {$name = "ctg7180000009703_2037-2765" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.172\./) {$name = "ctg7180000009263_11168-8438" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.173\./) {$name = "ctg7180000008261_10413-8881" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.174\./) {$name = "ctg7180000009838_18886-17875" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.175\./) {$name = "ctg7180000009610_15745-16525" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.176\./) {$name = "ctg7180000014184_5022-5702" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.177\./) {$name = "ctg7180000009938_18849-18360" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.178\./) {$name = "ctg7180000011488_13749-9895" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.179\./) {$name = "ctg7180000009491_5546-4648" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.180\./) {$name = "ctg7180000011765_1601-4604" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.181\./) {$name = "ctg7180000008643_44984-45683" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.182\./) {$name = "ctg7180000010095_4466-7553" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.183\./) {$name = "ctg7180000008754_13264-14076" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.184\./) {$name = "ctg7180000015404_4625-3362" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.185\./) {$name = "ctg7180000010193_39221-40766" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.186\./) {$name = "ctg7180000008733_20720-18900" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.187\./) {$name = "ctg7180000012486_5886-5016" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.188\./) {$name = "ctg7180000009715_10109-9706" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.189\./) {$name = "ctg7180000008811_33200-32201" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.190\./) {$name = "ctg7180000008811_31208-32068" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.191\./) {$name = "ctg7180000009121_20984-22145" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.192\./) {$name = "ctg7180000008736_32982-34047" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.193\./) {$name = "ctg7180000008736_7771-8880" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.194\./) {$name = "ctg7180000009103_33453-32202" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.195\./) {$name = "ctg7180000008852_44353-45304" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.196\./) {$name = "ctg7180000012412_1920-2308" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.197\./) {$name = "ctg7180000009772_6275-4845" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.198\./) {$name = "ctg7180000008509_4071-5922" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.199\./) {$name = "ctg7180000010309_8463-11952" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.200\./) {$name = "ctg7180000009416_908-2536" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.201\./) {$name = "ctg7180000011540_12239-11114" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.202\./) {$name = "ctg7180000008285_10486-8906" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.203\./) {$name = "ctg7180000008615_5884-6928" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.205\./) {$name = "ctg7180000011333_6192-7320" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.206\./) {$name = "ctg7180000008249_4209-1993" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.207\./) {$name = "ctg7180000008249_4487-10736" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.208\./) {$name = "ctg7180000011165_6936-7693" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.209\./) {$name = "ctg7180000011529_18992-16587" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.210\./) {$name = "ctg7180000011529_3583-2898" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.211\./) {$name = "ctg7180000008450_32483-29733" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.213\./) {$name = "ctg7180000008301_1282-3031" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.214\./) {$name = "ctg7180000009713_2012-3720" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.215\./) {$name = "ctg7180000008275_57715-58334" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.217\./) {$name = "ctg7180000009989_18401-17583" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.218\./) {$name = "ctg7180000012158_12897-11476" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.219\./) {$name = "ctg7180000008454_50426-51119" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.220\./) {$name = "ctg7180000013395_8510-9262" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.221\./) {$name = "ctg7180000015353_2152-1489" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.222\./) {$name = "ctg7180000008415_6220-2983" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.223\./) {$name = "ctg7180000008439_6888-7644" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.224\./) {$name = "ctg7180000008189_29288-29713" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.226\./) {$name = "ctg7180000013773_3523-4922" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.227\./) {$name = "ctg7180000008252_21022-19585" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.228\./) {$name = "ctg7180000012393_1035-2542" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.229\./) {$name = "ctg7180000012525_3220-4193" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.230\./) {$name = "ctg7180000009220_30969-28590" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.231\./) {$name = "ctg7180000008813_32178-35877" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.232\./) {$name = "ctg7180000008510_30036-29484" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.233\./) {$name = "ctg7180000008153_45418-47227" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.234\./) {$name = "ctg7180000008510_28071-26470" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.235\./) {$name = "ctg7180000008167_23768-24668" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.236\./) {$name = "ctg7180000008189_21190-22009" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.238\./) {$name = "ctg7180000010159_10204-11453" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.239\./) {$name = "ctg7180000014167_2102-923" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.240\./) {$name = "ctg7180000008297_73019-72290" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.242\./) {$name = "ctg7180000014221_5296-6392" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.243\./) {$name = "ctg7180000008156_74675-66404" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.244\./) {$name = "ctg7180000012178_6394-1934" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.245\./) {$name = "ctg7180000009239_12879-13692" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.247\./) {$name = "ctg7180000012788_4435-3349" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.248\./) {$name = "ctg7180000008653_12818-14628" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.249\./) {$name = "ctg7180000010662_11958-11329" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.250\./) {$name = "ctg7180000012214_2950-817" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.251\./) {$name = "ctg7180000008325_47967-50058" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.252\./) {$name = "ctg7180000008824_25682-27893" }
elsif($startLine =~ /Bremia_UCOS\.protein\.fasta\.253\./) {$name = "ctg7180000008400_31197-27657" }


		while (my $line = <INPUT>) {

			if ($line =~ /-- completed exonerate analysis/) {last}

			if ($line =~ /--- START OF GFF DUMP ---/) {

				while(my $gffLine = <INPUT>) {

					if ($gffLine =~ /--- END OF GFF DUMP ---/) {last}

					print GFF $gffLine;

				} #end while for gff

			} elsif ($line =~ /--- START OF GENOMIC SEQUENCE DUMP ---/) {

				while(my $genomicLine = <INPUT>) {

					if ($genomicLine =~ /--- STOP OF GENOMIC SEQUENCE DUMP ---/) {last}

					print GENOMIC $genomicLine;

				} #end while for genomic sequence

			} elsif ($line =~ /--- START OF CODING SEQUENCE DUMP ---/) {

				my $codingSequence;

				my $header;

				while(my $codingLine = <INPUT>) {

					if ($codingLine =~ /--- STOP OF CONDING SEQUENCE DUMP ---/) {

						$header =~ s/Coding/Protein/;

						print PROTEIN $header;

						my $codonStart = 0;

						my $aa = "gibberish";

						while ($aa ne "") {

							my $codon = substr($codingSequence, $codonStart, 3);

							$aa = $CodonTable->translate_strict($codon);

							print PROTEIN $aa;

							$codonStart += 3;

						} #end while

						print PROTEIN "\n";

						my @info = split(" ", $header);

						print STATS substr($info[0],1,300),"_",substr($info[2],1,30),"-",substr($info[4],0,-1), "\t", $name, "\t", substr($info[2],1,30), "\t", substr($info[4],0,-1), "\t", $info[6], "\t", $info[8], "\t", $info[10], "\t", length($codingSequence), "\n";

						last

					} elsif ( !($codingLine =~ /^>/) && $line ne "") {chomp($codingLine); $codingSequence .= $codingLine

					} elsif ($codingLine =~ /^>/) {$header = $codingLine} #end for elsif

					print CODING $codingLine;

				} #end while for coding sequence

			} # End of looking for printing blocks

		} #end if reading block

	} #end if command line to print

} #end while

close(INPUT);
close(GENOMIC);
close(CODING);
close(PROTEIN);
close(GFF);

exit;

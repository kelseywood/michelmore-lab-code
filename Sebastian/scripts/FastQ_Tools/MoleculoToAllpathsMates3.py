#!/usr/bin/env python
# for python 2.x

#
# This script generates pseudo-Illumina mate-pair reads from long FASTQ files as for example generated by Moleculo technology.
# The resulting synthetic mate pair reads have the insert size indicate by the variable "mp" (e.g. 1400)  and are oriented like Illumina mate pair data, i. e.  reverse-forward.
# The synthetic reads are 100 bp long are generated every 200 bases for shorter Moleculo reads ( < 2000 bp), every 400 bases for longer Moleculo reads ( < 10kb),
# and every 800 bases for even longer Moleculo reads.
#
# Usage: python MoleculoToAllpathsMates2.py fastqfile.fq   MP-insertsize   minus  &
# e.g.: python MoleculoToAllpathsMates2.py   MoleculoReads.fq   1400   100  &

import os, sys, math

####################################
#reverse complement
def comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'R', 'Y': 'Y', 'S': 'S', 'W': 'W', 'M': 'M', 'K': 'K'}
    complseq = [complement[base] for base in seq]
    return complseq

def revcomp(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(comp(seq))

####################################

mp = int(sys.argv[2])  ### sets Mate Pair distance
minus = int(sys.argv[3])  ### basepairs removed from ends of Moleculo reads before processing

u = open(sys.argv[1])
sname = (sys.argv[1]).replace(".fq", "")
o1 = open(sname +'APMateR-' + sys.argv[2] + "minus" + sys.argv[3] + "-1.fq" ,'w')
o2 = open(sname +'APMateR-' + sys.argv[2] + "minus" + sys.argv[3] + "-2.fq" ,'w')

while True:
	name = u.readline()
	if name == '':
		break
	if name == '\n':
		break
	seq = u.readline()
	seq = seq[minus:-minus]
	fc = 0
	l = len(seq)
#	print l
	x = (l - mp - 1)
	q = x/200
	if minus == 0:
		seq = seq[:-1] # + "\n"
	plus = u.readline()
	plus = "+\n"
	qual = u.readline()
	qual = qual[minus:-minus]
	if minus == 0:
		qual = qual[:-1] # + "\n"

	if q < 10:
		for i in range(1, q+1):
			nn1 = name[:-1] + "_f1__" + str(i).zfill(4) + "\n"
			nn2 = name[:-1] + "_r1__" + str(i).zfill(4) + "\n"

			ns1 = seq[((i * 200) - 200):((i * 200) - 100)] +"\n"
			ns2 = seq[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"
			nq1 = qual[(i * 200 - 200):(i * 200 - 100)] +"\n"
			nq2 = qual[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"

			ns1 = ns1[:-1]
			nq1 = nq1[:-1]
			ns1 = revcomp(ns1) + "\n"
			nq1 = nq1[::-1] + "\n"

			o1.write(nn1 + ns1 + plus + nq1)
			o2.write(nn2 + ns2 + plus + nq2)

	if 10 < q < 5000 :
		for i in range(1, q+1, 2):
			nn1 = name[:-1] + "_f1__" + str(i).zfill(4) + "\n"
			nn2 = name[:-1] + "_r1__" + str(i).zfill(4) + "\n"

			ns1 = seq[((i * 200) - 200):((i * 200) - 100)] +"\n"
			ns2 = seq[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"
			nq1 = qual[(i * 200 - 200):(i * 200 - 100)] +"\n"
			nq2 = qual[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"

			ns1 = ns1[:-1]
			nq1 = nq1[:-1]
			ns1 = revcomp(ns1) + "\n"
			nq1 = nq1[::-1] + "\n"

			o1.write(nn1 + ns1 + plus + nq1)
			o2.write(nn2 + ns2 + plus + nq2)

	if q >= 5000 :
		for i in range(1, q+1, 4):
			nn1 = name[:-1] + "_f1__" + str(i).zfill(4) + "\n"
			nn2 = name[:-1] + "_r1__" + str(i).zfill(4) + "\n"

			ns1 = seq[((i * 200) - 200):((i * 200) - 100)] +"\n"
			ns2 = seq[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"
			nq1 = qual[(i * 200 - 200):(i * 200 - 100)] +"\n"
			nq2 = qual[(i * 200 - 100 + mp):(i * 200 + mp)]  +"\n"

			ns1 = ns1[:-1]
			nq1 = nq1[:-1]
			ns1 = revcomp(ns1) + "\n"
			nq1 = nq1[::-1] + "\n"

			o1.write(nn1 + ns1 + plus + nq1)
			o2.write(nn2 + ns2 + plus + nq2)


u.close()
o1.close()
o2.close()
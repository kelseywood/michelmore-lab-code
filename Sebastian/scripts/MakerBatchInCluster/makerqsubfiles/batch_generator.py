#! /usr/bin/env python
import os, sys


f1 = open(sys.argv[1], "r")
#f2 = open(sys.argv[2], "r")
filename = sys.argv[1]
filename = filename.replace(".txt", "")
while True: 
    qsub = f1.readline()
    if qsub == '':
         break
    qsub = qsub.replace('\n', '')
    f2 = open(sys.argv[2], "r")
    o1 = open("qsub_" + qsub + ".txt", "w")
    while True:
        script = f2.readline()
        if script == '':
            break
        script = script.replace('x-x', qsub)
        o1.write(script)
    f2.close()
    o1.close()

f1.close()

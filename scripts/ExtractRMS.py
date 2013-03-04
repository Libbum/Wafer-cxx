#!/usr/bin/python

import re
import os

fm = re.compile('zrun\d+.log')
cl = re.compile('Parameters from commandline')
tl = re.compile('T = ')
dl = re.compile('==> Ground State Binding Energy : ')

output = []

dirList = os.listdir('.')
for fname in dirList:
    m1 = fm.match(fname)
    if m1:
        file = open(fname,"r")
	fp = False
        for line in file:
            m2 = cl.match(line)
	    if m2:
		fp = True
	    m3 = tl.match(line)
	    if (m3 and fp):
	      (junk, t) = line.split('= ')
	      t = eval(t.strip())
	    m4 = dl.match(line)
            if m4:
                (junk, data) = line.split(": ")
                data = data.strip()
		(re,im) = data.split(',')
		re = -eval(re[1:])
		im = eval(im[:-1])
        output += [[t,re,im]]
        file.close()

soutput = sorted(output)

for item in soutput:
   print "%f\t%f\t%f" % (item[0],item[1],item[2])

#!/usr/bin/python

import re
import os

batch = ["_128","_256",""]

cl = re.compile('Parameters from commandline')
tl = re.compile('T = ')
dl = re.compile('==> Ground State Binding Energy : ')

for id in batch:
    output = []
    fm = re.compile("zrun\d+%s.log" % id)
    dirList = os.listdir('.')
    for fname in dirList:
        fd = False
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
		    (rep,imp) = data.split(',')
		    rep = -eval(rep[1:])
		    imp = eval(imp[:-1])
                    fd = True
            if fd:
                output += [[t,rep,imp]]
            file.close()

    soutput = sorted(output)

    of = open("be%s.txt" % id,"w")
    for item in soutput:
       of.write("%f\t%f\t%f\n" % (item[0],item[1],item[2]))
    of.close()

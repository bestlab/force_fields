#!/usr/bin/env python

import sys
import os

## The script uses parmed.py in AmberTools.
if len(sys.argv)!=3:
	print('Usage: ./ff03w-to-ff03ws.py <old prmtop> <new prmtop>\n')
	quit()


# input prmtop
fnin=sys.argv[1]
# output prmtop
fnout=sys.argv[2]

fscript=open('script_show','w')
fscript.write('printLJMatrix @%OW\n')
fscript.close()
os.system('parmed.py %s script_show>oldLJ.dat'%(fnin))

factor=1.1

fscript=open('script_change','w')
with open('oldLJ.dat','r') as fid:
	status=False
	for i in fid.readlines():
		line=i.rsplit()
		if len(line)!=0:
			if line[0]=="Atom" and line[1]=="Type":
				status=True
			if status:
				if len(line)==8:
					atomtype1=line[0].rsplit(',')[0]
					atomtype2=line[2].rsplit(',')[0]
					r12=float(line[6])
					eps12=float(line[7])
					if atomtype1!=atomtype2 and eps12!=0:
						fscript.write('changeLJPair @%%%s @%%%s %9.6f %9.6f\n'%(atomtype1,atomtype2,r12,eps12*factor))
fscript.write('outparm %s'%(fnout))
fscript.close()

os.system('parmed.py %s script_change'%(fnin))
os.system('parmed.py %s script_show>newLJ.dat'%(fnout))


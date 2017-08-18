#!/usr/bin/env python

import sys,os

nres=40
peptide="abeta%i"%(nres)
######################################################################
ff="amber99sbws_amberdye"
# alexa 488 -- N-terminal C5 oxime
a488pdb="../PDB_structures_patch/Alexa_488/A48_B1N.pdb"
# alexa 647 -- C-ter C2 maleimide
a647pdb="../PDB_structures_patch/Alexa_647/A64_C2C.pdb"
#
fileroot="%s_dyes_%s" % (peptide,ff)
boxlen=7.	# nm
np=20 
nn=12
######################################################################
pdbfile = "%s.pdb"%(fileroot)
grofileH = "data/%s_H.gro"%(fileroot)
pdbfileH = "data/%s_H.pdb"%(fileroot)
pdbtmp1 = "data/%s_tmp1.pdb"%(fileroot)
protdyepdb = "data/%s_protdye.pdb"%(fileroot)
protdyegro = "data/%s_vac.gro"%(fileroot)
protdyetop = "%s_vac.top"%(fileroot)
protdyevacmin = "data/%s_vac_min.gro"%(fileroot)
protdyevacdyn = "data/%s_vac_dyn.gro"%(fileroot)
gro_ed = "data/%s_ed.gro"%(fileroot)
protdyewattop = "%s_wat.top"%(fileroot)
protdyeionstop = "%s_ions.top"%(fileroot)
protdyewatgro = "data/%s_wat.gro"%(fileroot)
protdyewatgromin = "data/%s_wat_min.gro"%(fileroot)
protdyeionsgro = "data/%s_ions.gro"%(fileroot)
protdyeionsgromin = "data/%s_ions_min.gro"%(fileroot)
GMXBIN="~/gromacs465/bin"  # set to your gromacs installation...

#======================================================================
# 1. build peptide with ambertools
seq = "NALA ASP ALA GLU PHE ARG HIS ASP SER GLY TYR \
GLU VAL HIS HIS GLN LYS LEU VAL PHE PHE \
ALA GLU ASP VAL GLY SER ASN LYS GLY ALA \
ILE ILE GLY LEU MET VAL GLY GLY VAL VAL CCYS"

# This just provides a way to build the peptide in extended conformation
# using ambertools
#
#amberhome="/u/best/lib/amber12"
#leap_inp = """
#source leaprc.ff03.r1
#pep = sequence { %s } 
#saveAmberParm pep abeta40.prmtop abeta40.inpcrd
#savePdb pep %s
#quit 
#""" % ( seq, pdbfile )
#inp = open("leap.inp","w")
#inp.write(leap_inp)
#inp.close()
#
#os.system("%s/bin/tleap -s -f leap.inp"%(amberhome))

# ======================================================================
# 2. generate protein gro with gromacs hydrogen names, so we won't have to build H
# (can't directly make pdb because pdb2gmx generates mangled PDB files ;-)
os.system("export GMXLIB=${HOME}/gromacs/share/top")
os.system("""%s/pdb2gmx -ignh -ff %s -f %s \
	-o %s -p TEST.top  << EOF
1
EOF"""%(GMXBIN,ff,pdbfile,grofileH))

## now convert gro to pdb 
os.system("""%s/trjconv -s %s -f %s -o %s << EOF
0
EOF"""%(GMXBIN,grofileH,grofileH,pdbfileH))

## ======================================================================
## now add dyes and setup in gromacs
##
#
## a) 488 nter
os.system("../build_dyes_amberdye.py 1 %s %s %s"%(a488pdb,pdbfileH,pdbtmp1))

# b) 647 cter (+2 because the dyes are attached to residues added to the end
# of the native sequence, not replacing native residues in this case). 
os.system("../build_dyes_amberdye.py %i %s %s %s"%(nres+2,a647pdb,pdbtmp1,protdyepdb))
#
os.system("""%s/pdb2gmx -ff %s -f %s -o %s -p %s << EOF
1
EOF
"""%(GMXBIN,ff,protdyepdb,protdyegro,protdyetop))

# ======================================================================
# must do initial minimization because of the hackish way dye has been transplanted into protein
# (essentially because phi/psi of new residue may differ from the original)
# seems conjugate gradients works best for this.
os.system("%s/grompp  -p %s -f miniconj.mdp -c %s -o minivac.tpr"%(GMXBIN,protdyetop,protdyegro))

os.system("%s/mdrun -pd -s minivac.tpr -c %s -g minivac.log"%(GMXBIN,protdyevacmin))

os.system("%s/editconf -bt cubic -f %s -o tmp.gro -box 15" \
               %(GMXBIN,protdyevacmin))

# ======================================================================
# Some vacuum dynamics to collapse it a bit...
os.system("%s/grompp -f vac_dyn.mdp -c tmp.gro -o vacdyn.tpr -p %s"%(GMXBIN,protdyetop))
os.system("%s/mdrun  -nt 16 -s vacdyn.tpr -c %s -deffnm junk"%(GMXBIN,protdyevacdyn))

# ======================================================================
# setup and add water box
os.system("%s/editconf -f %s -o %s -box %f -bt octahedron "%(GMXBIN,protdyevacdyn,gro_ed,boxlen))

os.system("cp %s %s "%(protdyetop, protdyewattop))

os.system("%s/genbox -cp %s -cs tip4p2005.gro -o  %s -p %s"%(GMXBIN,gro_ed,protdyewatgro,protdyewattop))

os.system("%s/grompp -v -f mini.mdp -c %s -o mini.tpr -p %s"%(GMXBIN,protdyewatgro,protdyewattop))

os.system("%s/mdrun -pd -v -s mini.tpr -c %s"%(GMXBIN,protdyewatgromin))

# now add ions etc...
# ======================================================================
os.system("%s/grompp -v -f mini.mdp -c %s -o genion.tpr -p %s"%(GMXBIN,protdyewatgromin,protdyewattop))
os.system("cp %s %s "%(protdyewattop,protdyeionstop))

os.system("""%s/genion -s genion.tpr -nn %i -np %i -o %s -p %s<<EOF
13
EOF"""%(GMXBIN,nn,np,protdyeionsgro,protdyeionstop))

os.system("%s/grompp -v -f mini.mdp -c %s -o mini.tpr -p %s"%(GMXBIN,protdyeionsgro,protdyeionstop))

os.system("%s/mdrun -pd -v -s mini.tpr -c %s"%(GMXBIN,protdyeionsgromin))



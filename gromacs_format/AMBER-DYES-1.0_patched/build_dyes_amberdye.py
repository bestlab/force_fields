#!/usr/bin/env python

# quick hack to add dyes without external programs

# version for "AMBER-dyes" force field

import sys,os,math,numpy,random

# ============================================================
Usage = """

	build_dyes_amberdye.py residue_number dye.pdb protein.pdb out.pdb

	e.g. ./build_dyes.py 33 amberdye.pdb myprotein.pdb
	will take the template structure amberdye.pdb (from amber-dyes force
        field) and attempt to use it to replace residue 33 of the supplied 
        myprotein.pdb. Output to out.pdb

"""

# ============================================================

if len(sys.argv) != 5:
	print Usage

resnum = int(sys.argv[1])
dyepdb = sys.argv[2]
pdbinp = sys.argv[3]
pdbout = sys.argv[4]

sys.stdout.write("""

Adding chromophore structure from file:
    "%s"

to residue %i in pdb:
    "%s"

"""%(dyepdb,resnum,pdbinp))
sys.stdout.write("""...and writing result to pdb:
    "%s"

"""%(pdbout))

# ============================================================

def residue(line):
	return int(line[23:26])

def makevec(line):
	return numpy.array(line[30:54].split(),numpy.float64)

def dist(x1,x2):
    return numpy.linalg.norm(x1-x2)

def make_rotmat(CA_xyz,C_xyz,N_xyz):
    # X = CA-C
    # Z = X x (CA-N)
    # Y = X x Z
    X = C_xyz -  CA_xyz
    X  /= numpy.linalg.norm(X)
    CAN = N_xyz -  CA_xyz
    Z = numpy.cross(X,CAN)
    Z  /= numpy.linalg.norm(Z)
    Y = numpy.cross(Z,X)
    M = numpy.array( [X,Y,Z], numpy.float64 )
    return M

def calc_rotmat(CA_old_xyz, C_old_xyz, N_old_xyz, CA_new_xyz, C_new_xyz, N_new_xyz):

    old_X = C_old_xyz -  CA_old_xyz
    old_CAN = N_old_xyz -  CA_old_xyz
    old_X  /= numpy.linalg.norm(old_X)
    old_Z = numpy.cross(x,xnew)
    old_Z  /= numpy.linalg.norm(old_Z)
    old_Y = numpy.cross(old_Z,old_X)

# ============================================================

dye_pdblines = filter(lambda x: x[0:6] == "ATOM  ", open(dyepdb).readlines())
prot_pdblines = filter(lambda x: x[0:6] == "ATOM  ", open(pdbinp).readlines())

firstres = residue(prot_pdblines[0])
lastres = residue(prot_pdblines[-1])

pre = filter(lambda x: residue(x)<resnum, prot_pdblines)
oldres = filter(lambda x: residue(x) == resnum, prot_pdblines)
post = filter(lambda x: residue(x)>resnum, prot_pdblines)

CA_old = filter(lambda x: x[13:16] == "CA ", oldres)[0]
C_old = filter(lambda x: x[13:16] == "C  ", oldres)[0]
N_old = filter(lambda x: x[13:16] == "N  ", oldres)[0]

# This one is easy, it seems to always be present...
CA_new = filter(lambda x: x[13:16] == "CA ", dye_pdblines )[0]
N_list = filter(lambda x: x[13:16] == "N  ", dye_pdblines )
C_list = filter(lambda x: x[13:16] in "C  ", dye_pdblines )

if len(N_list) >0:
    N_new = N_list[0]
else:
    # we must search for N within 0.2 nm of CA...
    CA_new_xyz = makevec(CA_new)
    for line in dye_pdblines:
        tmp_xyz = makevec(line)
        dCAX = dist(tmp_xyz,CA_new_xyz)
        if dCAX < 2.0 and "N" in line[13:16]:
            N_new = line
    sys.stdout.write("Found N in dye : %s\n"%(N_new[13:16]))

if len(C_list) >0:
    C_new = C_list[0]
else:
    # we must search for C within 0.2 nm of CA, but which is not CB!
    CA_new_xyz = makevec(CA_new)
    for line in dye_pdblines:
        tmp_xyz = makevec(line)
        dCAX = dist(tmp_xyz,CA_new_xyz)
        if dCAX < 2.0 and "C" in line[13:16]:
            for line2 in dye_pdblines:
                tmp2_xyz = makevec(line)
                dCX = dist(tmp_xyz,tmp2_xyz)
                # must have oxygens close by...
                if dCX < 2. and "O" in line2[13:16]:
                    C_new =  line
    sys.stdout.write("Found C in dye : %s\n"%(C_new[13:16]))

CA_new_xyz = makevec(CA_new)
C_new_xyz = makevec(C_new)
N_new_xyz = makevec(N_new)
CA_old_xyz = makevec(CA_old)
C_old_xyz = makevec(C_old)
N_old_xyz = makevec(N_old)

M_old = make_rotmat(CA_old_xyz,C_old_xyz,N_old_xyz)
M_new = make_rotmat(CA_new_xyz,C_new_xyz,N_new_xyz)
M_old_T = numpy.transpose(M_old)

rot = numpy.dot( M_old_T, M_new ) # rotation matrix to apply to dye

dx = CA_old_xyz - CA_new_xyz # translation vector to apply to dye

#==================================================
# write first unchanged bit of protein ...
outp = open(pdbout,"w")

atom_idx = 0
res_idx = 0
pres = -1
for line in pre:
	tmp = int(line[21:26])
	if tmp>pres:
		pres = tmp
		res_idx += 1
	atom_idx+=1
	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))

#==================================================
# write dye coordinates including replaced amino acid residue
dyecrd  = []
for line in dye_pdblines:
	atom = line[12:16].strip()
	xyz = makevec(line)
	dr = xyz - CA_new_xyz 
	xyz_rottrans = CA_old_xyz + numpy.dot(rot,dr)
	dyecrd.append(xyz_rottrans)

res_idx += 1
for k in range(len(dye_pdblines)):
	line = dye_pdblines[k]
	x,y,z = dyecrd[k]
	atom_idx += 1
        outp.write("ATOM  %5i%6s%3s %5i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(atom_idx,line[11:17],line[17:20],res_idx,
		x,y,z,1.,0.))


#==================================================
# write C-ter remaining part of protein
pres = -1
for line in post:
	#outp.write(line)
	tmp = int(line[21:26])
	if tmp>pres:
		pres = tmp
		res_idx += 1
	atom_idx+=1
	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))

outp.close()

#==================================================


# first we check for clashes

#precrd = map(makevec,pre)
#postcrd = map(makevec,post)
#protcrd = precrd+postcrd
##dyecrd_raw = map(lambda x: x-SG_dye_xyz, map(makevec,dyelines))
##dyecrd = map(lambda x: SG_cys_xyz + numpy.dot(rot,x), dyecrd_raw)
#
#def isclash(protxyz,dyexyz):
#	rmin = 3.0
#	for k in dyexyz:
#		for l in protxyz:
#			if numpy.linalg.norm( l-k ) < rmin:
#				print l
#				print k
#				return 1
#	return 0
#
##==================================================
#outp = open(pdbout,"w")
#
#atom_idx = 0
#res_idx = 0
#pres = -1
#for line in pre:
#	tmp = int(line[21:26])
#	if tmp>pres:
#		pres = tmp
#		res_idx += 1
#	atom_idx+=1
#	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))
#
#cys_atoms = []
#
#
#res_idx += 1
##for line in cys:
##	atom = line[12:16].strip()
##	if atom == "HG":
##		continue
##	cys_atoms.append(atom)
##	atom_idx+=1
##	outp.write("ATOM  %5i%6s%4s%5i%s\n"%(atom_idx,line[11:17],dye,res_idx,line[26:66]))
#
#dyecrd_nocys  = []
#dyelines_nocys = []
#for line in dyelines:
#	atom = line[12:16].strip()
#	if atom in cys_atoms:
#		continue
#	#if atom in [ "O
#
#	dyelines_nocys.append(line)
#	xyz_old = makevec(line)
#	dr = xyz_old - SG_dye_xyz 
#	xyz_rottrans = SG_cys_xyz + numpy.dot(rot,dr)
#	dyecrd_nocys.append(xyz_rottrans)
#
#	#x,y,z = xyz_rottrans
#	#outp.write("ATOM  %5i%10s%5i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(atom_idx,line[11:21],res_idx,
#	#	x,y,z,1.,0.))
#
#ntry = 0
#dyecrd_nocys_save = []
#for xyz in dyecrd_nocys:
#	dyecrd_nocys_save.append(xyz)
#n_nocys = len(dyelines_nocys)
#
#maxtry = 20
#while ntry<maxtry and isclash(protcrd,dyecrd_nocys):
#	print "Trying to remove clashes by rotating chromophore"
#	ntry += 1
#	print "Try number %i of %i"%(ntry,maxtry)
#	zz = xnew
#	crossp = numpy.cross(sb_vec,zz)
#	xx = crossp / numpy.linalg.norm(crossp)
#	yy = numpy.cross(zz,xx)
#	print numpy.linalg.norm(xx)
#	print numpy.linalg.norm(yy)
#	print numpy.linalg.norm(zz)
#	rand_theta = 2.* math.pi * random.random()
#	costheta = math.cos(rand_theta)
#	sintheta = math.sin(rand_theta)
#	MM = numpy.array( [xx,yy,zz], numpy.float64 )
#	MM_T = numpy.transpose(MM)
#	rot1  = numpy.array( [ [ costheta, -sintheta, 0. ], \
#				[ sintheta, costheta, 0. ], \
#				[ 0., 0., 1. ] ], numpy.float64 )
#	
#	rot2 = numpy.dot( MM_T, rot1 )
#	rot = numpy.dot( rot2, MM )
#
#	dyecrd_nocys = []
#
#	for xyz_old in dyecrd_nocys_save:
#		dr = xyz_old - SG_cys_xyz 
#		xyz_rottrans = SG_cys_xyz + numpy.dot(rot,dr)
#		dyecrd_nocys.append(xyz_rottrans)
#	
#
#for k in range(n_nocys):
#	line = dyelines_nocys[k]
#	x,y,z = dyecrd_nocys[k]
#	atom_idx += 1
#	outp.write("ATOM  %5i%6s%4s%5i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(atom_idx,line[11:17],dye,res_idx,
#		x,y,z,1.,0.))
#
#pres = -1
#for line in post:
#	#outp.write(line)
#	tmp = int(line[21:26])
#	if tmp>pres:
#		pres = tmp
#		res_idx += 1
#	atom_idx+=1
#	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))
#
#
#outp.close()
#
#

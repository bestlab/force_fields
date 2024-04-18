#!/usr/bin/env python3

import sys,os,numpy

#T49_pdblines = list(filter(lambda x: x[0:6]=="ATOM  ", open("T49_C2R.pdb").readlines()))
template_pdblines = list(filter(lambda x: x[0:6]=="ATOM  ", open("A64_C2R.pdb").readlines()))

C2R_lines = list(filter(lambda x: x.find(' C2R ')>0, template_pdblines))

# T49
# C99 - O91 - C20 

# Cy3B
# C17 - O3 - C18

dye = "CF6"

dye_pdblines = list(filter(lambda x: x[0:6]=="ATOM  ", open(f"{dye}_GMX.pdb").readlines()))

outp = open(f"{dye}_C2R.pdb","w")

def get_xyz(line):
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    return numpy.array([x,y,z], numpy.float64)

#for line in dye_pdblines:
#    if line.find("C17")>0:
#        o_d = get_xyz(line)
#    elif line.find("O3")>0:
#        a_d = get_xyz(line)
#    elif line.find("C18")>0:
#        b_d = get_xyz(line)
#    outp.write(line)

#for line in template_pdblines:
#    if line.find("C99") > 0:
#        o_t = get_xyz(line)
#    elif line.find("O91")>0:
#        a_t = get_xyz(line)
#    elif line.find("C20")>0:
#        b_t = get_xyz(line)
for line in template_pdblines:
    if line.find("C99") > 0:
        o_t = get_xyz(line)
    elif line.find("O7")>0:
        a_t = get_xyz(line)
    elif line.find("C34")>0:
        b_t = get_xyz(line)

for line in dye_pdblines:
    if line.find("C10")>0:
        o_d = get_xyz(line)
    elif line.find("O4")>0:
        a_d = get_xyz(line)
    elif line.find("C23")>0:
        b_d = get_xyz(line)
    outp.write(line)



def define_basis(r_o,r_a,r_b):
    # origin =  r_o
    # x = r_oa/|r_oa|
    # y = r_oa x r_ob (normalized)
    # z = x x y
    r_oa = r_a - r_o
    r_ob = r_b - r_o
    x = r_oa / numpy.linalg.norm(r_oa)
    yy = numpy.cross(r_oa,r_ob)
    y = yy / numpy.linalg.norm(yy)
    z = numpy.cross(x,y)
    return numpy.array([x,y,z], numpy.float64 )

M_d = define_basis(o_d,a_d,b_d)
M_t = define_basis(o_t,a_t,b_t)

T = numpy.dot(M_d.transpose(),M_t)

index  = len(dye_pdblines)

for line in C2R_lines:
    index += 1
    r = get_xyz(line)
    rt = numpy.dot(T,r-o_t)+o_d
    outp.write("ATOM  %5i%19s%8.3f%8.3f%8.3f%s"%(index,line[11:30],rt[0],rt[1],rt[2],line[54:]))

outp.write("END\n")
outp.close()




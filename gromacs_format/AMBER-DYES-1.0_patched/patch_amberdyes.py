#!/usr/bin/env python

# fix incorrect structure of maleimide adducts in AMBER-DYES force field

import sys,os,numpy

oldffdir = "original/amber99sb_dyes.ff/"
newffdir = "amber99sb_dyes_patch.ff/"
oldPDBdir = "original/PDB_structures"
newPDBdir = "PDB_structures_patch"

if os.path.exists(newffdir):
    os.system("rm -r %s"%(newffdir))

if os.path.exists(newPDBdir):
    os.system("rm -r %s"%(newPDBdir))

os.system("cp -r %s %s "% (oldffdir,newffdir))
os.system("cp %s %s/ffbonded.itp "% ("ffbonded_patched.itp",newffdir))
os.mkdir(newPDBdir)

old_top = "%s/aminoacids.rtp" % (oldffdir)
new_top = "%s/aminoacids.rtp" % (newffdir)

# ======================================================================
# fix topology

old_top_raw = open(old_top).read()

old_top_blocks = old_top_raw.split("\n\n")

#sys.stdout.write("%sa\n"%(old_top_blocks[3]))

# improper_excl: delete any impropers containing these atoms
top_patches = {
        "C1C": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C3": ("c3g",-0.18), "H4": ("hcg",0.09), "C2": ("c3g",-0.09) }, \
                "new_bonds": [ ("C2","HX2"),("C3","HX3") ],
                "improper_excl": [ "C2", "C3", "C4", "N1", "O3", "O92" ] \
                },
        "C1R": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C8": ("c3g",-0.18), "H12": ("hcg",0.09), "C7": ("c3g",-0.09) }, \
                "new_bonds": [ ("C7","HX2"),("C8","HX3") ],
                "improper_excl": [ "C7", "C8", "C9", "N3", "O3", "O4" ] \
                },
        "E1N": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C6": ("c3g",-0.18), "H8": ("hcg",0.09), "C5": ("c3g",-0.09) }, \
                "new_bonds": [ ("C5","HX2"),("C6","HX3") ],
                "improper_excl": [ "C5", "C6", "C8", "N2", "O3", "O92" ] \
                },
        "C2C": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C7": ("c3g",-0.18), "H8": ("hcg",0.09), "C6": ("c3g",-0.09) }, \
                "new_bonds": [ ("C6","HX2"),("C7","HX3") ],
                "improper_excl": [ "C6", "C7", "C8", "N2", "O3", "O4" ] \
                },
        "C2R": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C8": ("c3g",-0.18), "H12": ("hcg",0.09), "C7": ("c3g",-0.09) }, \
                "new_bonds": [ ("C7","HX2"),("C8","HX3") ],
                "improper_excl": [ "C7", "C8", "C9", "N3", "O3", "O4" ] \
                },
        "E2N": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C6": ("c3g",-0.18), "H8": ("hcg",0.09), "C5": ("c3g",-0.09) }, \
                "new_bonds": [ ("C5","HX2"),("C6","HX3") ],
                "improper_excl": [ "C5", "C6", "C7", "N2", "O3", "O92" ] \
                },
        "C3C": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C7": ("c3g",-0.18), "H8": ("hcg",0.09), "C6": ("c3g",-0.09) }, \
                "new_bonds": [ ("C6","HX2"),("C7","HX3") ],
                "improper_excl": [ "C6", "C7", "C8", "N2", "O3", "O4" ] \
                },
        "C3R": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C8": ("c3g",-0.18), "H12": ("hcg",0.09), "C7": ("c3g",-0.09) }, \
                "new_bonds": [ ("C7","HX2"),("C8","HX3") ],
                "improper_excl": [ "C7", "C8", "C9", "N3", "O3", "O4" ] \
                },
        "E3N": { "new_atoms": { "HX2": ("hcg",0.09), "HX3": ("hcg",0.09), \
                "C6": ("c3g",-0.18), "H8": ("hcg",0.09), "C5": ("c3g",-0.09) }, \
                "new_bonds": [ ("C5","HX2"),("C6","HX3") ],
                "improper_excl": [ "C5", "C6", "C7", "N2", "O3", "O92" ] \
                }
        }

atomlist = []

# ======================================================================
def patch_residue(block,patch):
    bs = block.split('\n')
    #outp = [ bs[0], ]
    outp = bs[0]+"\n"
    #print outp
    atom_list = []
    section = "none"
    natom = 0
    qsum_orig = 0.
    existing_atoms = []
    for line in bs[1:]:
        #print line
        if line.strip()[0] == ';':
            #outp.append(line)
            outp+=line+"\n"
        elif line.strip()[0] == '[':
            #outp.append(line)
            section = line[line.find("[")+1:line[1:].find("]")].strip()
            if section == "bonds": # now dump new atoms
                print existing_atoms
                qsum_new = 0;
                print patch["new_atoms"].keys()
                pkeys = patch["new_atoms"].keys()
                pkeys.sort()
                for a in pkeys:   #patch["new_atoms"].keys():
                    qsum_new += patch["new_atoms"][a][1]
                    if a not in existing_atoms:
                        natom+=1
                        atom_list.append((a,patch["new_atoms"][a][0],patch["new_atoms"][a][1],natom))
                dqtot = qsum_orig - qsum_new
                print dqtot
                dq = dqtot/float(len(patch["new_atoms"].keys()))
                for atom in atom_list:
                    a,t,q,i = atom
                    if a in patch["new_atoms"].keys():
                        q+=dq
                    #outp.append("%6s %6s %12.6f %6i\n"%(a,t,q,i))
                    outp += "%6s %6s %12.6f %6i\n"%(a,t,q,i)
            if section == "impropers": # now dump new bonds
                for b in patch["new_bonds"]:
                    outp += "%6s%6s\n"%(b)
            print section
            outp+=line+"\n"

        else:
            if section == "atoms":
                natom+=1
                atom = line.split()[0]
                atype = line.split()[1]
                charge = float(line.split()[2])
                index = int(line.split()[3])
                if atom in patch["new_atoms"].keys():
                    existing_atoms.append(atom)
                    qsum_orig += charge
                    atom_list.append( (atom,patch["new_atoms"][atom][0],patch["new_atoms"][atom][1],natom) )
                else: 
                    atom_list.append( (atom,atype,charge,natom) )
            elif section == "impropers":
                a,b,c,d = line.split()
                if not ( a in patch["improper_excl"] or b in patch["improper_excl"] \
                        or c in patch["improper_excl"] or d in patch["improper_excl"]):
                            outp+=line+"\n"
            else:
                #outp.append(line)
                outp+=line+"\n"
    return outp

def patch_pdb(pdb_in, linker, patch,pdb_out):
    new_atoms = patch["new_atoms"]
    new_bonds = patch["new_bonds"]
    pdblines = filter(lambda x: x[0:6] == "ATOM  ", open(pdb_in).readlines())
    found = False
    new_written = False
    index = 0
    outp = open(pdb_out,"w")
    found_atoms = {}
    for line in pdblines:
        #index += 1
        res = line.split()[3]
        resnum = int(line.split()[4])
        print res, linker
        if res == linker:
            found = True
            atom = line[11:17].strip()
            #for bond in new_bonds:
            if atom in new_atoms.keys():
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                found_atoms[atom] = (x,y,z)
                # save crds and get name of other atom in bond to append!
            index += 1
            outp.write("ATOM  %5i%s\n"%(index,line[11:-1]))
        else:
            if found == True and new_written == False:
                pkeys = patch["new_atoms"].keys()
                pkeys.sort()
                for a in pkeys:
                    if a not in found_atoms.keys():
                        print a
                        for bond in new_bonds:
                            if a in bond:
                                for aa in bond:
                                    if aa != a:
                                        other=aa
                        x,y,z = found_atoms[other]
                        x+=.05
                        y+=.05
                        z+=.05
                        index+=1
                        outp.write("ATOM  %5i  %3s %3s%6i    %8.3f%8.3f%8.3f  1.00  0.00\n"\
                                %(index,a,linker,resnum-1,x,y,z))
                new_written = True
            index+=1
            outp.write("ATOM  %5i%s\n"%(index,line[11:-1]))
        print found_atoms
    if new_written == False:
            pkeys = patch["new_atoms"].keys()
            pkeys.sort()
            for a in pkeys:
                if a not in found_atoms.keys():
                    print a
                    for bond in new_bonds:
                        if a in bond:
                            for aa in bond:
                                if aa != a:
                                    other=aa
                    x,y,z = found_atoms[other]
                    x+=.5
                    y+=.5
                    z+=.5
                    index+=1
                    outp.write("ATOM  %5i  %3s %3s%6i    %8.3f%8.3f%8.3f  1.00  0.00\n"\
                            %(index,a,linker,resnum,x,y,z))
            new_written = True
    return





# ======================================================================

# patch residues in rtf


problems = top_patches.keys()
top_out = open(new_top,"w")
for block in old_top_blocks:
    bs = block.split('\n')
    key = bs[0][bs[0].find("[")+1:bs[0].find("]")].strip()
    if key in problems:
        # add HX2, HX3
        print key
        newblock = patch_residue(block,top_patches[key])
        #top_out.write("%s\n\n"%(block))
        top_out.write("%s\n\n"%(newblock))
    else:
        top_out.write("%s\n\n"%(block))

top_out.close()
        

# ======================================================================
# patch PDB files


## make Ala3
#amberhome="/u/best/lib/amber12"
#leap_inp = """
#source leaprc.ff03.r1
#pep = sequence { NALA ALA CALA } 
#saveAmberParm pep ala3.prmtop ala3.inpcrd
#savePdb pep ala3.pdb
#quit 
#""" 
#inp = open("leap.inp","w")
#inp.write(leap_inp)
#inp.close()
#
#os.system("%s/bin/tleap -s -f leap.inp"%(amberhome))

dyes = os.listdir(oldPDBdir)

if not os.path.exists("tmp_pdb"):
    os.mkdir("tmp_pdb")

for dye in dyes:
    os.mkdir("%s/%s"%(newPDBdir,dye))
    pdbs = os.listdir("%s/%s"%(oldPDBdir,dye))
    for pdb in pdbs:
        linker = pdb[-7:-4]
        dyex = pdb[0:3]
        if linker in problems:
            print dyex, linker
            if linker[2] =="N":
                pres = 1
            elif linker[2] =="R":
                pres = 2
            elif linker[2] =="C":
                pres = 3
            os.system("./build_dyes_amberdye.py %i %s/%s/%s ala3.pdb tmp_pdb/ala3_%s"%(pres,oldPDBdir,dye,pdb,pdb))
            #new_bonds = top_patches[linker]["new_bonds"]
            patch_pdb("tmp_pdb/ala3_%s"%pdb, linker, top_patches[linker],"tmp_pdb/pala3_%s"%pdb)

            command = """pdb2gmx -f tmp_pdb/pala3_%s -o tmp_pdb/pala3_%s.gro -p pala3_%s.top -ff %s <<EOF
1
EOF"""\
                    %(pdb,pdb,pdb,"amber99sb_dyes_patch")
            os.system(command)
            
            command = "grompp -p pala3_%s.top -c tmp_pdb/pala3_%s.gro -o tmp_pdb/pala3_%s.tpr -f miniconj.mdp" \
                    %(pdb,pdb,pdb)
            os.system(command)

            command = "mdrun -c tmp_pdb/pala3_%s_min.pdb -s tmp_pdb/pala3_%s.tpr" %(pdb,pdb)
            os.system(command)
            raw_crd = filter(lambda x: x[0:6]=="ATOM  ", open("tmp_pdb/pala3_%s_min.pdb"%(pdb)).readlines())
            dye_out = open("%s/%s/%s"%(newPDBdir,dye,pdb),"w")
            resind =1
            atom =0
            for line in raw_crd:
                res = line.split()[3]
                if res in [ dyex, linker ]:
                    atom+=1
                    dye_out.write("ATOM  %5i%10s%5i%s"%(atom,line[11:21],resind,line[26:]))
            dye_out.write("TER\nEND\n")
            dye_out.close()

        else:
            os.system("cp %s/%s/%s %s/%s/%s"%(oldPDBdir,dye,pdb, newPDBdir,dye,pdb))



# clean up

os.system("rm *pdb.top ener.edr md.log posre.itp traj.trr")

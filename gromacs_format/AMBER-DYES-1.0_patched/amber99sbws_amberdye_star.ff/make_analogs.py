#!/usr/bin/env python

import sys,math

# make amino acid analogs by replacing CA with H and
# adjusting charges using the usual rules...

inp_raw = sys.stdin.readlines()

noncomment = filter(lambda x: x[0] != ';', inp_raw)

newres = 1
atoms = []
bonds = []
impropers = []
dihedrals = []
backbone_atoms = [ 'N', 'H', 'CA', 'HA', 'C', 'O' ]
qsum_deleted = 0.

for line in noncomment:
	if len(line.strip())==0:
		newres = 1
		if len(atoms) > 0 and resname != 'GLY':
			sys.stdout.write('[ %s ]\n'%(analog))
			sys.stdout.write('[ atoms ]\n')
			for atom in atoms:
				if atom[0].find('HB')>=0:
					hbcharge = atom[2]
					lasthb = atom[0]
			idx = 0
			for x in atoms:
				idx += 1
				atom_name, atom_type, charge = x
				if atom_name == 'CB':
					charge += ( qsum_deleted - hbcharge )
					x = (atom_name, atom_type, charge)
				s=  x+(idx,)
				sys.stdout.write("%6s %6s   %12.6f  %5i\n"% s)
				if atom_name == lasthb:
					idx +=1
					if atom_name[-1].isdigit():
						hbi = int(atom_name[-1])
					else:
						hbi = 1
					extra_hb_name = "HB%i"%( hbi+1)
					s = (extra_hb_name, atom_type, charge, idx)
					sys.stdout.write("%6s %6s   %12.6f  %5i\n"% s)

			if len(bonds)>0:
				sys.stdout.write('[ bonds ]\n')
				for b in bonds:
					sys.stdout.write('%5s %5s\n'%b)
				sys.stdout.write('%5s %5s\n'%('CA', extra_hb_name))

			if len(impropers)>0:
				sys.stdout.write('[ impropers ]\n')
				for i in impropers:
					sys.stdout.write('%5s %5s %5s %5s'%(i[0:4]))
					if len(i)>4:
						for extra in i[4:]:
							sys.stdout.write(' %s '%(extra))
					sys.stdout.write('\n')
			if len(dihedrals)>0:
				sys.stdout.write('[ dihedrals ]\n')
				for i in dihedrals:
					sys.stdout.write('%5s %5s %5s %5s'%(i[0:4]))
					if len(i)>4:
						for extra in i[4:]:
							sys.stdout.write(' %s '%(extra))
					sys.stdout.write('\n')

			sys.stdout.write('\n\n')

				
				
				
				
		atoms = []
		bonds = []
		impropers = []
		dihedrals = []
		qsum_deleted = 0.
		continue
	if newres and len(line.strip())>0:
		resname = line[line.find('[')+1:line.find(']')].strip()
		analog = "A"+resname
		#sys.stdout.write('\n\n[ %s ]\n'%(analog))
		newres = 0
		continue
	if line.find('atoms')>0 or line.find('impropers')>0 or line.find('bonds') >0\
		or line.find('dihedrals') > 0:
		section = line[line.find('[')+1:line.find(']')].strip()
		continue
	if section == 'atoms':
		atom_name, atom_type, acharge, idx = line.split()[0:4]
		charge = float(acharge)
		if atom_name in backbone_atoms:
			qsum_deleted += charge
		else:
			atoms.append((atom_name,atom_type,charge))
	elif section == 'bonds':
		atom_a, atom_b = line.split()
		if atom_a in backbone_atoms or atom_b in backbone_atoms:
			continue
		else:
			bonds.append((atom_a,atom_b))
	elif section == 'impropers':
		#print line
		atom_a, atom_b, atom_c, atom_d = line.split()[0:4]
		if atom_a in backbone_atoms or atom_b in backbone_atoms \
			or atom_c in backbone_atoms or atom_d in backbone_atoms:
			continue
		else:
			if len(line.split())==4:
				impropers.append((atom_a,atom_b,atom_c,atom_d))
			else:
				#extras = line.split()[4:]
				impropers.append(tuple(line.split()))
	elif section == 'dihedrals':
		atom_a, atom_b, atom_c, atom_d = line.split()[0:4]
		if atom_a in backbone_atoms or atom_b in backbone_atoms \
			or atom_c in backbone_atoms or atom_d in backbone_atoms:
			continue
		else:
			if len(line.split())==4:
				dihedrals.append((atom_a,atom_b,atom_c,atom_d))
			else:
				#extras = line.split()[4:]
				dihedrals.append(tuple(line.split()))


				
		
		



		

	

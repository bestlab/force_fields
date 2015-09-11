#!/bin/csh -f

set pdbfile = "csp_cys_amb.pdb"

set ff = "amber03ws"
set water = "tip4p2005"

source ${HOME}/gromacs453/bin/GMXRC.csh

echo $pdbfile

#======================================================================
# step 1: import pdb into gromacs with ff03ws
#======================================================================

${GMXBIN}/pdb2gmx -f ${pdbfile} \
-p csp_vac.top \
-o csp_vac.pdb -ignh \
-ff ${ff} -nortpres <<EOF
1
EOF

#======================================================================
# step 2: add dyes with build_dyes.py script
#======================================================================

## put 488 on res 2
./build_dyes.py 2 488 csp_vac.pdb csp_s488.pdb
## put 488 on res 2 and 594 on res 68
./build_dyes.py 68 594 csp_s488.pdb csp_s488s594.pdb
## put 594 on res 2
./build_dyes.py 2 594 csp_vac.pdb csp_s594.pdb

set fileroot = "csp_dyes"

set dyepdb = "csp_s488s594.pdb"

#======================================================================
# step 3: create ff from generated pdb using pdbgmx
#======================================================================

${GMXBIN}/pdb2gmx -f ${dyepdb} \
-p ${fileroot}_vac.top \
-o ${fileroot}_vac.gro  \
-ff ${ff} -ter << EOF
1
EOF

#======================================================================
# step 4: set up water box and minimize
#======================================================================

set boxlen = 6.5
echo "writing box length"
${GMXBIN}/editconf -f ${fileroot}_vac.gro \
	-o ${fileroot}_ed.gro \
	-box $boxlen -bt octahedron

cp ${fileroot}_vac.top ${fileroot}.top

echo "generating water box"
${GMXBIN}/genbox -cp ${fileroot}_ed.gro \
	-cs amber03w.ff/tip4p2005.gro -o ${fileroot}_wat.gro \
	-p ${fileroot}.top

echo "preparing for water box minimization"
${GMXBIN}/grompp -v -f mini.mdp -c ${fileroot}_wat.gro \
	-o ${fileroot}_mini.tpr -p ${fileroot}.top

echo "minimizing water box"
${GMXBIN}/mdrun -pd -v -s ${fileroot}_mini.tpr \
	-o ${fileroot}_mini.trr \
	-c ${fileroot}_wat_min.gro \
	-e ${fileroot}_mini.ene \
	-g ${fileroot}_mini.log


#======================================================================
# if you get this far there shouldn't be any further issues to worry
# about!
#======================================================================



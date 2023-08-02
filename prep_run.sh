if [ $# -ne 1 ]; then
  echo "Usage: $0 protein.pdb"
  exit 1
fi

pdb_in=$1

rm *itp eq* md*

### PREPARING THE SYSTEM ###

# proceed the input pdb with pdbfixer to resolve possible 
pdbfixer $pdb_in --add-atoms heavy --replace-nonstandard --add-residues --keep-heterogens=none --output=aa_fixed.pdb

# get secondary structure using DSSP from MDtraj
ss=$(python dssp.py aa_fixed.pdb)
# martinize the system
martinize2 -maxwarn 100 -f aa_fixed.pdb -x _cg.pdb -o _topol.top -scfix -cys auto  -elastic -p backbone -ss $ss -nt -merge A,B,C,D,E,F,G,H,I
# create the cubic box
gmx editconf -f _cg.pdb -o _cg_box.gro -d 2 -c -bt cubic
# add water and probe molecules
python2 insane_probes.py -f _cg_box.gro -o _cg_sol.gro -p _insane.top -salt 0 -charge auto -sol W:0.952 -sol PHEN:0.003 -sol ACET:0.009 -sol IPA:0.009 -sol DMAD:0.009 -sol IPO:0.009 -sol PPN:0.009 -pbc keep
# create the topology file
echo '#include "_martini/martini_v3.0.0.itp"
#include "_martini/martini_v3.0.0_ions_v1.itp"
#include "_martini/martini_v3.0.0_probes_v1.itp"
#include "_martini/martini_v3.0.0_solvents_v1.itp"
#include "molecule_0.itp"

[ system ]
Title of the system

[ molecules ]
molecule_0    1' > topol.top

sed -n -E '/^(CLBZ|PHEN|BENZ|ACE|IPA|DMAD|IPO|PPN|W)/p' _insane.top >> topol.top

### RUNNING SIMULATIONS ###

# run EM
gmx grompp -f _mdp/cg_em.mdp -p topol.top -c _cg_sol.gro -r _cg_sol.gro -o em.tpr -maxwarn 10
gmx mdrun -deffnm em -v
# run NVT eq
gmx grompp -f _mdp/cg_eq_1.mdp -p topol.top -c em.gro -r em.gro -o eq_1.tpr -maxwarn 10
gmx mdrun -deffnm eq_1 -v
# run NpT eq
gmx grompp -f _mdp/cg_eq_2.mdp -p topol.top -c eq_1.gro -r em.gro -o eq_2.tpr -maxwarn 10
gmx mdrun -deffnm eq_2 -v
# run production
gmx grompp -f _mdp/cg_md.mdp -p topol.top -c eq_2.gro -r em.gro -o md.tpr -maxwarn 10
gmx mdrun -deffnm md

### PROCESSING TRAJECTORY ###

# create index file
gmx select -s md.tpr -select "not resname W" -on not_w.ndx
# perform rot+trans fit
gmx trjconv -f md.xtc -o md_fixed_fit.xtc -fit rot+trans -s md.tpr << EOF
1
0
EOF
# fix PBC
gmx trjconv -f md_fixed_fit.xtc -o md_fixed_pbc.xtc -pbc mol -s md.tpr -n not_w.ndx << EOF
0
EOF
# dump the fist frame
gmx trjconv -f md_fixed_fit.xtc -o md_fixed_pbc.pdb -pbc mol -s md.tpr -n not_w.ndx -dump 0 << EOF
0
EOF

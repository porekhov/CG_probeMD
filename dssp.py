import mdtraj as md
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('pdb', help='PDB file')
args = parser.parse_args()

str_ini = md.load(args.pdb, top=args.pdb)
ss = md.compute_dssp(str_ini, simplified=True)[0]
print(''.join(ss))

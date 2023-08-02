import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.density import DensityAnalysis
from gridData import Grid
from sklearn.cluster import MeanShift
import numpy as np

# bandwidth for MeanShift
bandwidth=6

# list with elements: [ probe specie(s) for density calculation ], bulk_density
probes = [
[['PHEN', 'ACET', 'IPA', 'DMAD', 'IPO', 'PPN'], 0.0003486],
[['PHEN'], 5.718e-05],
[['ACET'], 5.827e-05],
[['IPA'], 5.827e-05],
[['DMAD'], 5.827e-05],
[['IPO'], 5.827e-05],
[['PPN'], 5.827e-05]
]

Tref = 303.15 # temperature
grid_d = 2 # grid spacing
kt = 0.0083188*Tref # kT value for cut-off

def pdb_line(at_num, at_type, res_type, chain_id, res_num, xyz, occ, temp):
    # returns a pdb-formatted line
    return "%6s%5s %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
            % ("ATOM  ", at_num, at_type, res_type, chain_id, res_num, xyz[0], xyz[1], xyz[2], occ, temp)

# align the AA model to the CG model

ref = mda.Universe('md_fixed_pbc.pdb')
ref.select_atoms('name BB').names = 'CA'
mobile = mda.Universe('aa_fixed.pdb')  

align.alignto(mobile, ref, select='name CA and not altloc B C D E', match_atoms=False)
mobile.select_atoms('protein').write('aa_aligned.pdb')

# density calculation

u = mda.Universe('md_fixed_pbc.pdb', 'md_fixed_pbc.xtc')

protein = u.select_atoms('protein')
protein_COM = protein.center_of_mass()

for p in probes:
    p_name = p[0]
    n0 = p[1]
    p_sel = u.select_atoms('not name VS* and (resname ' + ' '.join(p_name) + ' and around 6 protein)', updating=True)

    D = DensityAnalysis(p_sel, delta=grid_d, gridcenter=protein_COM, \
                        xdim=u.dimensions[0], ydim=u.dimensions[1], zdim=u.dimensions[2])
    D.run()

    grid = D.results.density.grid
    grid = np.log(grid/n0)
    grid[np.isneginf(grid)] = 0 # remove inf due to log(0)
    grid = (-8.314*Tref/1000)*grid # density -> kJ/mol
    grid = grid*(grid < -1*kt)
    D.results.density.grid = grid
    D.results.density.export('dens_' + '_'.join(p_name) + '.dx')
    print('Wrote file for: ', ' '.join(p_name), ', min/max values:', grid.min(), grid.max())

for p in probes:
    p_name = p[0]
    
    G = Grid('dens_' + '_'.join(p_name) + '.dx')
    
    grid = G.grid
    edges = G.edges
    hotspots, energies = [], []

    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            for k in range(grid.shape[2]):
                if grid[i,j,k] < 0:
                    hotspots.append([edges[0][i], edges[1][j], edges[2][k]])
                    energies.append(grid[i,j,k])
                
    hotspots = np.array(hotspots)
    energies = np.array(energies)
    
    clustering = MeanShift(bandwidth=bandwidth).fit(hotspots)
    
    clusters = {}
    for i in np.unique(clustering.labels_):
        clusters[i] = [hotspots[clustering.labels_ == i].shape[0], \
                   energies[clustering.labels_ == i].mean(axis = 0), \
                   energies[clustering.labels_ == i].min(axis = 0), \
                   hotspots[clustering.labels_ == i][energies[clustering.labels_ == i].argmin()]]


    clusters_ene_sort = [i[0] for i in sorted(clusters.items(), \
                         key=lambda item: item[1][1], reverse=False)]
    clusters_size_sort = [i[0] for i in sorted(clusters.items(), \
                          key=lambda item: item[1][0], reverse=True)]

    # pdb_line(at_num, at_type, res_type, chain_id, res_num, x, y, z, occ, temp)
    f_out = open('dens_'+ '_'.join(p_name) +'_clusters_ene.pdb', 'w')
    for i, j in enumerate(clusters_ene_sort):
        f_out.write(pdb_line(i, 'CLR', 'CLR', 'A', i + 1, clusters[j][3], clusters[j][1], clusters[j][2]) + '\n')
    f_out.close()

    f_out = open('dens_'+ '_'.join(p_name) +'_clusters_size.pdb', 'w')
    for i, j in enumerate(clusters_size_sort):
        f_out.write(pdb_line(i, 'CLR', 'CLR', 'A', i + 1, clusters[j][3], clusters[j][0], clusters[j][1]) + '\n')
    f_out.close()
    
    f_out = open('dens_'+ '_'.join(p_name) +'_clusters.pdb', 'w')
    for i, j in enumerate(clusters_ene_sort):
        xyz = hotspots[clustering.labels_ == i]
        ene = energies[clustering.labels_ == i]
        for spot in range(len(ene)):
            f_out.write(pdb_line(i, 'CLR', 'CLR', 'A', i + 1, xyz[spot, :], clusters[j][0], clusters[j][1]) + '\n')
    f_out.close()

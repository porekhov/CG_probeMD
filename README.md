# CG_probeMD
a pipeline for running MD simulations in the presence of probe molecules for druggability assessment.

A fully automated Colab pipeline is coming soon.

Requirements:
- Gromacs (version 2022.2 and later were tested), https://www.gromacs.org/

- martinize2, https://github.com/marrink-lab/vermouth-martinize

- MDTraj, https://www.mdtraj.org/1.9.5/installation.html

- pdbfixer, https://github.com/openmm/pdbfixer

- MDAnalysis, https://www.mdanalysis.org/pages/installation_quick_start/

- griddataformats, https://www.mdanalysis.org/GridDataFormats/installation.html

- numpy, sklearn

Folders _martini and _mdp contain the forcefield parameters and topologies (including probe molecules) and Gromacs simulation parameters, respectively.

dssp.py is a simple wrapper for the mdtraj dssp function.

insane_probes.py is a slightly modified insane.py script with additional solvents (probes).

How to use:

1. Run prep_run.sh providing a PDB file with the protein of interest (may require pretreatment), ./run.sh protein.pdb

2. After the production simulation is finished, run calc_dens_clusters.py. If finished normally, several files should appear in the current directory:

 - dens_XXX.dx - densities for individual probes and all probes together. densities are recalculated into free energies of corresponding probe molecules, i.e., in kJ/mol units.

 - dens_XXX_clusters.pdb - clustered density, each grid dot is labeled with its cluster number in the residue field, cluster size (i.e., number of grid dots belonging to the cluster) in the occupancy field, and mean cluster energy in the temperature factor field.

 - dens_XXX_clusters_ene.pdb - clusters sorted by the mean free energy (from the lowest to the highest). Each cluster is represented by a single hotspot corresponding to the lowest energy observed for this cluster. Mean energy is in the occupancy field and the lowest observed energy is in the temperature factor field.

 - dens_XXX_clusters_size.pdb - clusters sorted by the cluster size (from the biggest to the smallest). Each cluster is represented by a single hotspot corresponding to the lowest energy observed for this cluster. Cluster size (i.e., the number of grid dots belonging to the cluster) is in the occupancy field, and mean cluster energy is in the temperature factor field.
 
3. Densities can be visualized by VMD, Pymol, or using the provided visualize.ipynb notebook.

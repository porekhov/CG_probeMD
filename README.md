# Colabind: A Cloud-Based Approach for Prediction of Binding Sites Using Coarse-Grained Simulations with Molecular Probes

a pipeline for running and analyzing co-solvent MD simulations in the presence of probe molecules for druggability assessment.

**Colab-based OpenMM version**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/porekhov/CG_probeMD/blob/main/colabind_probeMD.ipynb)
Run the pipeline in Google Colab using OpenMM.

**Local OpenMM version** 

Run the pipeline locally using OpenMM. Download and run notebook colabind_local.ipynb

Requirements:

- martinize2, https://github.com/marrink-lab/vermouth-martinize

- MDTraj, https://www.mdtraj.org/1.9.5/installation.html

- pdbfixer, https://github.com/openmm/pdbfixer

- MDAnalysis, https://www.mdanalysis.org/pages/installation_quick_start/

- griddataformats, https://www.mdanalysis.org/GridDataFormats/installation.html

- numpy, sklearn

The working directory should contain the following files from this repository:

- insane_probes.py, slightly modified insane.py script with additional solvents (probes)
- _martini, the forcefield parameters and topologies (including probe molecules)


**Local Gromacs version**

Alternatively, the same pipeline can be run using Gromacs. Download and run notebook colabind_local_gmx.ipynb



Requirements:

- Same as above + Gromacs (version 2022.2 and later were tested), https://www.gromacs.org/

The working directory should also contain the following files from this repository:

- Folders _martini and _mdp, forcefield parameters and topologies (including probe molecules) and Gromacs simulation parameters, respectively

- insane_probes.py, slightly modified insane.py script with additional solvents (probes)

Output files:

- dens_XXX.dx - densities for individual probes and all probes together. densities are recalculated into free energies of corresponding probe molecules (in kJ/mol units).

- centers_mean_ene.pdb Clusters sorted by their mean free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster minimal energy; beta factor = the cluster mean energy.

- centers_min_ene.pdb Clusters sorted by their minimal free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster mean energy; beta factor = the cluster minimal energy.

- centers_sum_ene.pdb Clusters sorted by their total free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster mean energy; beta factor = the cluster total energy.

- clusters_mean_ene.pdb Clusters sorted by their mean free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster mean energy.

- clusters_min_ene.pdb Clusters sorted by their minimal free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster mininal energy.

- clusters_sum_ene.pdb Clusters sorted by their total free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster total energy.

- File aa_aligned.pdb is an all-atom model aligned to the coarse-grained one (it is stored in _cg.pdb), which can be used as a reference structure for visualization and interpretation of the results.
 
3. Densities and clusters can be visualized by VMD, Pymol, or using the provided visualize.ipynb notebook.

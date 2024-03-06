# Colabind: A Cloud-Based Approach for Prediction of Binding Sites Using Coarse-Grained Simulations with Molecular Probes

a pipeline for running and analyzing co-solvent MD simulations in the presence of probe molecules for druggability assessment.

**Computational workflow**

a. The input structure of a target in the PDB format is pre-processed using pdbfixer36 to repair possible minor issues with missing atoms, etc.;
b. The coarse-grained topology in the GROMACS format (.top and .itp files) is built using martinize237 with the most recent MARTINI 3 version28 together with the ElNeDyn approach38 for maintaining protein spatial structure;
c. The coarse-grained target structure is placed in a cubic box with  the dimensions chosen such that the minimal distance between the periodic replicas of the solute is at least 3.0 nm, and solvated with the standard MARTINI water particles and the probe compounds at the desired concentrations using the slightly modified insane.py script;
d. The molecular dynamics simulations are run using the OpenMM toolkit augmented with the martini_openmm package: (1) minimization; (2) NVT equilibration; (3) NpT equilibration; (4) production simulation;
e. The molecular dynamics trajectory of the production simulation is processed in order to align the target to the first frame, wrap the probe molecules around it, and fix the periodic boundary conditions using MDAnalysis;
f. The densities are calculated at the 3D grid around the aligned target with the voxel size of 2 Å using the gridData module of MDAnalysis;
g. The obtained densities, Dprobe, (individual densities for each type of probe compounds and the integrated “joint” density calculated for all probes together considering all probe CG particles as equal) are converted to the “free energy” densities using the Boltzmann inversion and the bulk densities, Dbulk, estimated using one of two approaches:
from the bulk solution with the same composition and the water:probe ratio (default);
based on the relative number of solute molecules to water molecules with water assigned a concentration of 55 M (following 10.1073/pnas.2214024119);
h. The “free energy” densities are filtered using < −1 kT criterion as a threshold for “significant” voxels, which we term interaction “hotspots”;
i. The densities are clustered using the MeanShift algorithm with the bandwidth parameter by default set to 6.


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

**Output files**

- dens_XXX.dx - densities for individual probes and all probes together. densities are recalculated into free energies of corresponding probe molecules (in kJ/mol units).

- centers_mean_ene.pdb Clusters sorted by their mean free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster minimal energy; beta factor = the cluster mean energy.

- centers_min_ene.pdb Clusters sorted by their minimal free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster mean energy; beta factor = the cluster minimal energy.

- centers_sum_ene.pdb Clusters sorted by their total free energy (from the lowest to the highest). Each cluster is represented by a single hotspot with the mininal energy. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = the cluster mean energy; beta factor = the cluster total energy.

- clusters_mean_ene.pdb Clusters sorted by their mean free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster mean energy.

- clusters_min_ene.pdb Clusters sorted by their minimal free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster mininal energy.

- clusters_sum_ene.pdb Clusters sorted by their total free energy (from the lowest to the highest). Each cluster is represented by all belonging hotspots. Residue name = probe type (ALL for the joint density); Residue number = cluster rank; occupancy = hotspot energy; beta factor = the cluster total energy.

- File aa_aligned.pdb is an all-atom model aligned to the coarse-grained one (it is stored in _cg.pdb), which can be used as a reference structure for visualization and interpretation of the results.
 
3. Densities and clusters can be visualized by VMD, Pymol, or using the provided visualize.ipynb notebook.

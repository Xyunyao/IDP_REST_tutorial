# IDP Replica Exchange with Solute Scaling (REST2) simulation tutorial

<!-- <p align="center">
<br>
<a href=https://doi.org/10.5281/zenodo.14799045><img alt="DOI" src="https://zenodo.org/badge/DOI/10.5281/zenodo.14799045.svg"></a>
</p> -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14799045.svg)](https://doi.org/10.5281/zenodo.14799045)

## Overview

In this tutorial, you will learn how to implement replica exchange with solute scaling(REST2)[^1] enhanced sampling technique. This is an optimized implementation of replica exchange with solute tempering[^2] which is a variant of Hamiltonian replica exchange[^3].This repository provides the required files and helper scripts which will help you through the exciting journey. The basic outlook of the [file tree](#file-tree) is provided.

The repository is divided into three sections :

- **INSTALLATION_INSTRUCTIONS** : Before you start the simulations [this section](./INSTALLATION_INSTRUCTIONS/) helps in setting up and installing the required softwares.

- **TUTORIAL_FILES** : [This section](./TUTORIAL_FILES/) contains the files required for the tutorial. Reference files are also provided so that you can compare them against the files that you will generate.
  
- **POST_SIMULATION_ANALYSIS** : [This section](./POST_SIMULATION_ANALYSIS/) have the files required to perform the analysis on the REST2 simulations. The [simulated files](https://doi.org/10.5281/zenodo.14799045) can be downloaded and used for the analysis too.  

## Goals

- [ ] Setup python environment.
- [ ] Install Gromacs and Plumed
- [ ] Setup the simulation replicas.
- [ ] Equilibrate the replicas.
- [ ] Setup REST2 simulation replicas.
- [ ] Simulate replica exchange simulations.
- [ ] Demultiplexing and correction of periodic boundary conditions.
- [ ] Analysis of convergence of our simulation
- [ ] Cheers for the successful completion of REST2 simulation :tada:  

## File tree

<details>

<summary>File Tree:</summary>

```bash
.
├── INSTALLATION_INSTRUCTIONS
│   ├── README.md
│   ├── install.sh
│   ├── options.sh
│   └── setup_env.yml
├── LICENSE
├── POST_SIMULATION_ANALYSIS
│   ├── README.md
│   ├── analysis.ipynb
│   ├── demux.fix.pl
│   ├── make_demux.sh
│   ├── reference_files
│   │   ├── README.md
│   │   ├── replica_index.n50.s0.-80.xvg
│   │   ├── replica_index.xvg
│   │   └── replica_temp.xvg
│   └── scripts
│       ├── __pycache__
│       │   └── small_utilities.cpython-311.pyc
│       └── small_utilities.py
├── README.md
└── TUTORIAL_FILES
    ├── README.md
    ├── initial_input_files
    │   ├── README.md
    │   ├── processed.top
    │   ├── prot_only.pdb
    │   └── topol.top
    ├── mdp_files
    │   ├── NPT0.mdp
    │   ├── NPT1.mdp
    │   ├── NVT.mdp
    │   ├── README.md
    │   ├── minimz.mdp
    │   └── prod.mdp
    ├── reference_files
    │   ├── README.md
    │   ├── equilibrated_configurations
    │   │   ├── 0.gro
    │   │   ├── 1.gro
    │   │   ├── 2.gro
    │   │   ├── 3.gro
    │   │   ├── 4.gro
    │   │   ├── 5.gro
    │   │   ├── 6.gro
    │   │   ├── 7.gro
    │   │   ├── 8.gro
    │   │   └── 9.gro
    │   ├── og_topology
    │   │   ├── REST.top
    │   │   └── processed.top
    │   ├── scaled_topologies
    │   │   ├── topol_0.top
    │   │   ├── topol_1.top
    │   │   ├── topol_2.top
    │   │   ├── topol_3.top
    │   │   ├── topol_4.top
    │   │   ├── topol_5.top
    │   │   ├── topol_6.top
    │   │   ├── topol_7.top
    │   │   ├── topol_8.top
    │   │   └── topol_9.top
    │   └── solvated_protein
    │       ├── solv_ions_0.gro
    │       ├── solv_ions_1.gro
    │       ├── solv_ions_2.gro
    │       ├── solv_ions_3.gro
    │       ├── solv_ions_4.gro
    │       ├── solv_ions_5.gro
    │       ├── solv_ions_6.gro
    │       ├── solv_ions_7.gro
    │       ├── solv_ions_8.gro
    │       └── solv_ions_9.gro
    └── run_md.sh

14 directories, 61 files
```


</details>


## References

[^1]: Lingle Wang, Richard A. Friesner, and B. J. Berne, Replica Exchange with Solute Scaling: A More Efficient Version of Replica Exchange with Solute Tempering (REST2), The Journal of Physical Chemistry B 115 (30), 9431-9438, [DOI: 10.1021/jp204407d](https://doi.org/10.1021/jp204407d) (2011)

[^2]: P. Liu, B. Kim, R.A. Friesner, & B.J. Berne, Replica exchange with solute tempering: A method for sampling biological systems in explicit water, Proc. Natl. Acad. Sci. U.S.A. 102 (39) 13749-13754, [DOI: 10.1073/pnas.0506346102](https://doi.org/10.1073/pnas.0506346102) (2005).

[^3]: Hiroaki Fukunishi, Osamu Watanabe, Shoji Takada; On the Hamiltonian replica exchange method for efficient sampling of biomolecular systems: Application to protein structure prediction. J. Chem. Phys. 22 May 2002; 116 (20): 9058–9067. [DOI: 10.1063/1.1472510](https://doi.org/10.1063/1.1472510)


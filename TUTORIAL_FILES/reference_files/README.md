# Reference files

A set of reference files are provided here which will be used in the tutorial and also compare the files while performing the tutorial. 
- `equilibrated_configurations` contains the gromacs structure files of replicas after a successful equilibration at NVT ensemble.
    ```bash
    ./equilibrated_configurations/
    ├── 0.gro
    ├── 1.gro
    ├── 2.gro
    ├── 3.gro
    ├── 4.gro
    ├── 5.gro
    ├── 6.gro
    ├── 7.gro
    ├── 8.gro
    └── 9.gro
    ```
- `og_topology` has two files in it. `processed.top` topology file contains the information after the protein box is solvated and respective ions are added which you can see at `[ molecules ]` section.
    ```bash
    [ molecules ]
    ; Compound        #mols
    Protein_chain_A     1
    SOL         8763
    NA         8
    ```
    This topology file will be used as an input to `plumed` to generate scaled topology files for respective replicas.
    ```bash
    ./og_topology/
    ├── REST.top
    └── processed.top
    ```
- `scaled_topologies` contains the topology files which are scaled w.r.t to the temperature ladder.
    ```bash
    ./scaled_topologies/
    ├── topol_0.top
    ├── topol_1.top
    ├── topol_2.top
    ├── topol_3.top
    ├── topol_4.top
    ├── topol_5.top
    ├── topol_6.top
    ├── topol_7.top
    ├── topol_8.top
    └── topol_9.top
    ```
- `solvated_protein` contains the structure files with solvated protein and respective ions required for neutralization of the system.
    ```bash
    ./solvated_protein/
    ├── solv_ions_0.gro
    ├── solv_ions_1.gro
    ├── solv_ions_2.gro
    ├── solv_ions_3.gro
    ├── solv_ions_4.gro
    ├── solv_ions_5.gro
    ├── solv_ions_6.gro
    ├── solv_ions_7.gro
    ├── solv_ions_8.gro
    └── solv_ions_9.gro
    ```
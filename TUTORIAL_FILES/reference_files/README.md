# Reference files

We provide you with a set of reference files which you can use to compare the files while performing the tutorial. 
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
- `og_topology` has topology file after the protein box is solvated and respective ions are added and topology file used as a input to `plumed` to generate scaled topology files for respective replicas.
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
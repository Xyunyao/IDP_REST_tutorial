# Installation instructions

## Setup python environment

To setup the python environment required for the python scripts used for the analysis, use `tutorial_env.yml` which creates a python environment called `REST_tutorial` and installs all the modules required. Follow this [link](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) to install conda in your in your PC or follow [link](https://docs.python.org/3/library/venv.html) to create a python virtual env.

Following conda installation run the following commands to create and activate the tutorial environment:

```bash
conda env create -f tutorial_env.yml
conda activate REST_tutorial
```
or 

Following the creation and activation of a python=3.11 venv environment:
```bash
pip install -r tutorial_env.yml
 
```

## Installing gromacs patched with Plumed2

`install.sh` has a working script to install gromacs 2022.5 and plumed2 v2.9.1 in the users `$HOME` folder under `$HOME/opt`.
From this directory simply perform these two commands to install:
```bash
chmod +x install.sh
sh install.sh
```

Feel free to modify the installation directory. Be aware, if you change the installation directory to a location requiring elevated privileges the script will need to be performed via sudo, albiet not advised. 

```bash
chmod +x install.sh
sudo sh ./install.sh
```
To check if the installation is successful type `gmx -h`

## Goals

- [x] Setup python environment.
- [x] Install Gromacs and Plumed
- [ ] Setup the simulation replicas.
- [ ] Equilibrate the replicas.
- [ ] Setup REST2 simulation replicas.
- [ ] Simulate replica exchange simulations.
- [ ] Demultiplexing and correction of periodic boundary conditions.
- [ ] Analysis of convergence of our simulation
- [ ] Cheers for the successful completion of REST2 simulation :tada: 

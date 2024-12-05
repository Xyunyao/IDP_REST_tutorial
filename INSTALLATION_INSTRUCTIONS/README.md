# Setup instructions
## Installing gromacs patched with Plumed2

`install.sh` has a working script to install gromacs 2024.3 and plumed2 v2.9 in the users `$HOME` folder under `$HOME/opt`.
From this directory simply perform the two commands to install. 
```bash
chmod +x install.sh
sh install.sh
```

Feel free to modify the installation directory. Be aware, if you change the installation directory to a location requiring access above the users privileges the script will need to be performed via sudo. 

```bash
chmod +x install.sh
sudo sh ./install.sh
```

## Setup python environment

To setup the python environment required for post analysis python scripts use `setup_env.yml` which creates a python environment called `crazy_env` and installs all the modules required. Follow this [link](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) to install conda in your in your PC.
For that run the following commands:

```bash
conda env create -f setup_env.yml
```
or
```bash
pip install -r setup_env.yml
```

To activate the environment,

```bash
conda activate crazy_env
```

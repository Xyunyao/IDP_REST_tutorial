# IDP Replica Exchange with Solute Tempering (REST2)  Simulation Tutorial

## Required before proceeding
Before proceeding to running simulations, please make sure to install the required software as mentioned in [INSTALLATION_INSTRUCTIONS](./INSTALLATION_INSTRUCTIONS/). 
It is important to note, GROMACS must be compiled with mpi, patched with plumed and the executable must be:
```bash
gmx
```

## Simulation Helper Script
The provided helper script requires the executable to have the correct name. If you have a different naming scheme, such as a unique postfix, you are invited to modify the given script `run_md.sh`. For example,
```bash
sed -i "s/ gmx / gmx_mpi /g" run_md.sh
```


```bash
run_md.sh -h or --help   # provides all options 

Usage: usage [Options]
Options:
 -h, --help     Display this help message
 -v, --verbose  Enable verbosity
 -p, --topol  Topology File name for rest
 -s, --stage   Select current phase of simulation: Setup, Min, NVT, NPT, REST
 -l, --log      STDIO log File
 -strip, --strip Stripped TPR file excluding waters
```
If you desire to see verbose output, commands printed to screen and verbose output from gromacs use the `-v or --verbose` flags. 
Provided is a flag to log all output from the script, `-l <logname> or --log=<logname>`.

## Running Simulations
The preparation and simulations are separated into stages. 
Stages need to be performed in order:
*  Setup
*  Equilibration
  *  Minimization
  *  Thermalization and NVT equilibration
  *  NPT Equilibration
*  Replica Exchange Molecular Dynamics

### Setup
The helper script will setup the system with protein, water and ions with this command:
```bash
run_md.sh -s setup
```
This will produce the system topology and system.gro. 

### Equilibration
Equilibration stages:
*  Minimization
*  NVT
*  NPT
These stages can be performed with these series of commands:
```bash
# Minimize the system using system.gro and topol.top produced in the setup stage
run_md.sh -s Min 
# Thermalize and equilibrate system under the NVT ensemble using min.gro and topol.top
run_md.sh -s NVT
# Equilibrate the pressure and volume of the system using nvt.gro and topol.top
run_md.sh -s NPT
```
Test for pressure and volume convergence with `gmx energy`.

### Replica Exchange Simulations

#### Creating the Replica Exchange Base topology which will be used for scaling

Once well equilibrated creating a rest2 edited topology file is required to setup the rest directories and produce the run input files, `prod.tpr`. This edited topology will contain subscripts appended to the atom types of the molecule of interest, the protein. The `-s NPT` stage will produce the file `processed.top`. Before proceeding with the REST2 simulations modify the atom types with the command below. `start_line` and `end_line` are defined as the line number before and after the proteins `[atoms]` section, respectively. Using vim can make this process easier.
```bash
vim processed.top
```
Once `processed.top` is open in vim type the following command including the semicolon.
```vim
:set nl
```
Identify the start and endline and run the following bash command replacing the words start and end with the appropriate line numbers.
```bash
awk -v start_line=start -v end_line=end '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top
```

See the [REST.top](./reference_files/og_topology/REST.top) as a reference.

#### Running replica exchange simulations
You are strongly encouraged to have at least 2 GPUS; 5 or 10 GPUS would be best. If you are running **CPU only** simulations you are encouraged to have a minimum of 100 CPU cores, 10 cores per simulation, expect this to take a long time a month or two. If you are running **GPU** simulations you are encouraged to have a minimum of 2 cpus per replica, 20 CPUS total, 8 or 16 CPUS per simulation is better. With 2 GPUs the simulation will take on the order of 10 days, while using 10 GPUS will take approximately 3 days (tested on nvidia A5500 GPUs).

Running the REST simulations is simple:
```bash
run_md.sh -s REST -t REST.top
```

Questions to answer:
* Upon completing the basic tutorial answer these questions.
  * How are the acceptance ratios between replicas? How do you know if these acceptance ratios are good? 
  * Do you have good similarity between contact matrices? If not how do you correct for this issue? 
  * What differences do you observe between temperature replicas, as well as demultiplexed replicas? 
  * Are your simulations converged and how do you know? 
  * Try an observable of your interest, does this observable look converged, if not how can you overcome this discrepancy? 
  * Lastly, how do your simulations compare to those provided in the Book Chapter? If you do not have access the trajectories are provided on Zenodo. 

* After completing the basic tutorial, rerun the tutorial with 10 replicas each passing through the equilibration phases (you will need to modify the run_md.sh script or perform this by hand), you can start with your neutralized-solvated system for each replica, `system.gro`. 
  * What differences in final box size do you observe at the end of the NPT1 stage? 
  * Take the average box length for the last 10 ns of each replica's NPT1 stage. Use this average and change the box size of each replica, `gmx editconf`, thus making their volumes equal. Run rest again with these 10 equilibrated, equal volume systems. Now perform the same analyses and compare your results. What differences do you see and what contributes to these differences? 
  

# IDP Replica Exchange with Solute Tempering (REST2)  Simulation Tutorial

## Required before proceeding
Before proceeding to running simulations, please make sure to install the required softwares as mentioned in [INSTALLATION_INSTRUCTIONS](./INSTALLATION_INSTRUCTIONS/). 
<!-- It is important to note, GROMACS must be compiled with mpi, patched with plumed and the executable must be `gmx`. -->

> [!NOTE] 
>  GROMACS must be compiled with mpi, patched with plumed and the executable must be `gmx`.
 
<!-- ```bash
gmx
``` -->

## Performing Simulations with the helper Script
Below are two guides for producing input replicas and running REST2 simulations:
* A helper script ([run_md.sh](./run_md.sh)) with flags to handle gromacs commands and perform simulation, an aid for those with less experience. 
* Command-by-command with descriptions of each step for the more advanced user and should be followed when studying new systems. 

Before continuing forward with the provided helper script the user must note the name of the gromacs executable. If the executable is not simply `gmx`, the naming in the script must be substituted with the correct name. If you have a different naming scheme, such as a unique postfix as described in the GROMACS manual's Installation Section, you are invited to modify the given script `run_md.sh`. For example you compiled the gromacs `gmx` command to have the postfix `_mpi` such that the binary compiled is `gmx_mpi` you can substitute the `gmx` command in the script using sed inplace (`sed -i ...`) like so:

```bash
$ sed -i "s/ gmx / gmx_mpi /g" run_md.sh
```

Additionally, it is required to make the shell script executable with the following command:
```bash
$ chmod +x run_md.sh
```

### When using the flags provided in the helper script as described below you will be able to:
- Generate extended structure
- Run a short vacuum simulation and extracting $N_{replica}$ frames as starting structures
- Solvate and add ions for neutralization.
- Energy minimize each system  simultaneously by running them in parallel.
- Thermalize each system to 300 K using _NVT_ ensemble using V-rescale thermostat.
- Equilibrate the box size in two phases using _NPT_ ensemble : 
  - **Phase 1 :** Using Berendsen barostat, which helps in attaining equilibrium pressure quickly.
  - **Phase 2 :** Using Parrinello-Rahman barostat, which helps in generating biologically relavant ensembles.  
- Run REST2 Simulations from the final snapshot of each replica. 

> [!NOTE] 
>  With the implementation of C-rescale barostat in gromacs, we can use it during both equilibration and production.

```bash
$ run_md.sh -h or --help   # provides all options 

Usage: usage [Options]
Options:
 -h, --help     Display this help message
 -v, --verbose  Enable verbosity
 -g, --generate Generate extended structure with pmx and extract 10 frames from a short vacuum simulation
 -p, --topol  Topology File name for rest
 -s, --stage   Select current phase of simulation: Setup, Min, NVT, NPT, REST
 -l, --log      STDIO log File
 -strip, --strip Stripped TPR file excluding waters
```
If you desire to see verbose output, commands printed to screen and verbose output from gromacs use the `-v or --verbose` flags. 
Provided is a flag to log all output from the script, `-l <logname> or --log=<logname>`.

## Running Simulations
<!-- The preparation and simulations are separated into stages. 
Stages need to be performed in order: -->
The preparation and execution of REST2 are performed in the following order : 
* Generate extended chain from fasta sequence and extract N (=10) frames, one for each replica, from a short vacuum simulation. 
* Setup each replica by solvating and neutralizing the protein box.
* Energy minimization.
* Equilibration containing the following stages:
  *  Thermalization and _NVT_ equilibration.
  *  _NPT_ equilibration with the Berendsen barostat.
  *  _NPT_ equilibration with the Parrinello-Rahman barostat.
*  Replica Exchange with Solute Scaling (REST2) molecular dynamics simulations.

> [!NOTE] 
>  C-rescale barostat can be used during _NPT_ equilibration and production.

### Generating starting structures
For this tutorial we are simulating a 20-residue protein fragment from $\alpha$-synuclien, specifically the last 20 residues of the C-terminal. The 20-residue sequence in question will be referred to ad $\alpha$-synuclien for the remainder of the tutorial. The FASTA sequence for the terminal 20-residues of the C-terminal domain is:
```bash
DMPVDPDNEAYEMPSEEGYQDYEPEA
```

To produce the 10 initial starting structures, enter the following command and take note of the residue sequence appended:
```bash
$ run_md.sh -g DMPVDPDNEAYEMPSEEGYQDYEPEA
$ run_md.sh --generate DMPVDPDNEAYEMPSEEGYQDYEPEA
```
>[!NOTE]
> Both command entries are equivalent, either will do and you only need to select one to generate starting structures.
<!-- Both command entries are equivalent, either will do and you only need to select one to generate starting structures. -->

### Simulation setup
The helper script will setup the system with protein, water and ions with this command:
```bash
$ run_md.sh -s setup
```
This will produce the system topology and system.gro. 

### Energy minimization
Energy minimization on the well solvated structures can be performed by:
```bash
# Minimize the system using system.gro and topol.top produced in the setup stage
$ run_md.sh -s Min
```
### Equilibration
Equilibration stages:
* NVT
  ```bash
  # Thermalize and equilibrate system under the _NVT_ ensemble using min.gro and topol.top
  $ run_md.sh -s NVT
  ```
* NPT
  ```bash 
  # Equilibrate the pressure and volume of the system using nvt.gro and topol.top
  run_md.sh -s NPT
  ```
  
<!-- These stages can be performed with these series of commands: -->
Test for pressure and volume convergence with `gmx energy`.

### Replica Exchange Simulations

#### Creating the Replica Exchange base topology which will be used for scaling

Once well equilibrated, creating a topology file containing the respective scaled information is required to setup the rest directories and generate the structure input files, `prod.tpr`. This edited topology will contain subscripts appended to the atom types of the molecule of interest, the protein. The `-s NPT` stage will produce the file `processed.top`. Before proceeding with the REST2 simulations modify the respective atom types with the command below. `start_line` and `end_line` are the line numbers before and after the lines containing the atom type information of proteins in `[atoms]` section, respectively. For example, open the file using vim text editor.
```bash
$ vim processed.top
```
Once `processed.top` is open in vim type the following command including the semicolon to enable the line numbers.
```vim
:set nu
```
`:set nu!` will disable it. The file near the `[ atoms ]` sections looks like :
```bash
$ vim processed.top
.
.
.
964 [ atoms ]
965 ;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
966 ; residue 121 NASP rtp NASP q  0.0 <---START LINE
967      1         N3    121   NASP      N      1     0.0782      14.01
968      2          H    121   NASP     H1      2       0.22      1.008
969      3          H    121   NASP     H2      3       0.22      1.008
970      4          H    121   NASP     H3      4       0.22      1.008
971      5         CT    121   NASP     CA      5     0.0292      12.01
.        .
.        .
.        .
1276    290         HC    140   CALA    HB2    290     0.0764      1.008
1277    291         HC    140   CALA    HB3    291     0.0764      1.008
1278    292          C    140   CALA      C    292     0.7731      12.01
1279    293         O2    140   CALA    OC1    293    -0.8055         16
1280    294         O2    140   CALA    OC2    294    -0.8055         16   ; qtot -7.573
1281 <---END LINE
1282 [ bonds ]
.
.
.

```

Identify the start and endline and run the following bash command replacing the words start and end with the appropriate line numbers.
```bash
awk -v start_line=start -v end_line=end '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top

```
See the [REST.top](./reference_files/og_topology/REST.top) as a reference.

>[!NOTE]
>One can use any text editor or one is encouraged to write a small code to find the line numbers and add the subscripts to atom types. For us using vim made this process easier.

<!-- `start_line` and `end_line` are defined as the line numbers before and after the proteins `[atoms]` section, respectively -->



#### Running replica exchange simulations
You are strongly encouraged to have at least 2 GPUS; 5 or 10 GPUS would be best. If you are running **CPU only** simulations you are encouraged to have a minimum of 100 CPU cores, 10 cores per simulation, expect this to take a long time a month or two. If you are running **GPU** simulations you are encouraged to have a minimum of 2 cpus per replica, 20 CPUS total, 8 or 16 CPUS per simulation is better. With 2 GPUs the simulation will take on the order of 10 days, while using 10 GPUS will take approximately 3 days (tested on NVIDIA A5500 GPUs).

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
  

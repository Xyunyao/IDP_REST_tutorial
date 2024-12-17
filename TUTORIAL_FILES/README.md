# System setup

Before proceeding to running simulations, please make sure to install the required software as mentioned in [INSTALLATION_INSTRUCTIONS](./INSTALLATION_INSTRUCTIONS/). 
It is important to note, GROMACS must be compiled with mpi, patched with plumed and the executable must be:
```bash
gmx
```
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
```
Stages need to be performed in order. Once well equilibrated, creating a rest2 edited topology file is required. This edited topology will contain subscripts appended to the atom types of the molecule of interest, the protein. The `-s NPT` stage will produce the file `processed.top`, copy this file to a unique name, e.g. REST.top, before proceeding with the REST2 simulations. Modify the atom types with the command below. `start_line` and `end_line` are the line numbers where the `[atoms]` section start and end respectively. See the [REST.top](./reference_files/og_topology/REST.top) as a reference.

```bash
awk '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top
```

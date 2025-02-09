# Analysis 

Hola! You are here means REST2 simulations are successfully simulated. Let's learn how or what kind of analysis one can perform on the simulated trajectories. Pre-simulated trajectories can be downloaded from this [here](https://doi.org/10.5281/zenodo.14799045).


## Demultiplexing of trajectories.

First step in the analysis workflow is to extract the absolute replica path followed by respective replicas. This will be helpful as a sanity check and debugging of our REST2 simulations. `scripts/demux.fix.pl` is used to generate a set of data files which contains the required information useful for demultiplexing. We use the `prod.log` file of the base replica as an input. When prompted input the time step value of the simulation which is 0.002 ps in this case.

```bash
perl demux.fix.pl <path to base replica prod.log>
```

This generates `replica_temp.xvg` and `replica_index.xvg` as an output. But since the frequency of saving co-ordinates in `prod.mdp` is different from the proposed frequency of exchange attempts this need to be corrected. The correction is performed using:

```bash
awk '{if($1==0){print} if(n==50){$1=$1-80.0; print;n=0} n++;}' replica_index.xvg > replica_index.n50.s0.-80.xvg
```
- n50 = every 50th frame
- s0 = start at frame 0
- -80 shift time index by -80

If `.trr` is used where the frame deposition frequency is for every 800.0ps, we shift the time index by -800.0ps,

```bash
awk '{if($1==0){print} if(n==500){$1=$1-800.0; print;n=0} n++;}' replica_index.xvg > replica_index.n500.s0.-800.xvg
```
- n500 = every 500th frame
- s0 = start at frame 0
- -800 shift time index by -800

Now index file corrected, it is used to generate the demultiplexed trajectories. This can be done using `scripts/make_demux.sh` as :

```bash
% ./make_demux.sh -h
Usage: ./make_demux.sh [Options]
Options:
 -h, --help       Display this help message
 -v, --verbose    Enable verbosity
 -i, --index      replica index xvg file
 -d, --restdir    replica exchange directory
 -n, --name       trajectory name, i.e. whole.xtc
 -l, --log        STDIO log File
```

```bash
./make_demux.sh -i <number of replicas> -i <path to replica index xvg file> -d <path to replica directories> -n <[whole] prefix name of replica xtc file>
```

Now trajectories containing the absolute path taken by the respective replicas are generated.

## Periodic Boundary corrections

Before proceeding with analysis correction of the trajectories for the periodic boundary conditions which was implement during the simulations. To do that a `.tpr` file with only protein and ions information is required. This can be done using the following command :

```bash
gmx convert-tpr -s <production tpr> -o <protein-ions tpr> -n index.ndx
```

`index.ndx` file is optional. When prompted select `non-Water` for selection.

Gromacs can be used to perform the corrections. If wanted there are many other wrapper programs which can be used. Following the below commands one can perform PBC corrections.
- To make the molecules whole :
    ```bash
    gmx trjconv -f <trajectory file(.xtc)> -s <structure file(.tpr)> -o whole.xtc -pbc whole
    ```
- To correct for the jumps : (optional)
    ```bash
    gmx trjconv -f whole.xtc -s <structure file(.tpr)> -o nojump.xtc -pbc nojump
    ```
- To center the protein : (optional)
    ```bash
    gmx trjconv -f nojump.xtc -s <structure file(.tpr)> -o center.xtc -pbc center
    ```
- To correct the periodic boundaries :
    ```bash
    gmx trjconv -f <whole or nojump or center xtc> -s <structure file(.tpr)> -o pbc.xtc -pbc mol -ur compact
    ```

This should help in correcting the periodic boundary corrections to the trajectories which we will be analyzing now.

## Preliminary analysis

Now the trajectories are ready to compute the macroscopic parameters. In the notebook `analysis.ipynb` we will compute the the parameters that we discussed in the chapter such as Radius of gyration, inter-residue contact map etc. `replica_temp.xvg` and `replica_index.xvg` are used to accessing the convergence of the replica exchange simulation.  

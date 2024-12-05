# Analysis 

Hola! You are here means you successfully simulated an set of trajectories using REST2 method. Let's learn how or what kind of analysis we can perform on the simulated trajectories.

## Demultiplexing of trajectories.

First step of our analysis workflow is to extract the absolute replica path followed by respective replicas. This will be helpful as a sanity check and debugging of our REST2 simulations. `demux.fix.pl` is used to to generate a set of data files which contains the required information useful for demultiplexing. We use the [prod.log](https://dartmouth-my.sharepoint.com/:f:/r/personal/f006f50_dartmouth_edu/Documents/trajectories_for_book_chapter?csf=1&web=1&e=h22nwg) file of the base replica as an input. When prompted we need to mention the time step of our simulations which is 0.002 ps in our case.

```bash
perl demux.fix.pl prod.log
```

This gives us `replica_temp.xvg` and `replica_index.xvg` as an output. But since the frequency of saving co-ordinates in `prod.mdp` is different from the proposed frequency of exchange attempts we need to correct this. The correction is performed using:

```bash
awk '{if($1==0){print} if(n==50){$1=$1-80.0; print;n=0} n++;}' replica_index.xvg > replica_index.n50.s0.-80.xvg
```
- n==50 = every 50th frame
- s0 = start at frame 0
- -80 shift time index by -80

If we are taking the replica index every 800.0ps, we shift the time index by -800.0ps,

```bash
awk '{if($1==0){print} if(n==500){$1=$1-800.0; print;n=0} n++;}' replica_index.xvg > replica_index.n500.s0.-800.xvg
```
- n500 = every 500th frame
- s0 = start at frame 0
- -800 shift time index by -800

Now we have a corrected index file, we use it generate the demultiplexed trajectories. This can be done using `make_demux.sh` as :

```bash
sh make_demux.sh <path/to/replica/index/xvg/file/> <path/to/replica/directories/>
```

Now we have a set of trajectories which has the absolute path taken by the respective replicas.

## Periodic Boundary corrections

Before we proceed with the analysis we need to correct the trajectories for the periodic boundary conditions which we implement during the simulations. To do that let's create a `.tpr` file with only protein and ions information. This can be done using the following command :

```bash
gmx convert-tpr -s <production tpr> -o <protein-ions tpr> -n index.ndx
```

`index.ndx` file is optional. When prompted select `non-Water` for selection.

Gromacs can be used to perform the corrections. If you want you can use other wrapper programs. Following the below commands one can perform PBC corrections.
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


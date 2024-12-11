# System setup

Before proceeding to start setting up the simulations please make sure to install required softwares as mentioned in [INSTALLATION_INSTRUCTIONS](./INSTALLATION_INSTRUCTIONS/). Once installed then yay let's start your first enhanced sampling simulation. User will be provided with required files and a script `run_md.sh` which will be used to setup and run simulations step by step.

```bash
run_md.sh explaination and step by step options
```

Once well equilibrated solvated structures were generated, one small edit needs to be done to `processed.top` before proceeding with the REST2 simulations. `_` needs to be added to the atoms of interest which needs to be scaled. Use the below command to tag the atoms. `start_line` and `end_line` are the line numbers where the `[atoms]` section start and end respectively. See the [REST.top](./reference_files/og_topology/REST.top) as a reference.

```bash
awk '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top
```

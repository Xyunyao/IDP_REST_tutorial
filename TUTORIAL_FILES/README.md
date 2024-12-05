# System setup

Before we proceed to star setting up our simulations make sure to install required softwares as mentioned in [INSTALLATION_INSTRUCTIONS](./INSTALLATION_INSTRUCTIONS/). If you already installed them then yay you are ready to start your first enhanced sampling simulation. We provided user with required files and a script `run_md.sh` which we will use to setup and run simulations step by step.

```bash
run_md.sh explaination and step by step options
```

Once we have well equilibrated solvated structures, we need to make a small edit to the `processed.top`.  We need to add `_` to the atoms of interest which will be scaled. Use the below command to tag the atoms. `start_line` and `end_line` are the line numbers where the `[atoms]` section start and end respectively. See the [REST.top](./reference_files/og_topology/REST.top) as a reference.

```bash
awk '{n++; if(n>start_line) if(n<end_line) if(NF>7) if($1+0==$1){$2=$2"_"}; print;}' processed.top > REST.top
```


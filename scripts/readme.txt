This is an updated instruction to run REST2 with Slurm based on original run_md.sh file
Main difference from run_md.file:

1. setup for run in HPC with slurm scheduler
2. break run_md.sh into several sections: 

    The orignal run_md.sh file can be run with -s stage flag

    (-s: gen, min, nvt, npt, rest)

    I divide the whole script into 4 scripts:

    * change the .sh file permission to executable by (chmod +x name_of_script.sh)

    1. script to generate the initial structre, i.e. extended.pdb (optional)
    There are multiways to generate this (discussed in  paper:https://arxiv.org/pdf/2505.01860): 1. alphafold model 2. same as the tutorial using generate_chain.py
    for 2, you need to install pmx (not inclued in the setup_env.yml) 
    (I created a new conda environment to install it due to some verions conflict to softeare in REST_tutorial environment)
   
    2. script to run gen min nvt npt as parallel jobs:
    You can run in -s stage flag model (good for troubleshooting)
    But I recommend running without -s (run all stages all together) since the scipt copy the job to scrapt folder to run the jobs.
    Run them all together avoids repeating copying documents
    A little computationally expensive due to npt1 run (40 ns for each replica)
    run as: sbatch  

    3. script to prepare rest2 runnign files
    It generates the replica folds and necessay files for run REST2.
    The most important is generating the scaled topology files. 
    * The lamda_value function is wrong, it calculate the Tm temperature, not the scaling factor
    * It automatically find the first [atom] sections in processed.pp file and mark the hotatoms (This part is missing in the original run_md.sh script)

    It is computationally light. So run as simple as:
    ./   .sh  

    4. script to run actual rest2 production
    warning: compressed very expensive (1us for each replica, change num of steps in the mdp file in prod.mdp for study purpose)
    using slurm scheduler (4 gpu + 2*10=20 cpu to run the job; modify it accordingly based on your system)
    run as: sbatch ./    .sh


If you want to know more about the change in the scaled topology file, read scaled_topology.txt file


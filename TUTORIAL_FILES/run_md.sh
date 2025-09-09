#!/bin/bash

# Created by Korey Reid
# updated by Yunyao for running with Slurm using GPU
# Equilibrium simulation for use with mpi (e.g. mpich, openmpi, etc.)
# Requires mpirun to be in the PATH shell global variable
# If you want to use GPUs you *must* compile gromacs with GPU support (e.g. CUDA) 
# The code provided has flags to perform *Equilibration* and REST2 Simulations.
# 
# MIN -> NVT -> NPT -> REST2
# Ask yourself a simple question regarding equilibration, are things converged?
# If you have not already done so, parse over the Lysozyme tutorial for the basics.
#
# If performing on a cluster, it is advised to purge the environment and only load 
# necessary modules
# For example:
# % module purge
# % module load openmpi cuda
# % module list
#   Currently Loaded Modulefiles:
#       1) cuda/12           2) openmpi/4.1.3
#
# It is entirely possible that an academic HPC environment has gromacs with plumed,
# I leave this up to you to make sure/edit this script to work in that environment. 
# It is entirely to difficult to predict the variety in cluster environments when 
# writing this script and it is assumed gromacs is compiled with our gromacs/plumed 
# installation script provided with this tutorial on the computer/compute node the simulations are performed. 
# 

# Default variables

verbose=""
log=""
stage=""
topol=""
tpr=""
fasta=""
equil_dir="Equilibration"
initstructs="initial_structures"
tutdir=`pwd`
ntomp=4  # number of threads per task

# =====================
# Usage function
# =====================
usage() {
    echo "Usage: $0 [Options]"
    echo "Options:"
    echo " -h, --help       Display this help message"
    echo " -v, --verbose    Enable verbosity"
    echo " -g, --generate   Generate extended structure with pmx and extract frames from a short vacuum simulation"
    echo " -p, --topol      Topology File name for REST"
    echo " -s, --stage      Select simulation stage: Gen, Setup, Min, NVT, NPT, REST, PBC"
    echo " -l, --log        STDIO log File"
    echo " -strip, --strip  Stripped TPR file excluding waters"
}

has_argument() { [[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]]; }
get_argument() { echo "${2:-${1#*=}}"; }

lambda_value () {
    Ti=300
    Tf=450
    awk -v tmin=$Ti -v tmax=$Tf -v n=10 -v i=$@ 'BEGIN{ t=tmin*exp(i*log(tmax/tmin)/(n)); printf(t); }'
}

# =====================
# Flag handler
# =====================
handle_flags() {
    while [ $# -gt 0 ]; do
        case $1 in
            -h|--help) usage; exit 0 ;;
            -v|--verbose) verbose=" -v" ;;
            -g|--generate) fasta=$(get_argument $@); shift ;;
            -p|--topol*) topol=$(get_argument $@); shift ;;
            -strip|--strip*) tpr=$(get_argument $@); shift ;;
            -s|--stage*) stage=$(get_argument $@); shift ;;
            -l|--log*) log=$(get_argument $@); shift ;;
            *) echo "Invalid option:$1"; usage; exit 1 ;;
        esac
        shift
    done
}

# =====================
# Function to run mdrun safely
# =====================
run_mdrun() {
    local deffnm=$1
    local tprfile=$2
    if [ ! -f "$tprfile" ]; then
        echo "Error: $tprfile not found. Generate it first with grompp."
        exit 1
    fi
    gmx_mpi mdrun ${verbose} -deffnm "$deffnm" -nb gpu -ntomp $ntomp
}

# =====================
# Main simulation function
# =====================
simulations () {
    case ${stage,,} in
        gen)
            eval "$(conda shell.bash hook)"
            conda activate REST_tutorial
            mkdir -p ${initstructs}
            cd ${initstructs}

            printf "15\n1\n" | gmx_mpi pdb2gmx -f ${tutdir}/extended.pdb -o protein.gro
            gmx_mpi editconf -f protein.gro -o vac_box.gro -box 15 15 15
            cp ${tutdir}/mdp_files/minimz.mdp .
            cp ${tutdir}/mdp_files/VAC.mdp .

            gmx_mpi grompp -f minimz.mdp -c vac_box.gro -p topol.top -o em.tpr -maxwarn 1
            run_mdrun em em.tpr

            gmx_mpi grompp -f VAC.mdp -c em.gro -p topol.top -o vac.tpr -maxwarn 1
            run_mdrun vac vac.tpr

            printf "1\n" | gmx_mpi trjconv -s vac.tpr -f vac.xtc -o starting_structure.pdb -b 20 -dt 20 -sep

            cd ${tutdir}
            echo "Starting structures prepared under ${tutdir}/${initstructs}"
            ;;

        setup)
            box=6.5
            mkdir -p ${equil_dir} && cd ${equil_dir}
            for i in {0..9}; do
                mkdir -p ${i}
                cp ../${initstructs}/starting_structure${i}.pdb ${i}/
                cp -r ${tutdir}/initial_input_files/a99SBdisp.ff ${i}/
                cp ${tutdir}/mdp_files/{NPT0.mdp,NPT1.mdp,NVT.mdp,minimz.mdp} ${i}/
            done

            for i in {0..9}; do
                cd ${i}
                printf "1\n1\n" | gmx_mpi pdb2gmx -f starting_structure${i}.pdb -o prot.gro 
                gmx_mpi editconf -f prot.gro -o box.gro -bt cubic -box ${box} ${box} ${box} -c
                gmx_mpi solvate -cp box.gro -cs a99SBdisp.ff/a99SBdisp_water.gro -p topol.top -o solvate.gro
                gmx_mpi grompp -f ../minimz.mdp -c solvate.gro -p topol.top -o ions.tpr
                printf "13\n" | gmx_mpi genion -s ions.tpr -pname NA -nname CL -neutral -p topol.top -o ions.gro 
                cd ..
            done
            cd ${tutdir}
            echo "Systems prepared."
            ;;

        min|nvt|npt)
            cd ${tutdir}/${equil_dir}
            # Determine mdp and prefix
            case ${stage,,} in
                min) mdp_file=minimz.mdp; prefix=min ;;
                nvt) mdp_file=NVT.mdp; prefix=nvt ;;
                npt) mdp_file=NPT0.mdp; prefix=npt0 ;;
            esac

            for i in {0..9}; do
                cd ${i}
                gmx_mpi grompp -f ${tutdir}/mdp_files/${mdp_file} -c system.gro -p topol.top -o ${prefix}.tpr -maxwarn 2
                cd ..
            done

            for i in {0..9}; do
                cd ${i}
                run_mdrun ${prefix} ${prefix}.tpr
                cd ..
            done
            cd ${tutdir}
            ;;

        rest)
            if [ ! -f $topol ]; then echo "REST topology ($topol) does not exist."; usage; exit 1; fi
            cd ${tutdir}/${equil_dir}
            mkdir -p REST
            cp ${tutdir}/mdp_files/prod.mdp REST/
            cd REST
            for i in {0..9}; do
                mkdir -p ${i}
                plumed partial_tempering $(lambda_value ${i}) < ${tutdir}/${topol} > ${i}/topol.top
                touch ${i}/plumed.dat
                cd ${i}
                gmx_mpi grompp -f ../prod.mdp -c ${tutdir}/${equil_dir}/${i}/npt1.gro -p topol.top -o prod.tpr -maxwarn 2
                run_mdrun prod prod.tpr
                cd ..
            done
            cd ${tutdir}
            ;;

        pbc)
            if [ -z "${tpr}" ] || [ ! -f "${tpr}" ]; then echo "Must provide a valid stripped TPR file."; usage; exit 1; fi
            for i in {0..9}; do
                cd ${i}
                printf "1\n0\n" | gmx_mpi trjconv -s ${tpr} -f prod.xtc -o whole.xtc -pbc mol -center
                printf "1\n0\n" | gmx_mpi trjconv -s ${tpr} -f prod.gro -o prod.pdb -pbc mol -center
                cd ..
            done
            ;;

        *)
            echo "stage name must match one of: Gen, Setup, Min, NVT, NPT, REST, PBC"
            usage
            exit 1
            ;;
    esac
}

# =====================
# Main
# =====================
handle_flags "$@"

if [ "${verbose}" = true ]; then
    set -x
    set -v
fi

if [ -n "${log}" ]; then
    exec > >(tee ${log}) 2>&1
fi

simulations
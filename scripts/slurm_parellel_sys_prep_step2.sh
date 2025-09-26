#!/bin/bash
#SBATCH --job-name=protein_md_array
#SBATCH --nodes=1
#SBATCH --ntasks=1              # MPI ranks per replica
#SBATCH --cpus-per-task=8       # threads per replica
#SBATCH -p fgpu
#SBATCH --gres=gpu:a4000:1        # 1 GPU shared by all tasks
#SBATCH --time=5:30:00         # walltime with buffer
#SBATCH --array=0-9           # 7 replicas
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

# -------------------------------
# Variables use npt1 as example
# -------------------------------
WORKDIR=$HOME/IDP_REST_tutorial/TUTORIAL_FILES
SCRATCH=/local/$USER/job_$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
MD_TPR=npt1.tpr
MD_CPT=npt1.cpt    # checkpoint for npt1 continuation
OUTPUT_DIR=$WORKDIR/Equilibration2/$SLURM_ARRAY_TASK_ID
#JOB_NAME=npt1
#equil_dir="Equilibration"
initstructs="initial_structures"
export ntomp=8
# -------------------------------
# Create scratch directory for this replica
# -------------------------------
mkdir -p $SCRATCH
cd $SCRATCH

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
# 
# copy starting structure and mdp files
echo "Copying input files for replica $SLURM_ARRAY_TASK_ID to scratch..."
cp $WORKDIR/initial_structures/starting_structure${SLURM_ARRAY_TASK_ID}.pdb .
cp $WORKDIR/mdp_files/*  .
cp -r $WORKDIR/initial_input_files/a99SBdisp.ff  ./  # copy force field files

# -------------------------------
# Function to copy results back
# -------------------------------
copy_back() {
    echo "Copying results of replica $SLURM_ARRAY_TASK_ID back to $OUTPUT_DIR..."
    mkdir -p $OUTPUT_DIR  # create output directory if it doesn't exist
    cp -r $SCRATCH/* $OUTPUT_DIR/
}

# -------------------------------
# Trap SIGTERM (SLURM time limit)
# -------------------------------
trap 'echo "Replica $SLURM_ARRAY_TASK_ID terminated due to time limit"; copy_back; exit 0' SIGTERM

# -------------------------------
# Load cuda module
# -------------------------------
module load cuda12.2

# -------------------------------
# Run GROMACS
# -------------------------------
prepare_simulations () {
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
            run_mdrun em em.tpr ""

            gmx_mpi grompp -f VAC.mdp -c em.gro -p topol.top -o vac.tpr -maxwarn 1
            run_mdrun vac vac.tpr ""

            printf "1\n" | gmx_mpi trjconv -s vac.tpr -f vac.xtc -o starting_structure.pdb -b 20 -dt 20 -sep

            cd ${tutdir}
            echo "Starting structures prepared under ${tutdir}/${initstructs}"
            ;;

        setup)
            box=6.5
        
            cd $SCRATCH
            printf "1\n1\n" | gmx_mpi pdb2gmx -f starting_structure${SLURM_ARRAY_TASK_ID}.pdb -o prot.gro
            gmx_mpi editconf -f prot.gro -o box.gro -bt cubic -box ${box} ${box} ${box} -c
            # add the same number of waters to each system (the minium of all replicas, need test run to find)
            gmx_mpi solvate -cp box.gro -cs a99SBdisp.ff/a99SBdisp_water.gro -p topol.top -o solvate.gro -maxsol 8890
            gmx_mpi grompp -f minimz.mdp -c solvate.gro -p topol.top -o ions.tpr -maxwarn 2
            printf "13\n" | gmx_mpi genion -s ions.tpr -pname NA -nname CL -neutral -p topol.top -o ions.gro 
            echo "Systems prepared. check ion.gro to confirm the number of waters and ions."
            ;;

        min)
            cd $SCRATCH
            if [ ! -f "ions.gro" ]; then
                echo "ions.gro missing from replica $SLURM_ARRAY_TASK_ID, run the setup stage"
                exit 1
            fi
            gmx_mpi grompp -f minimz.mdp-c ions.gro -p topol.top -o min.tpr -maxwarn 2
            echo "Performing System Minimization."
            gmx_mpi mdrun -s min.tpr -deffnm min -nb gpu -ntomp $ntomp -dlb no | tee min.log
            echo "Minimized system."
            ;;

        nvt)
            cd $SCRATCH
            if [ ! -f "min.gro" ]; then
                echo "min.gro missing from replica $SLURM_ARRAY_TASK_ID, run the min stage"
                exit 1
            fi
            gmx_mpi grompp -f NVT.mdp -p topol.top -c min.gro -o nvt.tpr -maxwarn 2
            echo "Thermalizing system under NVT ensemble."
            srun gmx_mpi mdrun -s nvt.tpr -deffnm nvt -nb gpu -ntomp $ntomp -dlb no | tee nvt.log
            echo "System NVT equilibration complete."
            ;;

        npt)
            cd $SCRATCH
            if [ ! -f "nvt.gro" ]; then
                echo "nvt.gro missing from replica $SLURM_ARRAY_TASK_ID, run the nvt stage"
                exit 1
            fi
            gmx_mpi grompp -f NPT0.mdp -p topol.top -c nvt.gro -o npt0.tpr -maxwarn 2
            echo "System NPT0 equilibration."
            srun gmx_mpi mdrun -s npt0.tpr -deffnm npt0 -nb gpu -ntomp $ntomp -dlb no | tee npt0.log
            echo "NPT0 equilibration complete."

            # NPT1 stage
            gmx_mpi grompp -f NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr -maxwarn 2
            echo "System NPT1 equilibration."

            if [ -f "$MD_CPT" ]; then
                echo "Replica $SLURM_ARRAY_TASK_ID: Restarting from checkpoint..."
                srun gmx_mpi mdrun -s $MD_TPR -cpi $MD_CPT -deffnm npt1 \
                    -ntomp $SLURM_CPUS_PER_TASK -nb gpu -dlb no -noappend | tee $JOB_NAME.log
            else
                echo "Replica $SLURM_ARRAY_TASK_ID: Starting new run..."
                srun gmx_mpi mdrun -s $MD_TPR -deffnm npt1 \
                    -ntomp $SLURM_CPUS_PER_TASK -nb gpu -dlb no | tee $JOB_NAME.log
            fi
            ;;
        *)
            echo "Unknown stage: $stage"
            exit 1
            ;;
    esac
}
# should run after checking number of waters in each system in ions.gro
prepare_simulations_one_go () {
    cd $SCRATCH
    cp $OUTPUT_DIR/ions.gro .
    cp $OUTPUT_DIR/topol.top .
    if [ ! -f "ions.gro" ]; then
        echo "ions.gro missing from replica $SLURM_ARRAY_TASK_ID, run the setup stage"
        exit 1
    fi
    gmx_mpi grompp -f minimz.mdp -c ions.gro -p topol.top -o min.tpr -maxwarn 2
    echo "Performing System Minimization."
    gmx_mpi mdrun -s min.tpr -deffnm min -nb gpu -ntomp $ntomp -dlb no | tee min.log
    echo "Minimized system."

    gmx_mpi grompp -f NVT.mdp -p topol.top -c min.gro -o nvt.tpr -maxwarn 2
    echo "Thermalizing system under NVT ensemble."
    srun gmx_mpi mdrun -s nvt.tpr -deffnm nvt -nb gpu -ntomp $ntomp -dlb no | tee nvt.log
    echo "System NVT equilibration complete."

    gmx_mpi grompp -f NPT0.mdp -p topol.top -c nvt.gro -o npt0.tpr -maxwarn 2
    echo "System NPT0 equilibration."
    srun gmx_mpi mdrun -s npt0.tpr -deffnm npt0 -nb gpu -ntomp $ntomp -dlb no | tee npt0.log
    echo "NPT0 equilibration complete."

    # NPT1 stage
    gmx_mpi grompp -f NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr -maxwarn 2
    echo "System NPT1 equilibration."

    if [ -f "$MD_CPT" ]; then
        echo "Replica $SLURM_ARRAY_TASK_ID: Restarting from checkpoint..."
        srun gmx_mpi mdrun -s $MD_TPR -cpi $MD_CPT -deffnm npt1 \
            -ntomp $SLURM_CPUS_PER_TASK -nb gpu -dlb no -noappend | tee npt1.log
    else
        echo "Replica $SLURM_ARRAY_TASK_ID: Starting new run..."
        srun gmx_mpi mdrun -s $MD_TPR -deffnm npt1 \
            -ntomp $SLURM_CPUS_PER_TASK -nb gpu -dlb no | tee npt1.log
    fi
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

# choose run stage
if [ -z "${stage}" ]; then
    echo "Run whole simulation in one go."
    prepare_simulations_one_go
else
    prepare_simulations
fi

# -------------------------------
# Copy results back at the end
# -------------------------------
copy_back

# -------------------------------
# Clean up scratch
# -------------------------------
rm -rf $SCRATCH


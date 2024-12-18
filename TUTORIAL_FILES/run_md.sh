#!/bin/bash

# Created by Korey Reid
#
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


# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
    echo " -p, --topol  Topology File name for rest"
    echo " -s, --stage   Select current phase of simulation: Setup, Min, NVT, NPT, REST, PBC"
	echo " -l, --log	STDIO log File"
    echo " -strip, --strip Stripped TPR file excluding waters"
}

has_argument() {
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]] ;
}

get_argument() {
	echo "${2:-${1#*=}}"
}

lambda_value () {
    Ti=300
    Tf=450
    awk -v tmin=$Ti \
    -v tmax=$Tf \
    -v n=10 \
    -v i=$@\
  'BEGIN{
    t=tmin*exp(i*log(tmax/tmin)/(n));
    printf(t);
}'
}

# Flag handler
handle_flags() {
	while [ $# -gt 0 ]
	do
		case $1 in 
			-h | --help)
				usage
				exit 0
				;;
			-v | --verbose)
				verbose=" -v"
				;;
			-p | --topol*)
				if ! has_argument $@
				then
					echo "File not specified." >&2
					usage
					exit 1
				fi

				topol=$(extract_argument $@)

				shift
				;;
            -strip | --strip*)
				if ! has_argument $@
				then
					echo "Stripped TPR File not specified." >&2
					usage
					exit 1
				fi

				tpr=$(extract_argument $@)

				shift
				;;
            -s | --stage*)
                if ! has_argument $@
                then
                    echo "When using the stage flag you must provide Min, NVT, NPT, or REST"
                    usage
                    exit 1
                fi

                stage=$(extract_argument $@)

                shift
                ;;
            -l | --log*)
                if ! has_argument $@
                then
                    echo "When using the log flag you must provide a file name"
                    usage
                    exit 1
                fi

                log=$(extract_argument $@)

                shift
                ;;
			*)
				echo "Invalid option:$1" >&2
				usage
				exit 1
				;;
		esac
		shift
	done
}

simulations () {
    case ${stage,,} in
        setup)
            box=6.5
            echo "Solvating with box size $box nm."
            printf "0\n1\n" | gmx pdb2gmx -f initial_input_files/prot_only.pdb -o prot.gro 
            gmx editconf -f prot.gro -o box.gro -bt cubic -box $box $box $box -c
            gmx solvate -cp box.gro -cs a99SBdisp.ff/a99SBdisp_water.gro -p topol.top -o solvate.gro
            gmx grompp -f mdp_files/minimz.mdp -c solvate.gro -p topol.top -o ions.tpr
            printf "13\n" | gmx genion -s ions.tpr -pname NA -nname CL -neutral -p topol.top -o system.gro 

            echo "System Prepared."
            ;;
        min)
            if [ ! -f "ions.gro" ]
            then
                echo "ions.gro missing, run the setup stage"
                exit 1
            fi
            echo "Performing System Minimization."
            gmx grompp -f mdp_files/minimz.mdp -c system.gro -p topol.top -o min.tpr -maxwarn 2
            gmx mdrun $verbose -s min.tpr -deffnm min
            echo "Minimized system."
            ;;
        nvt)
            if [ ! -f "min.gro" ]
            then
                echo "min.gro missing, run the min stage"
                exit 1
            fi
            echo "Thermalizing the system and Equilibrating to 300K under the NVT Ensemble."
            gmx grompp -f mdp_files/NVT.mdp -p topol.top -c min.gro -o nvt.tpr 
            gmx mdrun $verbose -s nvt.tpr -deffnm nvt
            echo "System equilibrated under NVT."
            ;;
        npt)
            if [ ! -f "nvt.gro" ]
            then
                echo "nvt.gro missing, run the nvt stage"
                exit 1
            fi
            echo "Pressure and Temperature Equilibration at 1 bar and 300K (NPT)"
            echo "stage 1: Short Berendsen pressure relaxation."
            gmx grompp -f mdp_files/NPT0.mdp -p topol.top -c nvt.gro -o npt0.tpr 
            gmx mdrun $verbose -s npt0.tpr -deffnm npt0
            echo "stage 2: Long Parrinello-Rahman pressure Equilibration."
            gmx grompp -f mdp_files/NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr -pp processed.top
            gmx mdrun $verbose -s npt1.tpr -deffnm npt1
            echo "System NPT equilibration complete, check volume and pressure convergence."
            echo "Before performing REST2 Simulations, copy the processed.top file and"
            echo "add underscores, e.g. HA_, follow the README guide." 
            ;;
        rest)
            if [ ! -f $topol ]
            then 
                echo "Rest topology does not exist."
                usage
                exit 1
            fi
            if [ ! -f "npt1.gro" ]
            then
                echo "npt1.gro missing, run the npt stage"
                exit 1
            fi

            echo "Creating rest directory and subdiretories."
            mkdir REST
            cdir=`pwd`

            cp mdp_files/prod.mdp REST/
            cd REST
            echo "Scaling topologies for each replica and creating run input file prod.tpr"
            for((i=0;i<10;i++))
            do
                mkdir $i
                plumed partial_tempering $(lambda_value $i) < $topol > $i/topol.top
                touch plumed.dat
                
                cd $i/
                gmx grompp -f ../prod.mdp -c $cdir/npt1.gro -p topol.top -o prod.tpr -maxwarn 2
                cd $cdir/REST
            done
            mpirun -np 10 gmx mdrun $verbose -deffnm npt1 -multidir {0..9} -replex 800 -plumed plumed.dat -hrex -dlb no
            ;;
        pbc)
            if [ ! -n "$tpr" ]
            then
                echo "Must provide a tpr file with water removed."
                usage
                exit 1
            fi
            if [ ! -f "$tpr" ]
            then
                echo "$tpr does not exist"
                usage
                exit 1
            fi

            for((i=0;i<10;i++))
            do       
                cd $i/
                printf "1\n0\n" | gmx trjconv -s $tpr -f prod.xtc -o whole.xtc -pbc mol -center || ( echo "There is an issue with the provided TPR file, did you strip waters." && usage && exit 1 )
                printf "1\n0\n" | gmx trjconv -s $tpr -f prod.gro -o prod.pdb -pbc mol -center || ( echo "There is an issue with producing prod.pdb." && usage && exit 1 )
                cd $cdir/REST
            done
            echo "Created pbc corrected, whole molecules and centered the protein in the box saved as whole.xtc in each replica directory."
            ;;
        *)
            echo "stage name must match one of: Setup, Min, NVT, NPT, REST"
            usage
            exit
            ;;
    esac

}


# Main

handle_flags "$@"

# Perform the desired actions

if [ "$verbose" = true ]
then
	echo "Verbose mode enable."
    set -x
	set -v
fi

if [ -n "$log" ]
then
	echo "Log file specified: $log"
    exec > >(tee ${log}) 2>&1
fi

simulations
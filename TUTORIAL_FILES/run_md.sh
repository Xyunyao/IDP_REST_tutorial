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
fasta=""
equil_dir="Equilibration"
initstructs="initial_structures"
tutdir=`pwd`


# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
    echo " -g, --generate Generate extended structure with pmx and extract 10 frames from a short vacuum simulation"
    echo " -p, --topol  Topology File name for rest"
    echo " -s, --stage   Select current phase of simulation: Gen, Setup, Min, NVT, NPT, REST, PBC"
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
            -g | --genereate)
                if ! has_argument $@
				then
					echo "Fasta sequence not specified." >&2
					usage
					exit 1
				fi
                fasta=$(extract_argument $@)

                shift
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
        gen)
        python generate_chain.py ${fasta}
        mkdir ${initstructs}
        cd ${initstructs}
        gmx pdb2gmx -f ${tutdir}/extended.pdb -o protein.gro
        gmx editconf -f protein.gro -o vac_box.gro -box 15 15 15
        cp ${tutdir}/mdp_files/VAC.mdp .
        gmx grompp -f VAC.mdp -c vac_box.gro -p topol.top -o vac.tpr
        gmx mdrun -v -deffnm vac
        gmx trjconv -s vac.tpr -f vac.xtc -o starting_structure.pdb -b 20 -dt 20 -sep

        cd ..
        echo "Starting structures prepared under ${tutdir}/${initstructs}"
        ;;
        setup)
            box=6.5
            echo "Solvating each structure with box size ${box} nm."
            mkdir ${equil_dir} && cd ${equil_dir}
            for i in {0..9}
            do
                cp ../${initstructs}/starting_structure${i}.pdb ${i}/
                cp -r ${tutdir}/initial_input_files/a99SBdisp.ff ${i}/
                cp ${tutdir}/mdp_files/{NPT0.mdp,NPT1.mdp,NVT.mdp,minimz.mdp} ${i}/
            done

            for i in {0..9}
            do
                cd ${i}
                printf "1\n1\n" | gmx pdb2gmx -f starting_structures${i}.pdb -o prot.gro 
                gmx editconf -f prot.gro -o box.gro -bt cubic -box ${box} ${box} ${box} -c
                gmx solvate -cp box.gro -cs a99SBdisp.ff/a99SBdisp_water.gro -p topol.top -o solvate.gro
                gmx grompp -f ../minimz.mdp -c solvate.gro -p topol.top -o ions.tpr
                printf "13\n" | gmx genion -s ions.tpr -pname NA -nname CL -neutral -p topol.top -o ions.gro 
                cd ${tutdir}/${equil_dir}
            done
            cd ${tutdir}
            echo "Systems Prepared."
            ;;
        min)
            cd ${tutdir}/${equil_dir}
            for i in {0..9}
            do
                cd ${i}
                if [ ! -f "ions.gro" ]
                then
                    echo "ions.gro missing from replica $i, run the setup stage"
                    exit 1
                fi
                gmx grompp -f mdp_files/minimz.mdp -c system.gro -p topol.top -o min.tpr -maxwarn 2
                cd cd ${tutdir}/${equil_dir}
            done

            echo "Performing System Minimization."
            gmx mdrun ${verbose} -s min.tpr -deffnm min -multidir {0..9}
            echo "Minimized system."
            ;;
        nvt)
            cd ${tutdir}/${equil_dir}
            for i in {0..9}
            do
                cd ${i}
                if [ ! -f "min.gro" ]
                then
                    echo "min.gro missing, run the min stage"
                    exit 1
                fi
                gmx grompp -f mdp_files/NVT.mdp -p topol.top -c min.gro -o nvt.tpr 
                cd ${tutdir}/${equil_dir}
            done
            echo "Thermalizing the system and Equilibrating to 300 K under the NVT Ensemble."
            gmx mdrun ${verbose} -s nvt.tpr -deffnm nvt -mutlidir {0..9}
            echo "System equilibrated under NVT."
            cd ${tutdir}
            ;;
        npt)
            cd ${tutdir}/${equil_dir}
            echo "Pressure and Temperature Equilibration at 1 bar and 300K (NPT)\n"
            echo "stage 1: Short Berendsen pressure relaxation.\n"
            for i in {0..9}
            do
                cd ${i}
                if [ ! -f "nvt.gro" ]
                then
                    echo "nvt.gro missing from replica $i, run the nvt stage"
                    exit 1
                fi
                gmx grompp -f mdp_files/NPT0.mdp -p topol.top -c nvt.gro -o npt0.tpr
                cd ${tutdir}/${equil_dir}
            done
            
            gmx mdrun ${verbose} -s npt0.tpr -deffnm npt0 -multidir {0..9}

            echo "stage 2: Long Parrinello-Rahman pressure Equilibration.\n"
            
            for i in {0..9}
            do
                cd ${i}
                if [ ! -f "npt0.gro" ]
                then
                    echo "nvt.gro missing from replica $i, maybe npt stage 1 failed?"
                    exit 1
                fi
                if [ $i == 0]
                then
                    gmx grompp -f mdp_files/NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr -pp ${tutdir}/processed.top
                else
                    gmx grompp -f mdp_files/NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr
                fi
                cd ${tutdir}/${equil_dir}
            done

            gmx mdrun ${verbose} -s npt1.tpr -deffnm npt1 -multidir {0..9}
            gmx_s select -on -s 0/npt1.tpr -select 'group "Protein" or group "Ion" '
            gmx_s convert-tpr -n index.ndx -o ../non-water.tpr -s 0/npt1.tpr
            echo "System NPT equilibration complete, check volume and pressure convergence."
            echo "Before performing REST2 Simulations, copy the processed.top file in"
            echo "${tutdir}"
            echo "to a file with the name REST.top for editing"
            echo "and add underscores, e.g. HA_, as shown in the README guide." 
            cd ${tutdir}
            ;;
        rest)
            if [ ! -f $topol ]
            then 
                echo "Rest topology does not exist."
                usage
                exit 1
            fi
            cd ${tutdir}/${equil_dir}
            for i in {0..9}
            do
                cd ${i}
                if [ ! -f "npt1.gro" ]
                then
                    echo "npt1.gro missing from replica $i, run the npt stage"
                    exit 1
                fi
                cd ${tutdir}/${equil_dir}
            done
            
            cd ${tutdir}
            echo "Creating rest directory and subdiretories."
            mkdir REST

            cp mdp_files/prod.mdp REST/
            cd REST
            echo "Scaling topologies for each replica and creating run input file prod.tpr"
            for((i=0;i<10;i++))
            do
                mkdir ${i}
                plumed partial_tempering $(lambda_value ${i}) < ${tutdir}/${topol} > ${i}/topol.top
                touch ${i}/plumed.dat
                
                cd $i/
                gmx grompp -f ../prod.mdp -c ${tutdir}/${equil_dir}/${i}/npt1.gro -p topol.top -o prod.tpr -maxwarn 2
                cd ${tutdir}/REST
            done
            mpirun -np 10 gmx mdrun ${verbose} -deffnm prod -multidir {0..9} -replex 800 -plumed plumed.dat -hrex -dlb no
            ;;
        pbc)
            if [ ! -n "${tpr}" ]
            then
                echo "Must provide a tpr file with water removed."
                usage
                exit 1
            fi
            if [ ! -f "${tpr}" ]
            then
                echo "${tpr} does not exist"
                usage
                exit 1
            fi

            for((i=0;i<10;i++))
            do       
                cd ${i}/
                printf "1\n0\n" | gmx trjconv -s ${tpr} -f prod.xtc -o whole.xtc -pbc mol -center || ( echo "There is an issue with the provided TPR file, did you strip waters." && usage && exit 1 )
                printf "1\n0\n" | gmx trjconv -s ${tpr} -f prod.gro -o prod.pdb -pbc mol -center || ( echo "There is an issue with producing prod.pdb." && usage && exit 1 )
                cd ${tutdir}/REST
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

if [ "${verbose}" = true ]
then
	echo "Verbose mode enable."
    set -x
	set -v
fi

if [ -n "${log}" ]
then
	echo "Log file specified: ${log}"
    exec > >(tee ${log}) 2>&1
fi

simulations
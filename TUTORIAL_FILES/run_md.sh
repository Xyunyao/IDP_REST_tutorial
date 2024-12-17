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

verbose_mode=""
output_file=""
ensemble=""


# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
    echo " -p, --topol  Topology File name for rest"
    echo " -e, --ensemble   Select current phase of simulation: Setup, Min, NVT, NPT, REST"
	echo " -l, --log	STDIO log File"
}

has_argument() {
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]] ;
}

get_argument() {
	echo "${2:-${1#*=}}"
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
				verbose_mode=" -v"
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
            -e | --ensemble*)
                if ! has_argument $@
                then
                    echo "When using the ensemble flag you must provide Min, NVT, NPT, or REST"
                    usage
                    exit 1
                fi

                ensemble=$(extract_argument $@)

                shift
                ;;
            -b, --box*)
                if ! has_argument $@
                then
                    echo "When using the box flag you must provide a box length (nm)"
                    usage
                    exit 1
                fi

                box=$(extract_argument $@)

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
    case ${ensemble,,} in
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
            echo "Performing System Minimization."
            gmx grompp -f mdp_files/minimz.mdp -c system.gro -p topol.top -o min.tpr -maxwarn 2
            gmx mdrun $verbose -s min.tpr -deffnm min
            echo "Minimized system."
            ;;
        nvt)
            echo "Thermalizing the system and Equilibrating to 300K under the NVT ensemble."
            gmx grompp -f mdp_files/NVT.mdp -p topol.top -c min.gro -o nvt.tpr 
            gmx mdrun $verbose -s nvt.tpr -deffnm nvt
            echo "System equilibrated under NVT."
            ;;
        npt)
            echo "Pressure and Temperature Equilibration at 1 bar and 300K (NPT)"
            echo "Stage 1: Short Berendsen pressure relaxation."
            gmx grompp -f mdp_files/NPT0.mdp -p topol.top -c nvt.gro -o npt0.tpr 
            gmx mdrun $verbose -s npt0.tpr -deffnm npt0
            echo "Stage 2: Long Parrinello-Rahman pressure Equilibration."
            gmx grompp -f mdp_files/NPT1.mdp -p topol.top -c npt0.gro -o npt1.tpr 
            gmx mdrun $verbose -s npt1.tpr -deffnm npt1
            echo "System NPT equilibration complete, check volume and pressure convergence."
            ;;
        rest)
            echo ""
            ;;
        *)
            echo "Ensemble name must match one of: Setup, Min, NVT, NPT, REST"
            usage
            exit
            ;;
    esac

}


# Main

handle_flags "$@"


# Perform the desired actions

if [ "$verbose_mode" = true ]
then
	echo "Verbose mode enable."
fi

if [ -n "$log" ]
then
	echo "Output file specified: $log"
fi

# Main

handle_flags "$@"


# Perform the desired actions

if [ "$verbose_mode" = true ]
then
	echo "Verbose mode enable."
fi

if [ -n "$output_file" ]
then
	echo "Output file specified: $output_file"
fi



source ${HOME}/.bashrc
module load openmpi/5.0
. .gmx_container.bash

mpi_threads 20
# module list
if [[ ${flags[0]} == '--packmol-populate' ]]
then
	if [[ ! -d "input_gro" ]]
	then
		mkdir input_gro
	fi

	for((i=1;i<=20;i++))
	do
		l=${flags[1]}
		gmx_s editconf -f setup/system_pdbs/system${i}.pdb -o input_gro/system${i}.gro -box $l $l $l -resnr 1
	done
fi

if [[ ${flags[0]} == '-em' ]]
then

for((i=1;i<=20;i++))
do

rm -rf $i
mkdir $i

cp topol_start.top $i/topol.top
cp input_gro/system${i}.gro $i/system.gro

cd $i

echo 1

#
# If packmol was used for input, skip creation of box (assumed already to be in a .gro file 
# with box lengths. Then perform solvation and ion addition (requires ions.mdp)
# 
#
if [[ ! ${flags[1]} == '-packmol' ]]
then
	gmx_s editconf -f system.gro -o system.box.gro -c -box 6.5 6.5 6.5 -bt cubic
	gmx_s solvate -cp system.box.gro -cs ../a99SBdisp.ff/a99SBdisp_water.gro -maxsol 8450 -o system.solv.gro -p topol.top
	gmx_s grompp -f ../mdp/ions.mdp -c system.solv.gro -p topol.top -o system.ions.tpr -maxwarn 2

# echo 13 for apo, and apo 15 for ligand

echo 13 | gmx_s genion -s system.ions.tpr -o system.ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.02

echo -e "1|13\nq\n" | gmx_s make_ndx -f system.ions.gro -o index.ndx

else
	cp system.gro system.ions.gro
fi

gmx_s grompp -f ../mdp/minim.mdp -c system.ions.gro -p topol.top -o system.em.tpr -maxwarn 2

rm \#*
cd ../
done
# echo `pwd`

gmx_gpu_mpi mdrun -v -ntomp 2 -deffnm system.em -multidir {1..20} -dlb no 

fi

if [[ ${flags[0]} == '-eqnvt' ]]
then

for((i=1;i<=20;i++))
do 

   cd $i/
   gmx_s grompp -f ../mdp/nvt.mdp -c system.em.gro -p topol.top -o system.eqnvt.tpr 
   cd ../

done

gmx_gpu_mpi mdrun -v -ntomp 2 -deffnm system.eqnvt -multidir {1..20} -dlb no 

fi

if [[ ${flags[0]} == '-eqnpt' ]]
then

for((i=1;i<=20;i++))
do
	cd $i/
	gmx_s grompp -f ../mdp/npt.mdp -c system.eqnvt.gro -p topol.top -o system.eqnpt.tpr 
	cd ../
done

gmx_gpu_mpi mdrun -v -ntomp 2 -deffnm system.eqnpt -multidir {1..20} -dlb no 

fi

if [[ ${flags[0]} == '-eqnpt2' ]]
then

for((i=1;i<=20;i++))
do
	cd $i/
	gmx_s grompp -f ../mdp/npt2.mdp -c system.eqnpt.gro -p topol.top -o system.eqnpt2.tpr & 
	cd ../
done
wait

gmx_gpu_mpi mdrun -v -deffnm system.eqnpt2 -multidir {1..20} -dlb no 

fi

if [[ ${flags[0]} == '-rest' ]]
then
        if [[ -d structure ]]
        then
            rm -rf structure
        fi
        mkdir structure

for((i=1;i<=20;i++))
do

	cd $i/
	gmx_s editconf -f system.eqnpt2.gro -o ${i}.new.gro -box ${flags[1]} ${flags[1]} ${flags[1]} -bt cubic
	if [[ "$i" == "1" ]]
	then
		gmx_s grompp -f ../mdp/production.mdp -c ${i}.new.gro -p topol.top -pp ../processed.top -o del.tpr 
		rm del.tpr
	fi
        cd ../
	wait
        cp ${i}/${i}.new.gro structure/
done
fi



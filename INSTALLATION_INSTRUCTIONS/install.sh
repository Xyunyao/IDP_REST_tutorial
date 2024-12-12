#!/bin/bash

# Created by Korey Reid 2024/12/01
# This script is written to aid in compiling and installation. 
# This script will provide a local environment, either on a 
# personal computer or an HPC environment. If installing on 
# an HPC environment be carefull to compile this software on
# the compute node you intend to utilize or use flag --slow. 
# If you have a GPU dense node ( e.g. > 4 GPUs) and you want 
# to compile with support for parallel GPU support check 
# the help menu: install.sh <--help or -h>




# Default variables

verbose_mode=false
log_file=""
CMAKE_ARGS="-DGMX_MPI=on -DGMX_BUILD_OWN_FFTW=ON"
INSTALL_ROOT=$HOME/opt
MPICXX=""
make_flag=""
CUDA="-DGMX_USE_CUDA=off"
GMX_version=2022.5
Plumed_version=2.9.1
SIMD="-DGMX_SIMD=AVX2_256"


# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Required Options:"
	echo " -mpi $MPICXX, --mpi=$MPICXX mpi C compiler wrapper full path"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
	echo " -l, --log	Output FILE"
	echo " -j, --ncpu   number of cpu threads to supply make" 
	echo " -simd, --simd architecture type, refer to GMX Manual"
	echo " -c, --cuda	cuda "
    echo "     --installdir (Default) \$HOME/opt"
    echo "     --failsafe   Use SSE4.1"

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
				verbose_mode=true
				;;
			-mpi | --mpi*)
				if ! has_argument $@
				then
					echo "MPICXX variable cannot be empty" >&2
				fi

				MPICXX=$(get_argument $@)

				shift
				;;
			-j | --ncpu*)
				if ! has_argument $@
				then
					echo "-j, --ncpu requires an int value" >&2
				fi

				make_flag="-j $(get_argument $@)"

				shift
				;;
			-simd | --simd*)
				if ! has_argument $@
				then
					echo "-simd, --simd requires string" >&2
				fi

				SIMD=" -DGMX_SIMD=$(get_argument $@)"

				shift
				;;
            -c | --cuda*)
                if ! has_argument $@
                then
                    echo "CUDA_TOOLKIT_DIRECTORY not supplied, e.g. /usr/local/cuda" >&2
                fi

                CMAKE_ARGS+=" -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=$(get_argument $@)"
				CUDA=" -DGMX_USE_CUDA=on"

                shift
                ;;
            -l | --log*)
                if ! has_argument $@
                then
                    echo "Must include file with -l FILE or --log=FILE" >&2
                fi

				log_file="$(get_argument $@)"

                shift
                ;;
            --installdir*)
                if ! has_argument $@
                then
                    echo "When using --installdir you must supply a directory, --installdir=/usr/local/" >&2
                fi

                INSTALL_ROOT="$(get_argument $@)"

                shift
                ;;
            --failsafe)
                SIMD=" -DGMX_SIMD=SSE4.1"
				;;
			*)
				echo "Invalid option:$1" >&2
				usage
				exit 1
				;;
		esac
		shift
	done
	if [ ! -e "$MPICXX" ]
	then
		printf "\n####  Must supply -mpi MPICXX to proceed  ####\n\n" >&2
		usage
		exit 1
	fi

    CMAKE_ARGS+=" -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT}"
	CMAKE_ARGS+=" ${CUDA}"
	CMAKE_ARGS+=" ${SIMD}"
}


# Main

handle_flags "$@"


if [ -n "$log_file" ]
then
	echo "Log file specified: $log_file"
	exec > >(tee ${log_file}) 2>&1
fi

if [ -n "$verbose_mode" ]
then
	set -x
	set -v
fi

PS4='${LINENO}: '

eval "$(conda shell.bash hook)" 
conda activate REST_tutorial

mkdir src
cd src

# Clone both gromacs and plumed2
git clone https://github.com/plumed/plumed2.git
git clone https://github.com/gromacs/gromacs.git

# Checkout gromacs 2022.5
cd gromacs

git checkout -b v${GMX_version}

# Configure and install plumed2
# Plumed2 and gromacs require MPI and the MPICXX environment variable set
cd ../plumed2
git checkout -b v${Plumed_version}
./configure --prefix=${INSTALL_ROOT} --enable-modules=all --enable-mpi CXX="$MPICXX" 

# Compile Plumed2, this is a dirty install, normally following compiling 'make check' would follow 
make ${make_flag} && make install && echo "Plumed2 compiled successfully" || (echo "Plumed2 compilation failed" && exit 1)

# Add export commands to bashrc file and source bashrc file
echo 'export PLUMED_KERNEL=$HOME/opt/lib/libplumedKernel.so' > ../../.bashrc
echo 'export PATH=$HOME/opt/bin:$HOME/opt/lib:$PATH' >> ../../.bashrc
echo 'export PLUMED_ROOT=$HOME/opt/lib/plumed' >> ../../.bashrc
echo 'export LD_LIBRARY_PATH="$HOME/opt/lib/:$LD_LIBRARY_PATH"' >> ../../.bashrc
echo 'export LIBRARY_PATH="$HOME/opt/lib/:$LIBRARY_PATH"' >> ../../.bashrc
cd ..
source ../.bashrc

# Patch gromacs with plumed patch command
cd gromacs
${INSTALL_ROOT}/bin/plumed patch -p -e gromacs-${GMX_version} || (echo "patching failed, check log" && exit 1)

# Run cmake command and compile patched gromacs
mkdir build
cd build
cmake .. ${CMAKE_ARGS}
make ${make_flag} && make install && echo "gromacs successfully compiled with plumed2" || (echo "compilation of gromacs with plumed failed" >&2)

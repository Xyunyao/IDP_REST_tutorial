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
output_file=""
CMAKE_ARGS="-DGMX_MPI -DGMX_BUILD_OWN_FFTW=ON"
INSTALL_ROOT=$HOME/opt



# Function to display script usage

usage() {
	echo "Usage: $0 [Options]"
	echo "Options:"
	echo " -h, --help	Display this help message"
	echo " -v, --verbose	Enable verbosity"
	echo " -o, --ofile	Output FILE"
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
			-f | --file*)
				if ! has_argument $@
				then
					echo "File not specified." >&2
					usage
					exit 1
				fi

				output_file=$(extract_argument $@)

				shift
				;;
            -c | --cuda*)
                if ! has_argument $@
                then
                    echo "CUDA_TOOLKIT_DIRECTORY not supplied, e.g. /usr/local/cuda" >&2
                fi

                CMAKE_ARGS+=" -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=$(extract_argument $@)"

                shift
                ;;
            -l | --log*)
                if ! has_argument $@
                then
                    echo "Must include file with -l FILE or --log=FILE" >&2
                fi



                shift
                ;;
            --installdir*)
                if ! has_argument $@
                then
                    echo "When using --installdir you must supply a directory, --installdir=/usr/local/" >&2
                fi

                INSTALL_ROOT="$(extract_argument $@)"

                shift
                ;;
            --failsafe)
                CMAKE_ARGS+=" -DGMX_SIMD=SSE4.1"
			*)
				echo "Invalid option:$1" >&2
				usage
				exit 1
				;;
		esac
		shift
	done
    CMAKE_ARGS+=" -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT}"
}


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


mkdir src
cd src

# Clone both gromacs and plumed2
git clone https://github.com/plumed/plumed2.git
git clone https://github.com/gromacs/gromacs.git

# Checkout gromacs 2024.3
cd gromacs
GMX_version=2024.3
git checkout -b v${GMX_version}

# Configure and install plumed2
# Plumed2 and gromacs require MPI and the MPICXX environment variable set
cd ../plumed2
./configure --prefix=${INSTALL_ROOT}/opt --enable-modules=all CXX="$MPICXX" CXXFLAGS="-O3 -axSSE2,AVX" 

# Compile Plumed2
make && make check && make install && echo "Plumed2 compiled successfully" || echo "Plumed2 compilation failed" && exit 1

# Add export commands to bashrc file and source bashrc file
echo 'export PLUMED_KERNEL="$HOME/opt/lib/libplumedKernel.so"' >> ~/.bashrc
echo 'export PATH=$HOME/opt/bin:$PATH' >> ~/.bashrc
echo 'export PLUMED_ROOT=$HOME/opt'
source ~/.bashrc

# Patch gromacs with plumed patch command
cd ../gromacs
plumed patch -e gromacs-${GMX_version}

# Run cmake command and compile patched gromacs
mkdir build
cd build
cmake .. ${CMAKE_ARGS}
make && make check && make install && echo "gromacs successfully compiled with plumed2" || (echo "compilation of gromacs with plumed failed" >&2)

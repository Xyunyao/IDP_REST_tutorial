#!/bin/bash

# ===============================
# Cluster-ready GROMACS + Plumed2 installer using Conda GCC + MPI
# ===============================

# Default variables
INSTALL_ROOT=$HOME/opt
GMX_VERSION=2022.5
PLUMED_VERSION=2.9.1
CMAKE_ARGS="-DGMX_MPI=ON -DGMX_BUILD_OWN_FFTW=ON"
CUDA="-DGMX_USE_CUDA=OFF"
SIMD="-DGMX_SIMD=AVX2_256"
MAKE_JOBS=4
LOG_FILE=""
VERBOSE=false

# ===============================
# Functions
# ===============================
usage() {
    echo "Usage: $0 [Options]"
    echo "Options:"
    echo " -j <N>           Number of CPU threads for make"
    echo " -simd <ARCH>     SIMD type, e.g. AVX2_256"
    echo " -c <CUDA_PATH>   Enable CUDA with toolkit path"
    echo " -l <LOG_FILE>    Log output"
    echo " -h               Show help"
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -j)
            MAKE_JOBS=$2
            shift 2
            ;;
        -simd)
            SIMD="-DGMX_SIMD=$2"
            shift 2
            ;;
        -c)
            CUDA="-DGMX_USE_CUDA=ON -DCUDA_TOOLKIT_ROOT_DIR=$2"
            shift 2
            ;;
        -l)
            LOG_FILE=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# ===============================
# Logging
# ===============================
if [ -n "$LOG_FILE" ]; then
    exec > >(tee "$LOG_FILE") 2>&1
fi

[ "$VERBOSE" = true ] && set -x

# ===============================
# Activate Conda environment
# ===============================
eval "$(conda shell.bash hook)"
conda activate REST_tutorial

# ===============================
# Ensure Conda GCC + MPI + FFTW + pkg-config
# ===============================
conda install -y -c conda-forge gcc_linux-64 gxx_linux-64 mpich fftw pkg-config libgcc-ng libgfortran

export CC=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
export MPICXX=$CONDA_PREFIX/bin/mpicxx
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export PATH=$CONDA_PREFIX/bin:$PATH

# ===============================
# Prepare source directories
# ===============================
mkdir -p src
cd src

# Clone repositories if not already present
[ ! -d plumed2 ] && git clone https://github.com/plumed/plumed2.git
[ ! -d gromacs ] && git clone https://github.com/gromacs/gromacs.git

# ===============================
# Compile Plumed first
# ===============================
cd plumed2
git checkout -B v${PLUMED_VERSION}

# Clean previous build
make clean
git clean -xfd

# Configure Plumed
./configure --prefix=${INSTALL_ROOT} --enable-modules=all --enable-mpi CXX=$MPICXX

# Compile and install
make -j${MAKE_JOBS} && make install || { echo "Plumed2 compilation failed"; exit 1; }

# ===============================
# Set environment for Plumed
# ===============================
export PLUMED_ROOT=${INSTALL_ROOT}
export PLUMED_KERNEL=${INSTALL_ROOT}/lib/libplumedKernel.so


if [ ! -f "$PLUMED_ROOT/bin/plumed" ] || [ ! -f "$PLUMED_KERNEL" ]; then
    echo "ERROR: Plumed was not installed correctly. Check build logs."
    exit 1
fi

export PATH=${PLUMED_ROOT}/bin:$PATH
export LD_LIBRARY_PATH=${PLUMED_ROOT}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${PLUMED_ROOT}/lib:$LIBRARY_PATH

# ===============================
# Patch GROMACS with Plumed
# ===============================
cd ../gromacs
git checkout -B v${GMX_VERSION}


# ===============================
# Build GROMACS
# ===============================
mkdir -p build
cd build

CMAKE_ARGS+=" ${CUDA} ${SIMD} -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} \
             -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DPLUMED_ROOT=${PLUMED_ROOT}"

cmake .. ${CMAKE_ARGS} || { echo "CMake configuration failed"; exit 1; }

make -j${MAKE_JOBS} && make install || { echo "GROMACS compilation failed"; exit 1; }

echo "GROMACS + Plumed installation complete!"

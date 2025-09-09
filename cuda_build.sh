#!/bin/bash
# this is the script that successfully installed cuda_mpi_plumed patched gromacs. This script assume plumed has been installed.
# ===============================
# Cluster-ready GROMACS installer with CUDA + MPI (using existing Plumed)
# ===============================

# ===============================
# User settings
# ===============================
INSTALL_ROOT=$HOME/opt/gromacs_cuda   # new folder for CUDA build
GMX_VERSION=2022.5                     # or 2026 dev
MAKE_JOBS=8
#SIMD="-DGMX_SIMD=AVX2_256"
CUDA_PATH=/cm/shared/apps/cuda12.2/toolkit/12.2.2   # adjust to your CUDA toolkit

export PATH=${CUDA_PATH}/bin:$PATH
export LD_LIBRARY_PATH=${CUDA_PATH}/lib64:$LD_LIBRARY_PATH

LOG_FILE=""                            # optional log file
VERBOSE=false

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
# Ensure Conda GCC + MPI + FFTW
# ===============================
conda install -y -c conda-forge  mpich fftw pkg-config libgcc-ng libgfortran gcc_linux-64=11 gxx_linux-64=11


export CC=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
export CFLAGS="-mavx2"
export CXXFLAGS="-mavx2"
export MPICXX=$CONDA_PREFIX/bin/mpicxx
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export PATH=$CONDA_PREFIX/bin:$PATH

# ===============================
# Use existing Plumed installation
# ===============================
export PLUMED_ROOT=$HOME/opt
export PLUMED_KERNEL=$PLUMED_ROOT/lib/libplumedKernel.so
export PATH=$PLUMED_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$PLUMED_ROOT/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$PLUMED_ROOT/lib:$LIBRARY_PATH

if [ ! -f "$PLUMED_ROOT/bin/plumed" ] || [ ! -f "$PLUMED_KERNEL" ]; then
    echo "ERROR: Plumed not found in $PLUMED_ROOT. Exiting."
    exit 1
fi

# ===============================
# Prepare GROMACS source
# ===============================
mkdir -p src
cd src

#[ ! -d gromacs ] && git clone https://github.com/gromacs/gromacs.git
# ================================
# Download GROMACS source
# ================================
if [ ! -d gromacs-${GMX_VERSION} ]; then
    wget https://ftp.gromacs.org/gromacs/gromacs-${GMX_VERSION}.tar.gz
    tar -xvzf gromacs-${GMX_VERSION}.tar.gz
fi
cd gromacs-${GMX_VERSION}

${PLUMED_ROOT}/bin/plumed patch -p -e gromacs-${GMX_VERSION} || (echo "patching failed, check log" && exit 1)
#git checkout -B v${GMX_VERSION}

# ===============================
# Detect GPU architecture automatically
# ===============================
if command -v nvidia-smi >/dev/null 2>&1; then
    SM_MAJOR=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -n1 | cut -d'.' -f1)
    SM_MINOR=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader,nounits | head -n1 | cut -d'.' -f2)
    SM_VERSION="${SM_MAJOR}${SM_MINOR}"
    echo "Detected GPU architecture: SM $SM_VERSION"
else
    echo "Warning: nvidia-smi not found. Using default CUDA SM."
    SM_VERSION="Auto"
fi

# ===============================
# Build GROMACS
# ===============================
mkdir -p build_cuda
cd build_cuda


CMAKE_ARGS="-DGMX_MPI=ON \
            -DGMX_BUILD_OWN_FFTW=ON \
            -DGMX_GPU=CUDA \
            -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_PATH} \
            -DGMX_CUDA_TARGET_SM= "86" \
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
            -DPLUMED_ROOT=${PLUMED_ROOT} \
            -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} \
            -DCMAKE_C_COMPILER=$CC \
            -DGMX_SIMD=AUTO \
            -DCMAKE_CXX_COMPILER=$CXX"
            



cmake .. ${CMAKE_ARGS} || { echo "CMake configuration failed"; exit 1; }

make -j${MAKE_JOBS} && make install || { echo "GROMACS compilation failed"; exit 1; }

# ===============================
# Finished
# ===============================
echo "GROMACS CUDA + MPI installation complete!"
echo "Use: source $INSTALL_ROOT/bin/GMXRC to set environment variables"


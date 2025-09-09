#!/bin/bash

# ===============================
# Cluster-ready GROMACS installer with CUDA + MPI (using existing Plumed)
# ===============================

# ===============================
# User settings
# ===============================
INSTALL_ROOT=$HOME/opt/gromacs_cuda   # new folder for CUDA build
GMX_VERSION=2022.5                     # or 2026 dev
MAKE_JOBS=8
SIMD="-DGMX_SIMD=AVX2_256"
CUDA_PATH=/cm/shared/apps/cuda12.2/toolkit/12.2.2/bin/nvcc    # adjust to your CUDA toolkit
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
conda install -y -c conda-forge gcc_linux-64 gxx_linux-64 mpich fftw pkg-config libgcc-ng libgfortran

export CC=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++
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
[ ! -d gromacs ] && git clone https://github.com/gromacs/gromacs.git
cd gromacs
git checkout -B v${GMX_VERSION}

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
    SM_VERSION=""
fi

# ===============================
# Build GROMACS
# ===============================
mkdir -p build_cuda
cd build_cuda

CMAKE_ARGS="-DGMX_MPI=ON \
            -DGMX_BUILD_OWN_FFTW=ON \
            -DGMX_USE_CUDA=ON \
            -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_PATH} \
            -DPLUMED_ROOT=${PLUMED_ROOT} \
            ${SIMD} \
            -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT} \
            -DCMAKE_C_COMPILER=$CC \
            -DCMAKE_CXX_COMPILER=$CXX"

# Add SM version if detected
if [ -n "$SM_VERSION" ]; then
    CMAKE_ARGS+=" -DGMX_CUDA_TARGET_SM=${SM_VERSION}"
fi

cmake .. ${CMAKE_ARGS} || { echo "CMake configuration failed"; exit 1; }

make -j${MAKE_JOBS} && make install || { echo "GROMACS compilation failed"; exit 1; }

# ===============================
# Finished
# ===============================
echo "GROMACS CUDA + MPI installation complete!"
echo "Use: source $INSTALL_ROOT/bin/GMXRC to set environment variables"

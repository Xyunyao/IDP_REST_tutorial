#!/bin/bash

INSTALL_ROOT=$HOME/opt

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
echo 'export PLUMED_KERNEL="/usr/share/lib/libplumedKernel.so"' >> ~/.bashrc
echo 'export PATH=$HOME/opt/bin:$PATH' >> ~/.bashrc
echo 'export PLUMED_ROOT=$HOME/opt'
source ~/.bashrc

# Patch gromacs with plumed patch command
cd ../gromacs
plumed patch -e gromacs-${GMX_version}

# Run cmake command and compile patched gromacs
mkdir build
cd build
cmake .. -DGMX_THREAD_MPI=OFF -DGMX_MPI=ON -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=${INSTALL_ROOT}
make && make check && make install && echo "gromacs successfully compiled with plumed2" || echo "compilation of gromacs with plumed failed"

#!/usr/bin/env bash

module load local-fluidity-gcc


echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
echo ""
echo " BUILDING EVERYTHING ..."
echo ""
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="

if [ "$CONTRIB_DIR" != "" ]
then
    echo "------------------------------------"
    echo "Cleaning out old Zoltan and PETSc..."

    if [ "$PETSC_DIR" != "" ]
    then
	rm -rf "$PETSC_DIR"
    fi
    
    if [ "$ZOLTAN_DIR" != "" ]
    then
	rm -rf "$ZOLTAN_DIR"
    fi
fi

cd ~/src

# TCMALLOC
echo "------------------------------------"
echo " Building tcmalloc libraries"
echo "------------------------------------"

cd ~/src/gperftools-2.9.1
make clean
./configure --with-tcmalloc-alignment=16 --with-tcmalloc-pagesize=32 --with-pic --enable-minimal --prefix=$HOME/contrib/gcc73/base/1.0 
make -j 8
make install
cd ..

# HDF5 and NetCDF
echo "------------------------------------"
echo " Building HDF5 and NetCDF libraries"
echo "------------------------------------"

bash ~/store/build-hdf-netcdf.sh


# OpenBLAS
echo "------------------------------------"
echo " Building OpenBLAS"
echo "------------------------------------"

cd ~/src
rm -rf OpenBLAS
tar xzf ~/store/OpenBLAS-0.3.15.tar.gz
mv OpenBLAS-0.3.15 OpenBLAS
cd OpenBLAS
make TARGET=HASWELL CC=$CC FC=$FC BINARY=64 USE_THREAD=0 USE_OPENMP=0 USE_LOCKING=1 BUFFERSIZE=25
make PREFIX=$CONTRIB_DIR/openblas install
cd ..


# PETSc
echo "------------------------------------"
echo " Building PETSc libraries"
echo "------------------------------------"

cd petsc-3.6.4
bash ~/store/petsc-build-gcc.sh
cd ..



# Zoltan
echo "------------------------------------"
echo " Building Zoltan libraries"
echo "------------------------------------"

mkdir -p zolbuild
rm -rf zolbuild/*

cd zolbuild
bash ~/store/zoltan-build-gcc.sh ../Zoltan_v3.8
cd ..


# Fluidity
echo "------------------------------------"
echo " Building Fluidity"
echo "------------------------------------"

cd ~/fluidity
bash ~/store/fluidity-build-gcc.sh cto=/path/to/example-sim.flml
cd ..

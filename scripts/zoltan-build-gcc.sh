#!/bin/bash --login

if [ "$ZOLTAN_DIR" == "" ]
then
    echo "**** error: please define ZOLTAN_DIR environmental variable"
    exit 1
fi

# xflags="-I/usr/include/scotch"
xflags=""

# --prefix=$HOME/contrib/gcc49/zoltan/3.8 \

../Zoltan_v3.8/configure x86_64-linux-gnu \
--prefix=$ZOLTAN_DIR \
--enable-f90interface=yes \
--enable-zoltan-cppdriver \
--enable-mpi --with-mpi-compilers \
--disable-examples \
--with-scotch --with-scotch-incdir=$PETSC_DIR/include \
CC=mpicc CXX=mpicxx FC=mpif90 \
CFLAGS="$CFLAGS $xflags" \
CXXFLAGS="$CXXFLAGS $xflags" \
FCFLAGS="$FCFLAGS $xflags" \
LIBS="-lptscotch -lscotch -lscotcherr -lscotcherrexit"

# --with-parmetis=yes \
# --with-parmetis-incdir=$PETSC_DIR/include \
# --with-parmetis-libdir=$PETSC_DIR/lib \


make clean
make -j 4
make install

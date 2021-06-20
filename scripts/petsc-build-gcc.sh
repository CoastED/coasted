PETSC_INSTALL_DIR="$PETSC_DIR"

petsc_prec_opts=""

if [ "$1" == "single" ]
then
    echo "Building single precision."
    petsc_prec_opts="--with-precision=single"
fi

export PETSC_DIR=`pwd`

EXTRA_LIBS="" # -lparmetis -lmetis" 

xflags="-O2 -fpic $CFLAGS"
# xflags="-O2 -ffast-math -march=native -fpic"
# xflags="-O2 -ffast-math -march=native -funroll-loops -fpic"
# -falign-labels=32 -falign-functions=32 \
# -mpreferred-stack-boundary=8"

export LIBRARY_PATH="$LIBRARY_PATH"
export LDFLAGS="$LDFLAGS"

export FC=mpif90

rm -rf arch-linux2-c-opt


./configure --prefix="$PETSC_INSTALL_DIR" "$petsc_prec_opts" \
--with-fortran-interfaces=1 \
--with-fortran-kernels=1 \
--with-shared-libraries=1 \
--with-ssl=0 \
--with-pthread=0 \
--with-openmp=0 \
--with-batch=0 \
--with-single-library=0 \
--with-mpi-compilers=1 \
--with-pic=1 --with-debugging=0 \
--with-blas-lapack-lib=$BLAS_LIB_DIR/libopenblas.so \
--with-metis --download-metis=1 \
--with-parmetis --download-parmetis=1 \
--with-hypre=1 --download-hypre=1 \
--with-scotch=1 --download-scotch=1 \
--with-ptscotch=1 --download-ptscotch=1 \
--CC=mpicc --CXX=mpicxx --FC=mpif90 \
--CFLAGS="$CFLAGS $xflags" \
--CXXFLAGS="$CXXFLAGS $xflags" \
--COPTFLAGS="$CFLAGS $xflags" \
--CXXOPTFLAGS="$CXXFLAGS $xflags" \
--FFLAGS="$FCFLAGS $xflags" \
--FOPTFLAGS="$FCFLAGS $xflags"
--LIBS="$LIBS $EXTRA_LIBS -lopenblas" 


# --with-hypre-dir=$CONTRIB_DIR/hypre \
# --with-blas-lapack-dir=$CONTRIB_DIR/base/default/lib \
# --with-single-library=1 
# --with-parmetis=1 \
# --with-metis=1 \ 


# --with-blas-lapack=$CONTRIB_DIR/base/default/lib \
# --with-blas-lapack-lib=$CONTRIB_DIR/base/default/lib/libopenblas.a \	    
# --with-fortran-kernels=0 \

make MAKE_NP=16 PETSC_DIR=$PETSC_DIR all

make install

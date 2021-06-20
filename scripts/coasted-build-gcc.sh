#!/bin/bash --login

module load local-fluidity-gcc

cd ~/fluidity

debuggingConfig=""
precisionFlag=""
debugCompile=""
fcheckCompile=""
# petscDoubleOnlyLibs="-lHYPRE -lptscotch -lscotch -lscotcherr -lscotcherrexit"
petscDoubleOnlyLibs="-lHYPRE -lpetscsnes -lptscotch -lscotcherrexit -lmetis -lpetscsys -lptscotcherr -lparmetis -lpetsctao -lptscotcherrexit -lpetscdm -lpetscts -lptscotchparmetis -lpetscksp -lpetscvec -lscotch -lpetscmat -lptesmumps -lscotcherr"

for arg in $@
do
    case $arg in
	"debug" )
	    echo "CONFIGURING FOR DEBUG"
    	    # --enable-debugging doesn't work with extruded meshes. Sheesh.
	    debuggingConfig="" 
	    debugCompile="-ggdb -fbacktrace -fstack-usage -fstack-check"
	    ;;
	"check" )
	    echo "Adding memory checks"
	    fcheckCompile="-fcheck=bounds -fcheck=all -fstack-usage -fstack-check"
            ;;
        cto=* )
            echo "Will try to use compile-time optimisations"
            flml_file="`echo $arg | sed 's/^cto=//g'`"
            if [ "$flml_file" != "" ] 
            then
                ctoConfig="--enable-cto=$flml_file"
            fi
            ;;

	"single" )
	    echo "Compiling in single-precision (32-bit) mode"
	    precisionFlag="--enable-dp=no"
	    petscDoubleOnlyLibs=""
	    ;;

    esac
done

echo "---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---"

# make clean
make distclean

mv -f Makefile Makefile.bak

cpu_type=native

xflags="$precisionFlag"
# xflags="-mpreferred-stack-boundary=4 $precisionFlag"
# xflags="$debugCompile $fcheckCompile \
# -fopenmp -march=$cpu_type -mtune=$cpu_type \
# -DHAVE_LIBNUMA -Wunused-result"

cp -f ~/libturb3/interfaceturb.o main/


./configure --enable-openmp \
--with-debugging=0 \
$ctoConfig \
--enable-shared \
--with-x=no \
--with-adjoint=no \
--enable-petsc-fortran-modules \
PYTHON_VERSION="2.7" \
CC=mpicc CXX=mpicxx FC=mpif90 \
MPICC=mpicc MPICXX=mpicxx MPIF77=mpif90 MPIF90=mpif90 \
CFLAGS="$xflags $fcheckCompile $CFLAGS" \
CXXFLAGS="$xflags $fcheckCompile $CXXFLAGS" \
FCFLAGS="$xflags $fcheckCompile $FCFLAGS" \
LIBS="$LIBS $petscDoubleOnlyLibs -lexpat -lnetcdf -lzoltan \
$PETSC_LIBS -lopenblas -lm" \
LDFLAGS="$LIBS"
#
#-lptscotch -lscotch -lscotcherr -lscotcherrexit -lzoltan" \

# --enable-2d-adaptivity \


if test -f Makefile
then
    cd libspud
    ./configure --prefix=`pwd`/.. --enable-shared=no --with-pic
    cd ..

    make -j 16
    make -j 16 fltools
fi


if [ "$1" == "debug" ] && [ -f $HOME/fluidity/bin/fluidity ]
then
    echo "Generating structure file with hpcstruct"
    hpcstruct -I$HOME/fluidity/+ -o $HOME/fluidity/bin/fluidity.hpcstruct \
	$HOME/fluidity/bin/fluidity
fi

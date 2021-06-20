bdirs="hdf5-1.8.15 netcdf-4.3.2 netcdf-c-4.7.3 netcdf-cxx-4.2 \
netcdf-fortran-4.2"

for bdir in $bdirs
do
    cd $bdir
    make distclean
    make clean
    ./configure --prefix=$CONTRIB_DIR/base/1.0
    make -j 8
    make install
    cd ..
done

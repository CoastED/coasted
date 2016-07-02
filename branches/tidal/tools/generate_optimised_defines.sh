#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Batch script to scan an FLML file and generate a file with #defines
# for compile-time optimisation
#
# A. Creech <a.creech@ed.ac.uk>, 7th Aug 2016
# ----------------------------------------------------------------------------

opt_flml_file="$1"

if ! test -f "$opt_flml_file"
then
	>&2 echo "error: cannot find FLML file '$opt_flml_file'"
	exit 1
fi

opt_exec_dir=`dirname "$opt_flml_file"`
opt_compile="yes"

# Pull out some parameters from the FLML file.

opt_dimension=`grep -A 2 dimension\> "$opt_flml_file" | grep integer_value | sed 's/.*<integer_value.*0\">\(.*\)<\/integer_value>/\1/g'`
opt_quad_degree=`grep -A 3 quadrature\> "$opt_flml_file" | grep integer_value | sed 's/.*<integer_value.*0\">\(.*\)<\/integer_value>/\1/g'`
opt_surface_degree=`grep -A 2 surface_degree\> "$opt_flml_file" | grep integer_value | sed 's/.*<integer_value.*0\">\(.*\)<\/integer_value>/\1/g'`

# This isn't always needed in a simulation, but needs to be set.

if [ "$opt_surface_degree" == "" ] 
then
	opt_surface_degree=$opt_quad_degree
fi

# Disable 1D and 2D for now.

if [ "$opt_dimension" != "3" ]
then
	>&2 echo "Error. Currently, only optimisation for 2D is supported."
	exit 1
fi

# We need to have Velocity and Pressure meshes.

if [ "`grep -A 40 geometry $opt_flml_file  | grep -A 20 VelocityMesh`" == "" ]
then
	>&2 echo "Error. VelocityMesh does not exist in $opt_flml_file."
	exit 1
fi		 
if [ "`grep -A 40 geometry $opt_flml_file  | grep -A 20 PressureMesh`" == "" ]
then
	>&2 echo "Error. PressureMesh does not exist in $opt_flml_file."
	exit 1
fi		 

# Pull out element degrees for velocity and pressure

vel_ele_degree=`grep -A 40 geometry $opt_flml_file  | grep -A 20 VelocityMesh | grep -A 10 polynomial_degree | grep integer_value |  sed 's/.*<integer_value.*0\">\(.*\)<\/integer_value>/\1/g'`
pres_ele_degree=`grep -A 40 geometry $opt_flml_file  | grep -A 20 PressureMesh | grep -A 10 polynomial_degree | grep integer_value |  sed 's/.*<integer_value.*0\">\(.*\)<\/integer_value>/\1/g'`




# We're assuming the velocity and pressure elements are geometrically
# based upon tetrahedra. Currently only first and second order elements
# make any sense. This really needs to be revisted, as the numbers are
# wholly dependent on element type and representation.

case "$opt_dimension" in
	1) opt_nfaces=1
		((opt_nloc=2+vel_ele_degree-1))
		opt_floc=1
		opt_p_floc=2
		;;
    2) opt_nfaces=3
    	((opt_nloc=3*vel_ele_degree))
    	((opt_floc=vel_ele_degree+1))
    	((opt_p_floc=pres_ele_degree+1))
    	;;
     3) opt_nfaces=4
    	((opt_nloc=4+(vel_ele_degree-1)*6))
    	((opt_floc=3*vel_ele_degree))
    	((opt_p_floc=3*pres_ele_degree))
		;;
     *) >&2 echo "Error. Only 1-3D supported"
esac

# Calculating NGI for tets

case "$opt_quad_degree" in
     1) opt_ele_ngi=1;;
     2) opt_ele_ngi=4;;
     3) opt_ele_ngi=5;;
     4) opt_ele_ngi=11;;
     5) opt_ele_ngi=14;;
     6) opt_ele_ngi=24;;
     7) opt_ele_ngi=31;;
     8) opt_ele_ngi=43;;
     *) >&2 echo "Error. quad_degree '$opt_quad_degree' unsupported."
     	exit 1;;
esac

# Now calculate the surface NGI (triangles).

case "$opt_surface_degree" in
     1) opt_fngi=1;;
     2) opt_fngi=3;;
     3) opt_fngi=4;;
     4) opt_fngi=6;;
     5) opt_fngi=7;;
     6) opt_fngi=12;;
     7) opt_fngi=12;;
     8) opt_fngi=16;;
     *)	>&2 echo "Error. quad_degree '$opt_surface_degree' unsupported."
     	exit 1;;
esac

# Now write defines to include file compile_opt_defs.h
# This will be used in Fortran only, so '!' comments allowed.

cat > include/compile_opt_defs.h << EOF
! ----------------------------------------------------------------------------
! Compile-time optimisation definitions file. Based upon FLML file
! $opt_flml_file
!
! Generated at `date`
! ----------------------------------------------------------------------------

#define OPT_NDIM $opt_dimension
#define OPT_NFACES $opt_nfaces
#define OPT_NLOC $opt_nloc
#define OPT_NGI $opt_ele_ngi
#define OPT_FLOC $opt_floc
#define OPT_FNGI $opt_fngi
#define OPT_P_FLOC $opt_p_floc
#define OPT_EFLOC ( OPT_NLOC + OPT_NFACES * OPT_FLOC )
EOF

exit 0


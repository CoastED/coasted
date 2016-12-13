!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

module calculate_mesh_metrics
    use elements
    use sparse_tools
    use fetools
    use fefields
    use fields
    use sparse_matrices_fields
    use state_module
    use shape_functions
    use transform_elements
    use vector_tools
    use fldebug
    use vtk_interfaces
    use Coordinates
    use spud
    use boundary_conditions
    use boundary_conditions_from_options
    use solvers
    use dgtools
    use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN, COLOURING_DG2, &
        COLOURING_DG0, new_mesh_geometry, new_mesh_connectivity
    use coriolis_module
    use halos
    use petsc_tools
    use diagnostic_fields
    use smoothing_module
    use fields_manipulation
    use field_derivatives
    use field_options

    implicit none

contains
    subroutine xxxx()
    end subroutine xxxx

end module calculate_mesh_metrics

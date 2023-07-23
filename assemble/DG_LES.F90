!    Copyright (C) 2008 Imperial College London and others.
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
#include "fdebug.h"
#include "compile_opt_defs.h"

! Filter types:
! 1 = original LES filter (element averages and local peaks (quite diffuse)
! 2 = box filter
! 3 = hat type filter

! Currently only implemented for anisotropic LES
! Box filter not fully implemented




module dg_les
  !!< This module contains several subroutines and functions used to implement LES models
  use state_module
  use fetools
  use fields
  use field_options
  use field_derivatives
  use smoothing_module
  use fefields
  use vector_tools
  use state_fields_module
  use solvers
  use global_parameters
  use spud
  use halos

  implicit none

  private

  public :: calc_dg_sgs_scalar_viscosity, calc_dg_sgs_amd_viscosity
  public :: calc_dg_sgs_vreman_viscosity, calc_dg_sgs_roman_viscosity

  ! This scales the Van Driest effect
  real, parameter :: van_scale=1.0
  integer :: smooth_called_count=0

#include "mpif.h"


contains


    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    subroutine calc_dg_sgs_amd_viscosity(state, x, u)
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        character(len=OPTION_PATH_LEN) :: scalar_eddy_visc_path
        character(len=256) :: mesh_name

        ! Currently goes straight to nodal AMD LES

        FLAbort("DG AMD LES is no longer supported.")

    end subroutine calc_dg_sgs_amd_viscosity



    subroutine calc_dg_sgs_scalar_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        type(scalar_field), pointer :: dist_to_wall
        type(scalar_field), pointer :: sgs_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(scalar_field), pointer :: artificial_visc

        type(vector_field), pointer :: elelen, nodelen, smoothlen
        type(tensor_field) :: u_grad
        real, dimension(u%dim, u%dim):: rate_of_strain, u_grad_node
        
        integer :: i, j, k, m, e, num_elements, n, num_nodes
        integer :: u_cg_ele(ele_loc(u,1))

        real :: Cs
        logical :: have_van_driest


        ! Standard LES vars
        real :: visc_turb, tmp_visc, Cs_length_sq
        real :: mu, rho
        integer :: state_flag

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_reference_density, have_lengths_field

        real :: sgs_visc_val
        logical :: have_artificial_visc

        ! Constants / variablesfor Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0
        real :: vd_damping, y_plus

        print*, "In calc_dg_sgs_scalar_viscosity()"


        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")


        if(have_van_driest) then
           dist_to_wall=>extract_scalar_field(state, "DistanceToWall", stat=state_flag)
           if (state_flag/=0) then
              FLAbort("DG_LES: Van Driest damping requested, but no DistanceToWall scalar field exists")
           end if
           
        end if


        ! Velocity projected to continuous Galerkin
        ! u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)
        ! Nonlinear velocity is better...?
        u_cg=>extract_vector_field(state, "ProjectedNonlinearVelocity", stat=state_flag)


        ! Allocate gradient field and calculate gradient
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")
        call grad(u_cg, x, u_grad)
        call halo_update(u_grad)

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for Smagorinsky DG LES. (This should not happen)")
        end if


        ! Length scales for filter
        elelen => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodelen => extract_vector_field(state, "NodeLengthScales", stat=state_flag)
        smoothlen => extract_vector_field(state, "SmoothedLengthScales", stat=state_flag)

        if(state_flag == 0) then
           print*, "**** Have ElementLengthScales and NodeLengthScales fields"
           have_lengths_field = .true.
        else
           FLAbort("Error: must have ElementLengthsScales and NodeLengthScales fields for Smagorinsky DG LES. This should have been automatically created.")
        end if


        ! We can use this in areas of insufficient resolution
        have_artificial_visc = .false.
        artificial_visc => extract_scalar_field(state, "ArtificialViscosity", &
             stat=state_flag)
        if(state_flag == 0) then
           print*, "ArtificialViscosity field detected."
           have_artificial_visc = .true.
        end if


        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

        if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
            .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
            mu=mviscosity%val(1,1,1)
        else
            FLAbort("DG_LES: must have constant or normal viscosity field")
        end if


        ! We only use the reference density. This assumes the variation in density will be
        ! low (ie < 10%)
        have_reference_density=have_option("/material_phase::"&
            //trim(state%name)//"/equation_of_state/fluids/linear/reference_density")

        if(have_reference_density) then
            call get_option("/material_phase::"&
                //trim(state%name)//"/equation_of_state/fluids/linear/reference_density", &
                rho)
        else
            FLAbort("DG_LES: missing reference density option in equation_of_state")
        end if


        ! Calculate the Poincare constant from the Smagorinsky coefficient
        call get_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/" &
            // "smagorinsky_coefficient", &
            Cs)
        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes


        ! Only calculate new lengths if mesh geometry is new (eg. after adapted mesh)
        if(new_mesh_geometry .or. smooth_called_count<2) then
            ! Calculate nodal filter lengths.

            do e=1, num_elements
                elelen%val(:, e) = element_volume(x, e)**(1./3.)
            end do
            call halo_update(elelen)

            ! project to node-based (1st order) field
            call project_field(elelen, nodelen, x)
            call halo_update(nodelen)

            ! Final lengthscale calculation - smooth
            call anisotropic_smooth_vector(nodelen, x, smoothlen, 3.5,  &
                trim(complete_field_path(smoothlen%option_path)) // "/algorithm")
            call halo_update(smoothlen)

            smooth_called_count = smooth_called_count+1
        end if

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

        ! calculate it at each node
        do n=1, num_nodes

           Cs_length_sq = (Cs* smoothlen%val(1,n))**2.0
           rate_of_strain = 0.5 * (u_grad%val(:,:, n) + transpose(u_grad%val(:,:, n)))

           visc_turb = Cs_length_sq * rho * norm2(2.0 * rate_of_strain)
           
           ! We may add terms to the viscosity here, so...
           tmp_visc = visc_turb

           if(have_van_driest) then
              
              u_grad_node = u_grad%val(:,:, n)
              y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
              vd_damping = (1-exp(-y_plus/A_plus))**pow_m
              tmp_visc = vd_damping * visc_turb
           else
              tmp_visc = visc_turb
           end if
           
           if(have_artificial_visc) then
              tmp_visc = tmp_visc + artificial_visc%val(n)
           end if

           call set(sgs_visc, n,  tmp_visc)

        end do


        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)


    end subroutine calc_dg_sgs_scalar_viscosity



    ! ========================================================================
    ! Scalar LES SGS viscosity
    ! ========================================================================
    subroutine calc_dg_sgs_scalar_viscosity_old(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc, dist_to_wall

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field) :: u_grad
        type(tensor_field), pointer :: mviscosity

        integer :: e, num_elements, n, num_nodes, ln

        real :: Cs, length, ele_vol, Cs_length_sq
        real, dimension(u%dim, u%dim):: rate_of_strain, u_grad_node
        real :: sgs_ele_av, visc_turb, mu, node_visc
        real :: rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density, use_dg_velocity
        type(vector_field), pointer :: elelen, nodelen

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0

        print*, "calc_dg_sgs_scalar_viscosity (Smag)"

        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")

        ! Using DG velocity field for LES calculations.
        ! This is SLOW! Only use for performance tests.
        use_dg_velocity = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/use_dg_velocity")

        ! Length scales for filter
        elelen => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodelen => extract_vector_field(state, "NodeLengthScales", stat=state_flag)

        if(state_flag == 0) then
           print*, "**** Have ScalarElementLengthScales and ScalarNodeLengthScales fields"
        else
           FLAbort("Error: must have ElementLengthsScales and NodeLengthScales fields for Scalar DG LES. This should have been automatically created.")
        end if


        ! Can either use projected CG or DG velocity field in LES calculations
        ! CG is preferred as it is much, much faster.
        ! Velocity projected to continuous Galerkin

        if(use_dg_velocity) then
            u_cg=>extract_vector_field(state, "Velocity", stat=state_flag)
        else
            u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)
        end if

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")
        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)

        ! Crucially, update halos for use
        call halo_update(u_grad)

        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        ! (This will be checked for elsewhere)

        mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

        if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
            .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
            mu=mviscosity%val(1,1,1)
        else
            FLAbort("DG_LES: must have constant dynamic viscosity field")
        end if

        if(have_van_driest) then
            dist_to_wall=>extract_scalar_field(state, "DistanceToWall", stat=state_flag)
            if (state_flag/=0) then
                FLAbort("DG_LES: Van Driest damping requested, but no DistanceToWall scalar field exists")
            end if

        end if

        ! We only use the reference density. This assumes the variation in density will be
        ! low (ie < 10%)
        have_reference_density=have_option("/material_phase::"&
            //trim(state%name)//"/equation_of_state/fluids/linear/reference_density")

        if(have_reference_density) then
            call get_option("/material_phase::"&
            //trim(state%name)//"/equation_of_state/fluids/linear/reference_density", &
            rho)
        else
            FLAbort("DG_LES: missing reference density option in equation_of_state")
        end if


        call get_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/" &
            // "smagorinsky_coefficient", &
            Cs)


        ! Calculate element-wise filter lengths
        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes
        do e=1, num_elements
            elelen%val(:,e) = (element_volume(x, e))**(1./3.)
        end do
        call halo_update(elelen)

        ! project to node-based (1st order) field
        call project_field(elelen, nodelen, x)
        call halo_update(nodelen)

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0



        ! Go round all the nodes, calculating SGS viscosity
        num_nodes = u_cg%mesh%nodes
        do n=1, num_nodes
            Cs_length_sq = (Cs* nodelen%val(1,n))**2.0
            rate_of_strain = 0.5 * (u_grad%val(:,:, n) + transpose(u_grad%val(:,:, n)))
            visc_turb = Cs_length_sq * rho * norm2(2.0 * rate_of_strain)

            if(have_van_driest) then

                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m
                node_visc = vd_damping * visc_turb
            else
                node_visc = visc_turb
            end if

            call set(sgs_visc, n, node_visc )
        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)


    end subroutine calc_dg_sgs_scalar_viscosity_old


    ! =========================================================================
    !  Vreman Model.
    ! =========================================================================

    subroutine calc_dg_sgs_vreman_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        type(scalar_field), pointer :: dist_to_wall
        type(scalar_field), pointer :: sgs_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(scalar_field), pointer :: artificial_visc

        type(vector_field), pointer :: elelen, nodelen, smoothlen
        type(tensor_field) :: u_grad
        real, allocatable :: dx_ele_raw(:,:)
        
        integer :: i, j, k, m, e, num_elements, n, num_nodes
        integer :: u_cg_ele(ele_loc(u,1))

        ! Vreman specific
        ! real, allocatable :: alpha(:,:), beta(:,:), delta(:)
        real :: Cs, Cpoin !, B_beta, alpha_sq_sum, beta_sum

        ! From Vreman source
!        real :: d1,d2,d3
!        real :: d1v1,d2v1,d3v1,d1v2,d2v2,d3v2,d1v3,d2v3,d3v3
!        real :: b11,b12,b13,b22,b23,b33
        real :: abeta, bbeta
        real, allocatable :: dudx(:,:), beta_ij(:,:), dx(:)

        ! Standard LES vars
        real :: visc_turb, tmp_visc, tmp_val, sgs_max
        real :: mu, rho
        integer :: state_flag

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_reference_density, have_lengths_field

        real :: sgs_visc_val
        logical :: have_artificial_visc, have_artificial_scaling
        real :: max_artificial_visc, av_alpha

        ! Wall stuff
        logical :: have_wall
!        real :: wdamp=1.0
        
        print*, "In calc_dg_sgs_vreman_viscosity()"

        ! allocate( alpha(opDim,opDim), beta(opDim,opDim), delta(opDim) )
        
        t1=mpi_wtime()

!        nullify(elelen)
!        nullify(nodelen)
        nullify(dist_to_wall)
        nullify(mviscosity)

        allocate( dudx(opDim, opDim), beta_ij(opDim, opDim), dx(opDim) )


        ! Velocity projected to continuous Galerkin
!        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)
        u_cg=>extract_vector_field(state, "ProjectedNonlinearVelocity", stat=state_flag)

        ! Allocate gradient field and calculate gradient
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")
        call grad(u_cg, x, u_grad)
        call halo_update(u_grad)

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for Vreman DG LES. (This should not happen)")
        end if


        ! Length scales for filter
        elelen => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodelen => extract_vector_field(state, "NodeLengthScales", stat=state_flag)
        smoothlen => extract_vector_field(state, "SmoothedLengthScales", stat=state_flag)

        if(state_flag == 0) then
           print*, "**** Have ElementLengthScales and NodeLengthScales fields"
           have_lengths_field = .true.
        else
           FLAbort("Error: must have ElementLengthsScales and NodeLengthScales fields for Vreman DG LES. This should have been automatically created.")
        end if

        ! Distance to wall used to shut off Vreman at wall
        dist_to_wall=>extract_scalar_field(state, "DistanceToWall", stat=state_flag)
        ! If there is no DistanceToWall field
        if (state_flag/=0) then
           have_wall=.false.
        else
           have_wall=.true.
        end if

        
        ! We can use this in areas of insufficient resolution
        have_artificial_visc = .false.
        artificial_visc => extract_scalar_field(state, "ArtificialViscosity", &
             stat=state_flag)

        max_artificial_visc = maxval(artificial_visc%val(:))

        if(state_flag == 0) then
           print*, "ArtificialViscosity field detected."

           have_artificial_visc = .true.
            have_artificial_scaling=have_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/" &
                // "scale_with_artificial_viscosity")

            if(have_artificial_scaling) &
                 print*, "Scaling between LES viscosity and ArtificialViscosity field"
        end if


        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

        if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
            .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
            mu=mviscosity%val(1,1,1)
        else
            FLAbort("DG_LES: must have constant or normal viscosity field")
        end if


        ! We only use the reference density. This assumes the variation in density will be
        ! low (ie < 10%)
        have_reference_density=have_option("/material_phase::"&
            //trim(state%name)//"/equation_of_state/fluids/linear/reference_density")

        if(have_reference_density) then
            call get_option("/material_phase::"&
                //trim(state%name)//"/equation_of_state/fluids/linear/reference_density", &
                rho)
        else
            FLAbort("DG_LES: missing reference density option in equation_of_state")
        end if


        ! Calculate the Poincare constant from the Smagorinsky coefficient
        call get_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/" &
            // "smagorinsky_coefficient", &
            Cs)
        Cpoin = 2.5 * (Cs**2)
        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        allocate(dx_ele_raw(3, num_elements))

        ! Only calculate new lengths if mesh geometry is new (eg. after adapted mesh)
        if(new_mesh_geometry .or. smooth_called_count<2) then
            ! Calculate nodal filter lengths.
            call aniso_filter_elelengths(x, dx_ele_raw)
            ! Put into field

            do e=1, num_elements
                elelen%val(:, e) = dx_ele_raw(:, e)
            end do
            call halo_update(elelen)

            ! project to node-based (1st order) field
            call project_field(elelen, nodelen, x)
            call halo_update(nodelen)

            ! Final lengthscale calculation - smooth
            call anisotropic_smooth_vector(nodelen, x, smoothlen, 3.5,  &
                trim(complete_field_path(smoothlen%option_path)) // "/algorithm")
            call halo_update(smoothlen)

            smooth_called_count = smooth_called_count+1
        end if

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

        ! calculate it at each node
        do n=1, num_nodes

           ! Wall stuff, if needed. This shuts off Vreman at wall boundary
!           if(have_wall) then
!              if(abs(dist_to_wall%val(n)) < 10e-10 ) then
!                 wdamp=0.0
!              else
!                 wdamp=1.0
!              end if
!           end if

!           wdamp=1.0
           
!           d1 = smoothlen%val(1, n)**2
!           d2 = smoothlen%val(2, n)**2
!           d3 = smoothlen%val(3, n)**2

!           Not entirely sure the *actual* Vreman code is any good. Commented
!           out for now
!
!           d1v1=u_grad%val(1,1,n) ! du1dx1(n)
!           d2v1=u_grad%val(1,2,n) ! du1dx2(n)
!           d3v1=u_grad%val(1,3,n) ! du1dx3(n)
!
!           d1v2=u_grad%val(2,1,n) ! du2dx1(n)
!           d2v2=u_grad%val(2,2,n) ! du2dx2(n)
!           d3v2=u_grad%val(2,3,n) ! du2dx3(n)
!
!           d1v3=u_grad%val(3,1,n) ! du3dx1(n)
!           d2v3=u_grad%val(3,2,n) ! du3dx2(n)
!           d3v3=u_grad%val(3,3,n) ! du3dx3(n)
!
!           b11=d1*d1v1*d1v1+d2*d2v1*d2v1+d3*d3v1*d3v1
!           b12=d1*d1v1*d1v2+d2*d2v1*d2v2+d3*d3v1*d3v2
!           b13=d1*d1v1*d1v3+d2*d2v1*d2v3+d3*d3v1*d3v3
!           b22=d1*d1v2*d1v2+d2*d2v2*d2v2+d3*d3v2*d3v2
!           b23=d1*d1v2*d1v3+d2*d2v2*d2v3+d3*d3v2*d3v3
!           b33=d1*d1v3*d1v3+d2*d2v3*d2v3+d3*d3v3*d3v3
!
!           abeta = d1v1**2 + d1v2**2 + d1v3**2 &
!                + d2v1**2 + d2v2**2 + d2v3**2 &
!                + d3v1**2 + d3v2**2 + d3v3**2
!
!           bbeta=b11*b22-(b12**2)+b11*b33-(b13**2)+b22*b33-(b23**2)

           dudx(:,1) = u_grad%val(1,:,n)
           dudx(:,2) = u_grad%val(2,:,n)
           dudx(:,3) = u_grad%val(3,:,n)

           dx = smoothlen%val(:,n)

           abeta = sum(u_grad%val(:,:,n)*u_grad%val(:,:,n))
           beta_ij = 0.0

           do i=1, opDim
              do j=1, opDim
                 do k=1, opDim
                    beta_ij(i,j) = beta_ij(i,j) + (dx(k)**2.)*dudx(k,i)*dudx(k,j)
                 end do
              end do
           end do
           
           bbeta = &
                  (beta_ij(1,1)*beta_ij(2,2)) &
                  + (beta_ij(1,1)*beta_ij(3,3)) &
                  + (beta_ij(2,2)*beta_ij(3,3)) &
                  - (beta_ij(1,2)**2.)-(beta_ij(1,3)**2.)-(beta_ij(2,3)**2.)


           ! If abeta is v. small, sgs turb visc must be zero also.
           if(abeta < 10e-10) then
              visc_turb = 0.0
           else 
              ! Otherwise calculate Vreman visc_turb.
              
              ! Limiter on B_beta
              if(bbeta<0) bbeta=0.0
              
              visc_turb = Cpoin * sqrt( bbeta / abeta )
           end if

           ! Limit on sgs viscosity
           sgs_max = min( Cpoin * sqrt(abeta/3) * max(dx(1), dx(2), dx(3)),  mu*10e4 )
           if(visc_turb > sgs_max) visc_turb=sgs_max
           
           ! Artificial viscosity nonsense
           if(have_artificial_visc) then
              if(.not. have_artificial_scaling) then
                  tmp_visc = visc_turb + artificial_visc%val(n)

!                  if(tmp_visc>max_artificial_visc) tmp_visc=max_artificial_visc
              else
                  ! Otherwise we scale between one and the other
                  av_alpha = artificial_visc%val(n)/max_artificial_visc
                  tmp_visc = (1-av_alpha)*visc_turb + av_alpha*artificial_visc%val(n)
               end if
           end if
           
           call set(sgs_visc, n,  tmp_visc)
           
        end do


        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)
        deallocate( dudx, beta_ij, dx )

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)


    end subroutine calc_dg_sgs_vreman_viscosity



    ! ========================================================================
    ! Use LES for anisotropic grids, using vorticity-derived filter lengths.
    ! Based upon Chauvet (2007)
    ! ========================================================================

    subroutine calc_dg_sgs_chauvet_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        FLAbort("DG Chauvet LES is no longer supported.")

    end subroutine calc_dg_sgs_chauvet_viscosity



    subroutine les_length_scales_squared_mk2(positions, ele, horzSq, vertSq)
        type(vector_field), intent(in) :: positions
        integer :: ele

        real :: horzSq, vertSq
        real, dimension(:,:), allocatable :: X_val
        real, dimension(:), allocatable :: dx

        integer :: i

        allocate( X_val(positions%dim, opNloc) )
        allocate( dx(positions%dim) )

        X_val=ele_val(positions, ele)

        ! Calculate largest dx, dy, dz for element nodes
        do i=1, opDim
            dx(i) = maxval( X_val(i,:)) - minval (X_val(i,:))
        end do

        ! Why squares? Used in LES visc calcs, no need for expensive square roots
        horzSq = dx(1)**2 + dx(2)**2
        vertSq = dx(3)**2

        deallocate(X_val, dx)

    end subroutine les_length_scales_squared_mk2



    ! ========================================================================
    ! Calculate vorticity-based filter lengths based upon Chauvet et al.
    ! 3D only.
    ! ========================================================================

    subroutine vorticity_filter_lengths(pos, ugrad, dx, vort_del)
      type(vector_field), intent(in) :: pos
      type(tensor_field), intent(in) :: ugrad
      real, dimension(:,:), intent(in) :: dx
      real, dimension(:), intent(out) :: vort_del

      real :: magvort
      real, dimension(:), allocatable :: vort, N
      real, dimension(:,:), allocatable :: gr

      integer :: i, num_nodes

      num_nodes = pos%mesh%nodes

      allocate( vort(pos%dim), N(pos%dim) )
      allocate( gr(pos%dim, pos%dim) )

      ! Go round all nodes, calculating vorticity-based filter length
      do i=1, num_nodes
         gr=ugrad%val(:,:,i)

         ! Calculate vorticity: \/ x u
         vort(1) = gr(2,3)-gr(3,2)
         vort(2) = gr(3,1)-gr(3,1)
         vort(3) = gr(2,1)-gr(1,2)

         magvort = sqrt(vort(1)**2 + vort(2)**2 + vort(3)**2)

         ! Chauvet et al.
         N = vort / magvort
         vort_del(i) = sqrt( dx(2, i) * dx(3,i) * N(1)**2 &
              + dx(3, i) * dx(1,i) * N(2)**2 &
              + dx(1, i) * dx(2,i) * N(3)**2 )

      end do

      deallocate(gr, vort, N)

    end subroutine vorticity_filter_lengths


    ! -----------------------------------------------------------------------
    ! Roman et al tensor viscosity. Split horizontal/vertical LES SGS
    ! viscosity. Requres tensor sgs viscosity field.
    ! ========================================================================

   subroutine calc_dg_sgs_roman_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        FLAbort("DG Roman LES is no longer supported.")

    end subroutine calc_dg_sgs_roman_viscosity


    ! ========================================================================
    ! Aniso lengthscales, but really only for pancake elements (small dz).
    ! ========================================================================

    subroutine aniso_length_ele(pos, ele, del, extruded_mesh)
        type(vector_field), intent(in) :: pos
        integer :: ele
        real, intent(inout) :: del(:)
        real :: dx(pos%dim)
        
        real :: X_val(pos%dim, opNloc)
        real :: X_mean(pos%dim), r, diffx, diffy, diffz, maxdz
        real :: X_tri(pos%dim-1, 3)
        real :: area, a, b, c, s, tmplen, vol
        real, parameter :: pi=3.14159265359
        integer :: i, n, m, trix
        integer :: stpair(2)

        ! Conditional to switch in/out specialised extruded mesh case logic
        logical :: extruded_mesh

        X_val=ele_val(pos, ele)
        vol = element_volume(pos, ele)

        ! What we do depends upon whether this is an extruded mesh or not.
!
!        if(extruded_mesh) then
!            ! First look for two points vertically aligned. These two will provide
!            ! dz metric. (There are always two in an extruded mesh)
!
!            stpair=0.0
!            maxdz=0.0
!
!            maxdz=abs(maxval(X_val(3,:))-minval(X_val(3,:)))
!
!            outer_loop: do n=1, opNloc
!
!                do m=1, opNloc
!                    if ( n /= m ) then
!                        diffx = abs(X_val(1, n)-X_val(1,m))
!                        diffy = abs(X_val(2, n)-X_val(2,m))
!
!                        ! If two points share x and y, one must be atop the other.
!                        if(diffx < 10e-10 .and. diffy < 10e-10) then
!                           stpair(1)=n
!                            stpair(2)=m
!                            exit outer_loop
!                        end if
!                    end if
!                end do
!
!            end do outer_loop
!            del(3) = abs( X_val(3,stpair(1))-X_val(3,stpair(2)) )
!
!            ! Catch all. If somehow dz left zero, use this instead.
!            if(del(3) < 10e-10) del(3) = maxdz
!
!            ! Find 2D points for horizontal triangle (easier/quicker than using
!            ! Fluidity framework)
!            trix=1
!            do n=1, opNloc
!               if(n /= stpair(1)) then
!                    X_tri(:, trix)=X_val(1:2, n)
!                    trix=trix+1
!                end if
!            end do
!
!            ! Heron's formula for area of triangle
!            a = sqrt( (X_tri(1,2)-X_tri(1,1))**2. + (X_tri(2,2)-X_tri(2,1))**2.)
!            b = sqrt( (X_tri(1,3)-X_tri(1,2))**2. + (X_tri(2,3)-X_tri(2,2))**2.)
!            c = sqrt( (X_tri(1,3)-X_tri(1,1))**2. + (X_tri(2,3)-X_tri(2,1))**2.)
!
!            s = 0.5*(a+b+c)
!
!            area = (s*(s-a)*(s-b)*(s-c))**0.5
!
!            if(area>10e5) print*, "WARNING: AMD LES metrics: large area"
!
!            ! Calculate radius of circle with same area
!            r = sqrt(area/ 3.141592653)
!
!            del(1) = 2.*r
!            del(2) = del(1)
!
!            ! Not quite done. Sanity check for very, very thin elements
!            if(del(3)/del(1) < 0.05) del(3)=0.05 * del(1)
!
!        else

       ! For unstructured meshes, do something slightly different.
        del(:)=0.0
        do m=1, opNloc
            do n=1, opNloc

                if (n/=m) then
                    diffx = abs(X_val(1, n)-X_val(1,m))
                    diffy = abs(X_val(2, n)-X_val(2,m))
                    diffz = abs(X_val(3, n)-X_val(3,m))

                    if(diffx > del(1)) del(1)=diffx
                    if(diffy > del(2)) del(2)=diffy
                    if(diffz > del(3)) del(3)=diffz
                end if
            end do
        end do

        if(extruded_mesh) then
           ! tmplen = sqrt(del(1)**2 + del(2)**2)
           ! del(1) = tmplen
           ! del(2) = tmplen

           ! Fitting to an ellipsoid
           del(1) = 2.0 * ( (6*vol) / (pi*del(3)) )**0.5
           del(2) = del(1)
        end if

    end subroutine aniso_length_ele



    ! This one is not per-element as above, but loops over all elements
    subroutine aniso_filter_lengths(pos, del_fin, lenrange)
        type(vector_field), intent(in) :: pos
        real, dimension(:,:), intent(inout) :: del_fin

        real, optional :: lenrange(2)
        real :: tmplen

        real, dimension(pos%dim, opNloc) :: X_val
        real, dimension(pos%dim) :: dx, del

        real, allocatable :: dx_sum(:,:), dx_ele_raw(:,:) !, dx_ele_filt(:,:)
        real, dimension(pos%dim) :: dx_neigh_sum, dx_neigh_average
        real :: del_max, mean_val, ele_filt_sum(opDim)
        integer, allocatable :: visits(:)
        integer, pointer, dimension(:) :: neighs
        integer :: local_gnodes(ele_loc(pos,1))

        integer :: e, n, i, gn, f
        integer :: num_elements, num_nodes, num_neighs

        integer :: ierr

        logical :: extruded_mesh

        num_elements = ele_count(pos)
        num_nodes = node_count(pos)

        allocate(dx_sum(pos%dim, num_nodes))
        allocate(visits(num_nodes))
        allocate(dx_ele_raw(pos%dim, num_elements))
        ! allocate(dx_ele_filt(pos%dim, num_elements))

        ! Reset counters and sums
        visits=0
        dx_sum=0.

        ! Is this an extruded mesh?
        extruded_mesh = option_count("/geometry/mesh/from_mesh/extrude") > 0

!        if( .not. extruded_mesh ) then
!            FLExit("Error: AMD LES is currently only works correctly on extruded meshes")
!        end if

        ! Raw sizes per element
        do e=1, num_elements
           if(element_owned(pos, e)) then
               call aniso_length_ele(pos, e, dx, extruded_mesh=extruded_mesh)

               dx_ele_raw(:,e) = dx(:)
            end if
        end do

!        ! Filter sizes per element
!         do e=1, num_elements
!             neighs=>ele_neigh(pos, e)
!             ele_filt_sum=dx_ele_raw(:,e)
!             num_neighs=0
!             do f=1, size(neighs)
!                 if(neighs(f)>0) then
!                     ele_filt_sum=ele_filt_sum+dx_ele_raw(:,neighs(f))
!                     num_neighs=num_neighs+1
!                 end if
!             end do
!             dx_ele_filt(:,e) = ele_filt_sum / num_neighs
!         end do


        ! Now create sizes per-node
        do e=1, num_elements
            if(element_owned(pos, e)) then
                local_gnodes = ele_nodes(pos, e)

                X_val=ele_val(pos, e)

                do n=1, opNloc
                   gn=local_gnodes(n)

                   if(gn>0 .and. gn<=num_nodes) then

                       ! Add element sizes to node size sum at each corner
                       dx_sum(:, gn) = dx_sum(:, gn) + dx_ele_raw(:,e)
                       visits(gn) = visits(gn)+1
                   else
                       print*, "gn:", gn
                   end if
                end do
            end if
        end do

        ! Now set node values, and calculate scalar min/max range
        lenrange(1) = 10e10
        lenrange(2) = 0.
        do n=1, num_nodes
            del_fin(:, n) = dx_sum(:, n) / visits(n)

            if(present(lenrange)) then
                tmplen = maxval(del_fin(:,n))

                if(tmplen<lenrange(1)) lenrange(1)=tmplen
                if(tmplen>lenrange(2)) lenrange(2)=tmplen
            end if
            print*, "del_fin: ", del_fin(:,n)
        end do

        ! Only output element length stats if we've asked for them
        if(present(lenrange)) then
            ! Calculate global min/max values
            call mpi_allreduce(lenrange(1), tmplen, 1, mpi_double_precision, &
                mpi_min, mpi_comm_world, ierr)
            lenrange(1) = tmplen

            call mpi_allreduce(lenrange(2), tmplen, 1, mpi_double_precision, &
                mpi_max, mpi_comm_world, ierr)
            lenrange(2) = tmplen
        end if

        deallocate(dx_sum)
        deallocate(dx_ele_raw)
        ! deallocate(dx_ele_filt)
        deallocate(visits)

    end subroutine aniso_filter_lengths





    ! This one is not per-element as above, but loops over all elements
    subroutine aniso_filter_elelengths(pos, dx_ele_raw)
        type(vector_field), intent(in) :: pos
        real, dimension(:,:), intent(inout) :: dx_ele_raw
        real :: tmplen

        real, dimension(pos%dim) :: dx

        integer :: e
        integer :: num_elements

        logical :: extruded_mesh

        num_elements = ele_count(pos)

        ! allocate(dx_ele_filt(pos%dim, num_elements))


        ! Is this an extruded mesh?
        extruded_mesh = option_count("/geometry/mesh/from_mesh/extrude") > 0

        ! Raw sizes per element
        do e=1, num_elements
           if(element_owned(pos, e)) then
               call aniso_length_ele(pos, e, dx, extruded_mesh=extruded_mesh)

               dx_ele_raw(:,e) = dx(:)
            end if
        end do

    end subroutine aniso_filter_elelengths




    ! This one is not per-element as above, but loops over all elements
    subroutine aniso_filter_lengths_field(state, pos)
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: pos

        real, dimension(pos%dim) :: dx
        type(vector_field), pointer :: elef, nodef
        integer :: e, num_elements, state_flag
        logical :: extruded_mesh

        num_elements = ele_count(pos)

        ! Is this an extruded mesh?
        extruded_mesh = option_count("/geometry/mesh/from_mesh/extrude") > 0

!        if( .not. extruded_mesh ) then
!            FLExit("Error: AMD LES is currently only works correctly on extruded meshes")
!        end if

        ! Length scales for filter
        nullify(elef)
        nullify(nodef)

        elef => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodef => extract_vector_field(state, "NodeLengthScales", stat=state_flag)


        ! Set element dimensions in ElementLengthscales field. (0th order)
        do e=1, num_elements
           call aniso_length_ele(pos, e, dx, extruded_mesh=extruded_mesh)

           elef%val(1, e) = dx(1)
           elef%val(2, e) = dx(2)
           elef%val(3, e) = dx(3)
        end do
        call halo_update(elef)

        call quick_smooth_ele_lengths_field(state, 10)

        ! Now project to NodeLengthscales field. (1st order)
        call project_field(elef, nodef, pos)
        call halo_update(nodef)

    end subroutine aniso_filter_lengths_field




    ! Very quick / dirty smoothing routine with averaging filter for
    ! element lenegths field

    subroutine quick_smooth_ele_lengths_field(state, niters)
        type(state_type), intent(in) :: state
        integer, intent(in) :: niters

        type(vector_field), pointer :: elef

        integer :: nix
        integer :: e, i, nh, num_elements, num_neighs
        real, allocatable :: ele_len_sum(:,:)
        integer, allocatable :: ele_sum_ct(:)

        integer, pointer, dimension(:) :: neighs
        integer :: state_flag

        nullify(elef)
        elef=>extract_vector_field(state, "ElementLengthScales", stat=state_flag)


        num_elements = ele_count(elef)


        allocate(ele_len_sum(elef%dim, num_elements))
        allocate(ele_sum_ct(num_elements))

        ! Smooth field using averaging filter
        do nix=1, niters
            do e=1, num_elements
                neighs => ele_neigh(elef, e)
                num_neighs = size(neighs)

                ele_sum_ct(e) = 1+num_neighs

                do i=1, 3
                    ele_len_sum(i, e) = elef%val(i, e)
                end do

                ! Doing an average of an elements lengthscales and its neighbours' - sum the lot
                do nh=1, num_neighs
                    do i=1, 3
                        ele_len_sum(i, e) = ele_len_sum(i, e) + elef%val(i, neighs(nh))
                    end do
                end do
            end do

            ! Now calculate the average of those summed lengths, and put in elef field.
            do e=1, num_elements
                do i=1, 3
                    elef%val(i, e) = ele_len_sum(i, e) / ele_sum_ct(e)
                end do
            end do

            call halo_update(elef)
        end do

        deallocate(ele_len_sum)
        deallocate(ele_sum_ct)

    end subroutine quick_smooth_ele_lengths_field



end module dg_les

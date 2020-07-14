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
  use vector_tools
  use state_fields_module
  use solvers
  use global_parameters
  use spud
  use halos

  implicit none

  private

  public :: calc_dg_sgs_scalar_viscosity
  public :: calc_dg_sgs_tensor_viscosity, calc_dg_sgs_amd_viscosity
  public :: calc_dg_sgs_chauvet_viscosity

  ! This scales the Van Driest effect
  real, parameter :: van_scale=1.0

#include "mpif.h"


contains


    ! ========================================================================
    ! Scalar LES SGS viscosity (dynamic)
    ! ========================================================================
    subroutine calc_dg_sgs_scalar_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc, dist_to_wall

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field) :: u_grad
        type(tensor_field), pointer :: mviscosity

        integer :: e, num_elements, n, num_nodes, ln
        integer :: u_cg_ele(ele_loc(u,1))

        real :: Cs, length, ele_vol, Cs_length_sq
        real, dimension(u%dim, u%dim) :: rate_of_strain, u_grad_node
        real :: sgs_ele_av, visc_turb, mu, node_visc
        real :: rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable, save :: node_vol_weighted_sum(:), &
             node_neigh_total_vol(:)
        
        real, allocatable, save :: node_sum(:)
        integer, allocatable, save :: node_visits(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density, use_dg_velocity

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0

        print*, "calc_dg_sgs_scalar_viscosity"

        t1=mpi_wtime()


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

!        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
!        ! (This will be checked for elsewhere)

        if(have_van_driest) then
            mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

            if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
                .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
                mu=mviscosity%val(1,1,1)
            else
                FLAbort("DG_LES: must have constant dynamic viscosity field")
            end if

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

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_sum)) then
                deallocate(node_sum)
                deallocate(node_visits)
                
                deallocate(node_vol_weighted_sum)
                deallocate(node_neigh_total_vol)
            end if

            allocate(node_sum(num_nodes))
            allocate(node_visits(num_nodes))

            allocate(node_vol_weighted_sum(num_nodes))
            allocate(node_neigh_total_vol(num_nodes))
            
        end if

        node_sum(:)=0.0
        node_visits(:)=0
        node_vol_weighted_sum(:)=0.0
        node_neigh_total_vol(:)=0.0
        
        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0


        do e=1, num_elements
            u_cg_ele=ele_nodes(u_cg, e)

            ele_vol = element_volume(x, e)

            length=(ele_vol**0.333333333333333333333)
            ! Factor of two included in length_scale_scalar
            Cs_length_sq = (Cs* length)**2.0

            ! This is the contribution to nu_sgs from each co-occupying node
            sgs_ele_av=0.0
            do ln=1, opNloc
                gnode = u_cg_ele(ln)
                rate_of_strain = 0.5 * (u_grad%val(:,:, gnode) + transpose(u_grad%val(:,:, gnode)))
                visc_turb = Cs_length_sq * rho * norm2(2.0 * rate_of_strain)

                sgs_ele_av = sgs_ele_av + visc_turb/opNloc

                node_sum(gnode) = node_sum(gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
            end do

        end do


        ! Set final values. Two options here: one with Van Driest damping, 
        ! one without.
        if(have_van_driest) then
            do n=1, num_nodes
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m

                node_visc =  vd_damping * rho*node_sum(n) / node_visits(n)

                call set(sgs_visc, n, node_visc)
            end do
        else
            do n=1, num_nodes
               node_visc = rho*node_sum(n) / node_visits(n)

               call set(sgs_visc, n, node_visc )
            end do
        end if

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)


    end subroutine calc_dg_sgs_scalar_viscosity


    ! Just a wrapper for switching between LES tensor subroutines.

    subroutine calc_dg_sgs_tensor_viscosity(state, x, u)
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        ! call calc_dg_sgs_tensor_viscosity_original(state, x, u)
        
        call calc_dg_sgs_tensor_viscosity_node(state, x, u)

    end subroutine calc_dg_sgs_tensor_viscosity


    ! ========================================================================
    ! Split horizontal/vertical LES SGS viscosity
    ! ========================================================================

    subroutine calc_dg_sgs_tensor_viscosity_original(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: dist_to_wall
        type(tensor_field), pointer :: sgs_visc
        type(scalar_field), pointer :: sgs_visc_mag

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes, ln
        integer :: u_cg_ele(ele_loc(u,1))

        real :: Cs_horz, Cs_length_horz_sq, Cs_vert, Cs_length_vert_sq
        real :: length_horz_sq, length_vert_sq, ele_vol
        real :: mag_strain_horz, mag_strain_vert, mag_strain_r
        real :: sgs_horz, sgs_vert, sgs_r
        real :: mu, rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable,save:: node_vol_weighted_sum(:,:,:),node_neigh_total_vol(:)
        real, allocatable,save:: node_sum(:, :,:), sgs_unfiltered_sum(:,:,:)
        integer, allocatable, save :: node_visits(:)

        real :: visc_norm2, visc_norm2_max

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0

        ! Control of mixture of sharp / averaged filter
        real :: blend

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val
        integer :: i

        real, dimension(:,:), allocatable :: u_grad_node, rate_of_strain
        real, dimension(:,:), allocatable :: sgs_ele_av, visc_turb, tmp_tensor

!        real, dimension(u%dim, u%dim) :: u_grad_node, rate_of_strain
!        real, dimension(u%dim, u%dim) :: sgs_ele_av, visc_turb, tmp_tensor

        
        print*, "In calc_dg_sgs_tensor_viscosity()"

        print*, "- using original LES filter"


        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Allocating work arrays
        allocate(u_grad_node(u%dim,u%dim), rate_of_strain(u%dim,u%dim))
        allocate(sgs_ele_av(u%dim,u%dim), visc_turb(u%dim,u%dim))
        allocate(tmp_tensor(u%dim,u%dim))
      
        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)
        call halo_update(u_grad)

        sgs_visc_mag => extract_scalar_field(state, "TensorEddyViscosityMagnitude", stat=state_flag)

        if (state_flag /= 0) then
            FLAbort("DG_LES: TensorEddyViscosityMagnitude absent for tensor DG LES. (This should not happen)")
         end if

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")

        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        ! (This will be checked for elsewhere)
        if(have_van_driest) then
            mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

            if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
                .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
                mu=mviscosity%val(1,1,1)
            else
                FLAbort("DG_LES: must have constant or normal viscosity field")
            end if

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
            // "smagorinsky_coefficient", Cs_horz)

        if(have_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/" &
            // "smagorinsky_coefficient_vertical")) then

            call get_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/" &
                // "smagorinsky_coefficient_vertical", Cs_vert)
        else
            FLAbort("DG_LES: you've requested anisotropic LES, but have not specified smagorinsky_coefficient_vertical")
        end if

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_sum)) then
                deallocate(node_sum)
                deallocate(node_visits)
                deallocate(node_vol_weighted_sum)
                deallocate(node_neigh_total_vol)
                deallocate(sgs_unfiltered_sum)
            end if

            allocate(node_sum(u%dim, u%dim, num_nodes))
            allocate(node_visits(num_nodes))
            allocate(node_vol_weighted_sum(u%dim, u%dim, num_nodes))
            allocate(node_neigh_total_vol(num_nodes))
            allocate(sgs_unfiltered_sum(u%dim, u%dim, num_nodes))
        end if

        node_sum(:,:,:)=0.0
        sgs_unfiltered_sum(:,:,:)=0.0
        node_visits(:)=0
        node_vol_weighted_sum(:,:,:)=0.0
        node_neigh_total_vol(:)=0.0

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:,:,:)=0.0

        ! Loop over elements
        do e=1, num_elements
            u_cg_ele=ele_nodes(u_cg, e)
            ele_vol = element_volume(x, e)

            call les_length_scales_squared_mk2(x, e, length_horz_sq, length_vert_sq)

            Cs_length_horz_sq = (Cs_horz**2.0)* length_horz_sq
            Cs_length_vert_sq = (Cs_vert**2.0) * length_vert_sq

            ! This is the contribution to nu_sgs from each co-occupying node
            sgs_ele_av=0.0
            do ln=1, opNloc
                gnode = u_cg_ele(ln)
                rate_of_strain = 0.5 * (u_grad%val(:,:, gnode) + transpose(u_grad%val(:,:, gnode)))

                mag_strain_horz = sqrt(2.0* rate_of_strain(1,1)**2.0 &
                                     + 2.0* rate_of_strain(2,2)**2.0 &
                                     + 4.0* rate_of_strain(1,2)**2.0 )

                mag_strain_vert = sqrt(4.0* rate_of_strain(1,3)**2.0 &
                                 + 2.0* rate_of_strain(3,3)**2.0 &
                                 + 4.0* rate_of_strain(3,2)**2.0 )

                mag_strain_r = sqrt(2.0 * rate_of_strain(3,3)**2.0)

                ! Note, this is without density. That comes later.
                sgs_horz = rho * Cs_length_horz_sq * mag_strain_horz
                sgs_vert = rho * Cs_length_vert_sq * mag_strain_vert
                sgs_r = rho * Cs_length_vert_sq * mag_strain_r

                ! As per Roman et al, 2010.
                visc_turb(1:2, 1:2) = sgs_horz
                visc_turb(1, 3) = sgs_vert
                visc_turb(3, 1) = sgs_vert
                visc_turb(2, 3) = sgs_vert
                visc_turb(3, 2) = sgs_vert
                ! visc_turb(3, 3) = sgs_horz + (-2.*sgs_vert) + 2.*sgs_r)
                ! Otherwise this is potentially negative (!)
                visc_turb(3, 3) = sgs_horz + 2.*sgs_vert + 2.*sgs_r

                ! Original and box filter needs the element average
                sgs_ele_av = sgs_ele_av + visc_turb
                sgs_unfiltered_sum(:,:,gnode) = sgs_unfiltered_sum(:,:,gnode) + visc_turb

            end do

            sgs_ele_av = sgs_ele_av / opNloc

            ! Extra calculations for original filter
            do ln=1, opNloc
                gnode = u_cg_ele(ln)

                node_sum(:,:, gnode) = node_sum(:,:, gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
                node_vol_weighted_sum(:,:, gnode) = node_vol_weighted_sum(:,:, gnode) + ele_vol*sgs_ele_av
                node_neigh_total_vol(gnode) = node_neigh_total_vol(gnode) + ele_vol

            end do
        end do


        ! Set final values. Two options here: one with Van Driest damping, one without.
        ! We multiply by rho here.

        ! This controls the mixture of unfiltered (sharp) filter, and averaged cell (soft)
        ! filter. 0.0-1.0 (soft-sharp)
        blend=0.5

        if(have_van_driest) then

            do n=1, num_nodes
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m

!                tmp_tensor = vd_damping * rho * node_sum(:,:,n)/node_visits(n)

                tmp_tensor = vd_damping * rho * &
                ( blend * sgs_unfiltered_sum(:,:,n) + (1-blend) * node_sum(:,:,n) ) &
                / node_visits(n)

!                Not using vol-weighted sums for now
!                tmp_tensor = vd_damping * rho * node_vol_weighted_sum(:,:,n) &
!                     / node_neigh_total_vol(n))

                ! L2 Norm of tensor (tensor magnitude)
                sgs_visc_val=0.0
                do i=1, opDim
                    sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
                end do
                sgs_visc_val = sqrt(sgs_visc_val)

                call set(sgs_visc, n,  tmp_tensor)
                call set(sgs_visc_mag, n, sgs_visc_val )
            end do

        else
            do n=1, num_nodes

                tmp_tensor = rho * &
                ( blend * sgs_unfiltered_sum(:,:,n) + (1-blend) * node_sum(:,:,n) ) &
                / node_visits(n)
!!                tmp_tensor = rho * node_vol_weighted_sum(:,:,n) &
!!                     / node_neigh_total_vol(n))
!!                tmp_tensor = rho * node_vol_weighted_sum(:,:,n) &
!!                     / node_neigh_total_vol(n))
!                tmp_tensor = rho * node_sum(:,:,n) / node_visits(n)

                ! L2 Norm of tensor (tensor magnitude)
                sgs_visc_val=0.0
                do i=1, opDim
                    sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
                end do
                sgs_visc_val = sqrt(sgs_visc_val)

                ! tmp_tensor=0.0

                call set(sgs_visc, n, tmp_tensor )
                call set(sgs_visc_mag, n, sgs_visc_val )

            end do
        end if
        ! End of code for original and hat filters

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        call halo_update(sgs_visc_mag)

        call deallocate(u_grad)
        deallocate(u_grad_node,rate_of_strain,sgs_ele_av,visc_turb,tmp_tensor)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_tensor_viscosity_original


    ! This calculates per-node tensor LES, ie. a hat filter.

    subroutine calc_dg_sgs_tensor_viscosity_node(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: dist_to_wall
        type(tensor_field), pointer :: sgs_visc
        type(scalar_field), pointer :: sgs_visc_mag

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes
        integer :: u_cg_ele(ele_loc(u,1))

        real :: Cs_horz, Cs_length_horz_sq, Cs_vert, Cs_length_vert_sq
        real :: length_horz_sq, length_vert_sq, ele_vol
        real :: mag_strain_horz, mag_strain_vert, mag_strain_r
        real :: sgs_horz, sgs_vert, sgs_r
        real :: mu, rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable,save:: node_dx(:,:)

        real :: visc_norm2, visc_norm2_max

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=25.0, pow_m=2.0

        ! Control of mixture of sharp / averaged filter
        real :: blend

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val
        integer :: i

        real, dimension(:,:), allocatable :: u_grad_node, rate_of_strain
        real, dimension(:,:), allocatable :: sgs_ele_av, visc_turb, tmp_tensor

!        real, dimension(u%dim, u%dim) :: u_grad_node, rate_of_strain
!        real, dimension(u%dim, u%dim) :: sgs_ele_av, visc_turb, tmp_tensor


        print*, "In calc_dg_sgs_tensor_viscosity()"

        print*, "- using original LES filter"


        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Allocating work arrays
        allocate(u_grad_node(u%dim,u%dim), rate_of_strain(u%dim,u%dim))
        allocate(sgs_ele_av(u%dim,u%dim), visc_turb(u%dim,u%dim))
        allocate(tmp_tensor(u%dim,u%dim))

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)
        call halo_update(u_grad)

        sgs_visc_mag => extract_scalar_field(state, "TensorEddyViscosityMagnitude", stat=state_flag)

        if (state_flag /= 0) then
            FLAbort("DG_LES: TensorEddyViscosityMagnitude absent for tensor DG LES. (This should not happen)")
         end if

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")


        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        ! (This will be checked for elsewhere)
        if(have_van_driest) then
            mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

            if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
                .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
                mu=mviscosity%val(1,1,1)
            else
                FLAbort("DG_LES: must have constant or normal viscosity field")
            end if

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
            // "smagorinsky_coefficient", Cs_horz)

        if(have_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/" &
            // "smagorinsky_coefficient_vertical")) then

            call get_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/" &
                // "smagorinsky_coefficient_vertical", Cs_vert)
        else
            FLAbort("DG_LES: you've requested anisotropic LES, but have not specified smagorinsky_coefficient_vertical")
        end if

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_dx)) then
                deallocate(node_dx)
            end if

            allocate(node_dx(u%dim, num_nodes))
            call aniso_filter_lengths(x, node_dx)
        end if

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:,:,:)=0.0

        ! Loop over nodes
        do n=1, num_nodes


            Cs_length_horz_sq = (Cs_horz**2.0)* node_dx(1,n)**2.
            Cs_length_vert_sq = (Cs_vert**2.0) * node_dx(3,n)**2.

            gnode = n
            rate_of_strain = 0.5 * (u_grad%val(:,:, gnode) + transpose(u_grad%val(:,:, gnode)))

            mag_strain_horz = sqrt(2.0* rate_of_strain(1,1)**2.0 &
                                 + 2.0* rate_of_strain(2,2)**2.0 &
                                 + 4.0* rate_of_strain(1,2)**2.0 )

            mag_strain_vert = sqrt(4.0* rate_of_strain(1,3)**2.0 &
                             + 2.0* rate_of_strain(3,3)**2.0 &
                             + 4.0* rate_of_strain(3,2)**2.0 )

            mag_strain_r = sqrt(2.0 * rate_of_strain(3,3)**2.0)

            ! Note, this is without density. That comes later.
            sgs_horz = rho * Cs_length_horz_sq * mag_strain_horz
            sgs_vert = rho * Cs_length_vert_sq * mag_strain_vert
            sgs_r = rho * Cs_length_vert_sq * mag_strain_r

            ! As per Roman et al, 2010.
            visc_turb(1:2, 1:2) = sgs_horz
            visc_turb(1, 3) = sgs_vert
            visc_turb(3, 1) = sgs_vert
            visc_turb(2, 3) = sgs_vert
            visc_turb(3, 2) = sgs_vert
            ! visc_turb(3, 3) = sgs_horz + (-2.*sgs_vert) + 2.*sgs_r)
            ! Otherwise this is potentially negative (!)
            !visc_turb(3, 3) = sgs_horz + 2.*sgs_vert + 2.*sgs_r
            visc_turb(3, 3) = abs(sgs_horz + (-2.*sgs_vert) + 2.*sgs_r)

            ! Wall damping with Van Driest
            if(have_van_driest) then
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m

                tmp_tensor = vd_damping * rho * visc_turb

            else
                tmp_tensor = rho * visc_turb

            end if

!            ! For when stuff goes bananas.
!            if(norm2(tmp_tensor) > mu*10e5) tmp_tensor=0.

             ! L2 Norm of tensor (tensor magnitude)
            sgs_visc_val=0.0
            do i=1, opDim
                sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
            end do
            sgs_visc_val = sqrt(sgs_visc_val)

            call set(sgs_visc, n,  tmp_tensor)
            call set(sgs_visc_mag, n, sgs_visc_val )


        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        call halo_update(sgs_visc_mag)

        call deallocate(u_grad)
        deallocate(u_grad_node,rate_of_strain,sgs_ele_av,visc_turb,tmp_tensor)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_tensor_viscosity_node



    ! ========================================================================
    ! Use Anisotropic Minimum Dissipation (AMD) LES filter.
    ! See Rozema et al, Computational Methods in Engineering, 2020.
    ! ========================================================================

    subroutine calc_dg_sgs_amd_viscosity_node(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg, filt_len
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes

        integer :: state_flag

        real :: blend

        real, allocatable, save:: node_sum(:), node_peaky_sum(:)
        real, allocatable, save :: node_filter_lengths(:,:)
        integer, allocatable, save :: node_visits(:)
        real, allocatable, save :: node_weight_sum(:)
        real, save :: minmaxlen(2)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_reference_density, have_filter_field

        ! Reference density
        real :: rho, mu

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val, sgs_ele_av
        integer :: i, j, ln, gnode

        integer, allocatable :: u_cg_ele(:)

        ! AMD stuff
        real, allocatable :: dx(:)
        real, dimension(:, :), allocatable :: S, dudx_n, del_gradu, B
        real :: BS, topbit, btmbit, Cpoin

        ! Cell-size adaptive viscosity
        real :: thismaxlen, scalelen, startlen, endlen
        real :: startvisc, endvisc, adaptvisc
        real :: valpha


        print*, "In calc_dg_sgs_amd_viscosity()"

        t1=mpi_wtime()

        allocate( S(opDim,opDim),       dudx_n(opDim,opDim), &
                 del_gradu(opDim,opDim), B(opDim,opDim), &
                 dx(opDim) )
        allocate( u_cg_ele(opNloc) )

        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)


        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)

        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for tensor DG LES. (This should not happen)")
        end if

        ! Extract FilterLengths field for debugging
        filt_len=>extract_vector_field(state, "FilterLengths", stat=state_flag)
        if(state_flag /= 0 ) then
            have_filter_field=.false.
        else
            have_filter_field=.true.
        end if


        ! Molecular viscosity
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

        ! The Poincare constant (default 0.3)
        if(have_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/amd/" &
            // "poincare_constant")) then

            call get_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/amd/" &
                // "poincare_constant", Cpoin)
        else
            Cpoin=0.3
        end if

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_sum)) then
                deallocate(node_sum)
                deallocate(node_peaky_sum)
                deallocate(node_visits)
                deallocate(node_weight_sum)
                deallocate(node_filter_lengths)
            end if

            allocate(node_sum(num_nodes))
            allocate(node_peaky_sum(num_nodes))
            allocate(node_weight_sum(num_nodes))
            allocate(node_visits(num_nodes))
            
            allocate(node_filter_lengths(opDim, num_nodes))
            call aniso_filter_lengths(x, node_filter_lengths, minmaxlen)

        end if

        node_sum(:)=0.
        node_peaky_sum(:)=0.
        node_visits(:)=0.
        node_weight_sum(:)=0.


        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

        do e=1, num_elements
            u_cg_ele=ele_nodes(u_cg, e)

            ! call amd_length_ele(x, e, dx)

            sgs_ele_av=0.0
            do ln=1, opNloc
                gnode = u_cg_ele(ln)
                dx(:)=node_filter_lengths(:,gnode)

                dudx_n = u_grad%val(:,:,gnode)

                S = 0.5 * (dudx_n + transpose(dudx_n))

                do i=1, opDim
                    do j=1, opDim
                        del_gradu(i,j) = dx(j) * dudx_n(j,i)
                    end do
                end do

                B = matmul(transpose(del_gradu), del_gradu)

                BS=0.0
                do i=1, opDim
                    do j=1, opDim
                        BS = BS+B(i,j)*S(i,j)
                    end do
                end do

                topbit = rho * Cpoin * max(-BS, 0.)


                btmbit=0.
                do i=1, opDim
                    do j=1, opDim
                        btmbit = btmbit + dudx_n(i,j)**2
                    end do
                end do

                ! If the dominator is vanishing small, then set the SGS viscosity
                ! to zero

                if(btmbit < 10e-10) then
                    sgs_visc_val = 0.
                else
                    sgs_visc_val = topbit/btmbit
                    if(sgs_visc_val > mu*10e6) sgs_visc_val=0.
                end if

                ! Contributions of local node to shared node value.
                node_peaky_sum(gnode) = node_peaky_sum(gnode) + sgs_visc_val

                sgs_ele_av = sgs_ele_av + sgs_visc_val
            end do

            sgs_ele_av = sgs_ele_av / opNloc

            ! Give each corner node the average
            do ln=1, opNloc
                gnode = u_cg_ele(ln)
                node_sum(gnode) = node_sum(gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
            end do

        end do

        ! Set final values.

        ! Blend of element-averaged value and local node averaged value
        ! blend=0...1 (0=peaky, 1=smoothed)
        blend=0.5
        do n=1, num_nodes
            if(have_filter_field) then
                call set(filt_len, n, node_filter_lengths(:,n))
            end if
           

            sgs_visc_val = (blend*node_sum(n) + (1-blend)*node_peaky_sum(n)) &
                         / node_visits(n)

            ! Hard limit again.
            ! sgs_visc_val = min(sgs_visc_val, mu*10e5)
            if(sgs_visc_val > mu*10e5) sgs_visc_val=0.

            ! Lastly, to adaptive mesh scaling for large cell sizes
            thismaxlen = maxval(node_filter_lengths(:,n))
            scalelen = minmaxlen(2)-minmaxlen(1)
            startvisc = 0.0
            endvisc = 100.0*mu

            startlen = 5*minmaxlen(1)
            endlen = 0.5*minmaxlen(2)

            if(thismaxlen < startlen) then
                adaptvisc=0.0
            elseif(thismaxlen >= startlen .and.  thismaxlen < endlen) then
                valpha = (thismaxlen-startlen)/scalelen

                adaptvisc=(1-valpha)*startvisc + valpha*endvisc
            else
                adaptvisc=endvisc
            end if

            call set(sgs_visc, n,  sgs_visc_val+adaptvisc)
        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        call halo_update(filt_len)

        call deallocate(u_grad)
        deallocate( S, dudx_n, del_gradu, B, u_cg_ele, dx )

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_amd_viscosity_node



    ! ========================================================================
    ! Use LES for anisotropic grids, using vorticity-derived filter lengths.
    ! Based upon Chauvet (2007)
    ! ========================================================================

    subroutine calc_dg_sgs_chauvet_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc, dist_to_wall

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg, filt_len
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes

        integer :: state_flag

        real, allocatable, save:: node_sum(:), node_peaky_sum(:)
        real, allocatable, save :: node_filter_lengths(:,:)
        real, allocatable, save :: node_vort_lengths(:)
        integer, allocatable, save :: node_visits(:)
        real, allocatable, save :: node_weight_sum(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_dist_to_wall, have_reference_density, have_filter_field

        ! Reference density
        real :: rho, mu

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val, sgs_ele_av, Cs, vortlen
        integer :: i, j, ln, gnode

        real :: blend

        integer, allocatable :: u_cg_ele(:)

        real, dimension(:, :), allocatable :: S, dudx_n

        print*, "In calc_dg_sgs_chauvet_viscosity()"

        t1=mpi_wtime()


        allocate( S(opDim,opDim), dudx_n(opDim,opDim) )
        allocate( u_cg_ele(opNloc) )

        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)


        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)

        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for tensor DG LES. (This should not happen)")
        end if


        ! Extract FilterLengths field for debugging
        filt_len=>extract_vector_field(state, "FilterLengths", stat=state_flag)
        if(state_flag /= 0 ) then
            have_filter_field=.false.
        else
            have_filter_field=.true.
        end if


        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
        ! (This will be checked for elsewhere)

        mviscosity => extract_tensor_field(state, "Viscosity", stat=state_flag)

        if(mviscosity%field_type == FIELD_TYPE_CONSTANT &
            .or. mviscosity%field_type == FIELD_TYPE_NORMAL) then
            mu=mviscosity%val(1,1,1)
        else
            FLAbort("DG_LES: must have constant dynamic viscosity field")
        end if

        dist_to_wall=>extract_scalar_field(state, "DistanceToWall", stat=state_flag)

        have_dist_to_wall=.false.
        if (state_flag==0) then
            print*, "Using DistanceToWall for Chauvet LES."
            have_dist_to_wall=.true.
        else
            dist_to_wall=>extract_scalar_field(state, "DistanceToBottom", stat=state_flag)
            print*, "Using DistanceToBottom as substitute for DistanceToWall for Chauvet LES."
            have_dist_to_wall=.true.
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
            // "smagorinsky_coefficient", Cs)

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_sum)) then
                deallocate(node_sum)
                deallocate(node_peaky_sum)
                deallocate(node_visits)
                deallocate(node_weight_sum)
                deallocate(node_filter_lengths)
                deallocate(node_vort_lengths)
            end if

            allocate(node_sum(num_nodes))
            allocate(node_peaky_sum(num_nodes))
            allocate(node_weight_sum(num_nodes))
            allocate(node_visits(num_nodes))

            allocate(node_filter_lengths(opDim, num_nodes))
            allocate(node_vort_lengths(num_nodes))
            call aniso_filter_lengths(x, node_filter_lengths)

        end if

        node_sum(:)=0.
        node_peaky_sum(:)=0.
        node_visits(:)=0.
        node_weight_sum(:)=0.

        call vorticity_filter_lengths(x, u_grad, &
                    node_filter_lengths, node_vort_lengths)

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

#ifdef OLD
        do e=1, num_elements
            u_cg_ele=ele_nodes(u_cg, e)

            sgs_ele_av=0.0
            do ln=1, opNloc
                gnode = u_cg_ele(ln)

                dudx_n = u_grad%val(:,:,gnode)
                S = 0.5 * (dudx_n + transpose(dudx_n))

                sgs_visc_val = ((Cs*node_vort_lengths(gnode))**2.) &
                            * rho *norm2(2.*S)

                ! Contributions of local node to shared node value.
                node_peaky_sum(gnode) = node_peaky_sum(gnode) + sgs_visc_val

                sgs_ele_av = sgs_ele_av + sgs_visc_val
            end do

            sgs_ele_av = sgs_ele_av / opNloc

            ! Give each corner node the average
            do ln=1, opNloc
                gnode = u_cg_ele(ln)
                node_sum(gnode) = node_sum(gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
            end do

        end do

        ! Set final values.

        ! Blend of element-averaged value and local node averaged value
        ! blend=0...1 (0=peaky, 1=smoothed)
        blend=0.5
        do n=1, num_nodes
            if(have_filter_field) then
                call set(filt_len, n, node_filter_lengths(:,n))
            end if


            sgs_visc_val = (blend*node_sum(n) + (1-blend)*node_peaky_sum(n)) &
                         / node_visits(n)

            ! Hard limit again.
            ! sgs_visc_val = min(sgs_visc_val, mu*10e5)
            if(sgs_visc_val > mu*10e5) sgs_visc_val=0.

            call set(sgs_visc, n,  sgs_visc_val)
        end do
#else
        do n=1, num_nodes
            dudx_n = u_grad%val(:,:,n)
            S = 0.5 * (dudx_n + transpose(dudx_n))

            vortlen = node_vort_lengths(n)

            sgs_visc_val = ((Cs*vortlen)**2.) &
                        * rho *norm2(2.*S)

            if(sgs_visc_val > mu*100) then
               print *, "sgs_visc_val:", sgs_visc_val
               print *, "vortlen:", vortlen
               print *, "|S|:", norm2(2.*s)

            end if

            call set(sgs_visc, n,  sgs_visc_val)

        end do
#endif
        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        !call halo_update(filt_len)

        call deallocate(u_grad)
        deallocate( S, dudx_n, u_cg_ele )

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_chauvet_viscosity



    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

    subroutine calc_dg_sgs_amd_viscosity(state, x, u)
        type(state_type), intent(in) :: state
        type(vector_field), intent(in) :: u, x

        character(len=OPTION_PATH_LEN) :: scalar_eddy_visc_path
        character(len=256) :: mesh_name

        ! It's not  bullet-proof, but consistent with the logic of DG_prep.F90
        ! Please replace with something better...

        if(have_option("/geometry/mesh::ZeroMesh")) then
            call calc_dg_sgs_amd_viscosity_ele(state, x, u)
        else
            call calc_dg_sgs_amd_viscosity_node(state, x, u)
        end if


    end subroutine calc_dg_sgs_amd_viscosity
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


    subroutine calc_dg_sgs_amd_viscosity_ele(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements
        integer :: state_flag

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_reference_density

        ! Reference density
        real :: rho, mu

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val, sgs_ele_av
        integer :: i, j, ln, gnode

        integer, allocatable :: u_cg_ele(:)

        ! AMD stuff
        real, allocatable, save :: del(:,:)
        real, dimension(:, :), allocatable :: S, dudx_n, del_gradu, B
        real :: BS, topbit, btmbit, Cpoin


        print*, "In calc_dg_sgs_amd_viscosity_ele()"

        t1=mpi_wtime()

        allocate( S(opDim,opDim),       dudx_n(opDim,opDim), &
                 del_gradu(opDim,opDim), B(opDim,opDim) )
        allocate( u_cg_ele(opNloc) )

        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)
        call halo_update(u_grad)

        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for tensor DG LES. (This should not happen)")
        end if

        ! Molecular viscosity
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

        ! The Poincare constant (default 0.3)
        if(have_option(trim(u%option_path)//"/prognostic/" &
            // "spatial_discretisation/discontinuous_galerkin/les_model/amd/" &
            // "poincare_constant")) then

            call get_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/amd/" &
                // "poincare_constant", Cpoin)
        else
            Cpoin=0.3
        end if

        num_elements = ele_count(u_cg)

        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(del)) then
                deallocate(del)
            end if

            allocate(del(u%dim, num_elements))

            call aniso_filter_lengths(X, del)
        end if

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

        do e=1, num_elements
            u_cg_ele=ele_nodes(u_cg, e)

            ! Go around
            sgs_ele_av=0.0
            do ln=1, opNloc
                gnode = u_cg_ele(ln)

                dudx_n = u_grad%val(:,:,gnode)

                S = 0.5 * (dudx_n + transpose(dudx_n))

                do i=1, opDim
                    do j=1, opDim
                        del_gradu(i,j) = del(j,gnode) * dudx_n(j,i)
                    end do
                end do

                B = matmul(transpose(del_gradu), del_gradu)

                BS=0.0
                do i=1, opDim
                    do j=1, opDim
                        BS = BS+B(i,j)*S(i,j)
                    end do
                end do

                topbit = rho * Cpoin * max(-BS, 0.)


                btmbit=0.
                do i=1, opDim
                    do j=1, opDim
                        btmbit = btmbit + dudx_n(i,j)**2
                    end do
                end do

                ! If the dominator is vanishing small, then set the SGS viscosity
                ! to zero

                if(btmbit < 10e-10) then
                    sgs_visc_val = 0.
                else
                    sgs_visc_val = min(topbit/btmbit, mu*10e5)
                end if

                sgs_ele_av = sgs_ele_av + sgs_visc_val
            end do

            sgs_ele_av = sgs_ele_av / opNloc

            ! If we're over this value, reset
            if(sgs_visc_val > mu*10e4) sgs_visc_val=0.

            call set(sgs_visc, e,  sgs_visc_val)
        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)
        deallocate( S, dudx_n, del_gradu, B, u_cg_ele )

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_amd_viscosity_ele



    ! ========================================================================
    ! Horizontal and vertical lengthscales for anisotropic meshes,
    ! from Roman (2010). Only works for 3D linear tets.
    ! ========================================================================

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
         vort(1) = gr(3,2)-gr(2,3)
         vort(2) = gr(1,3)-gr(3,1)
         vort(3) = gr(1,2)-gr(2,1)

         magvort = sqrt(vort(1)**2 + vort(2)**2 + vort(3)**2)

         ! Chauvet et al.
         N = vort / magvort
         vort_del(i) = sqrt( dx(2, i) * dx(3,i) * N(1)**2 &
              + dx(3, i) * dx(1,i) * N(2)**2 &
              + dx(1, i) * dx(2,i) * N(3)**2 )

      end do

      deallocate(gr, vort, N)

    end subroutine vorticity_filter_lengths



    ! ========================================================================
    ! Aniso lengthscales, but really only for pancake elements (small dz).
    ! ========================================================================


    subroutine aniso_length_ele(pos, ele, del, extruded_mesh)
        type(vector_field), intent(in) :: pos
        integer :: ele
        real, intent(inout) :: del(:)
        real :: dx(pos%dim)
        
        real :: X_val(pos%dim, opNloc)
        real :: X_mean(pos%dim), r, diffx, diffy
        real :: X_tri(pos%dim-1, 3)
        real :: area, a, b, c, s, tmplen
        integer :: i, n, m, trix
        integer :: stpair(2)

        ! Not used for now. But will later need conditional to switch in/out
        ! specialised extruded mesh case logic
        logical :: extruded_mesh

        X_val=ele_val(pos, ele)

        ! First look for two points vertically aligned. These two will provide
        ! dz metric. (There are always two in an extruded mesh)

        stpair=0
        outer_loop: do n=1, opNloc
            do m=1, opNloc
                if ( n /= m ) then
                    diffx = abs(X_val(1, n)-X_val(1,m))
                    diffy = abs(X_val(2, n)-X_val(2,m))

                    if(diffx < 10e-10 .and. diffy < 10e-10) then
                        stpair(1)=n
                        stpair(2)=m
                        exit outer_loop
                    end if
                end if
            end do
        end do outer_loop
        del(3) = abs( X_val(3,stpair(1))-X_val(3,stpair(2)) )

        ! Find 2D points for horizontal triangle (easier/quicker than using
        ! Fluidity framework)
        trix=1
        do n=1, opNloc
            if(n /= stpair(1)) then
                X_tri(:, trix)=X_val(1:2, n)
                trix=trix+1
            end if
        end do

        ! Heron's formula for area of triangle
        a = sqrt( (X_tri(1,2)-X_tri(1,1))**2. + (X_tri(2,2)-X_tri(2,1))**2.)
        b = sqrt( (X_tri(1,3)-X_tri(1,2))**2. + (X_tri(2,3)-X_tri(2,2))**2.)
        c = sqrt( (X_tri(1,3)-X_tri(1,1))**2. + (X_tri(2,3)-X_tri(2,1))**2.)

        s = 0.5*(a+b+c)

        area = 0.25 * sqrt( s*(s-a)*(s-b)*(s-c) )
        
        if(area>10e5) print*, "WARNING: AMD LES metrics: large area"

        ! Calculate radius of circle with same area
        r = sqrt(area/ 3.141592653)

        del(1) = 2.*r
        del(2) = del(1)

        ! Not quite done. Sanity check for very, very thin elements
        if(del(3)/del(1) < 0.05) del(3)=0.05 * del(1)

    end subroutine aniso_length_ele



    ! This one is not per-element as above, but loops over all elements
    subroutine aniso_filter_lengths(pos, del_fin, lenrange)
        type(vector_field), intent(in) :: pos
        real, dimension(:,:), intent(inout) :: del_fin

        real, optional :: lenrange(2)
        real :: tmplen

        real, dimension(pos%dim, opNloc) :: X_val
        real, dimension(pos%dim) :: dx, del

        real, allocatable :: dx_sum(:,:), dx_ele_raw(:,:), dx_ele_filt(:,:)
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
        allocate(dx_ele_filt(pos%dim, num_elements))

        ! Reset counters and sums
        visits=0
        dx_sum=0.

        ! Is this an extruded mesh?
        extruded_mesh = option_count("/geometry/mesh/from_mesh/extrude") > 0

        if( .not. extruded_mesh ) then
            FLExit("Error: AMD LES is currently only works correctly on extruded meshes")
        end if

        ! Raw sizes per element
        do e=1, num_elements
           call aniso_length_ele(pos, e, dx, extruded_mesh=extruded_mesh)

!           dx_ele_raw(:,e) = dx(:)
           dx_ele_filt(:,e) = dx(:)
        end do

        ! Filter sizes per element
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
            local_gnodes = ele_nodes(pos, e)

            X_val=ele_val(pos, e)

            do n=1, opNloc
               gn=local_gnodes(n)

               if(gn>0 .and. gn<=num_nodes) then

                   ! Add element sizes to node size sum at each corner
                   dx_sum(:, gn) = dx_sum(:, gn) + dx_ele_filt(:,e)
                   visits(gn) = visits(gn)+1
               else
                   print*, "gn:", gn
               end if
            end do
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
        deallocate(dx_ele_filt)
        deallocate(visits)

    end subroutine aniso_filter_lengths



    ! This one is per-element
    subroutine amd_filter_lengths_ele_old(pos, del_fin)
        type(vector_field), intent(in) :: pos
        real, dimension(:,:), intent(inout) :: del_fin

        real, allocatable :: del_ele_raw(:, :)

        integer, pointer :: neighs(:)

        real, dimension(pos%dim, opNloc) :: X_val
        real :: mean_val, dx_neigh_sum(pos%dim)

        integer :: e, n, i, gn, f
        integer :: num_elements, num_neighs

        num_elements = ele_count(pos)

        allocate(del_ele_raw(pos%dim, num_elements))

        ! Raw sizes per element
        do e=1, num_elements
           X_val=ele_val(pos, e)
           do i=1, opDim
              mean_val = sum(X_val(i,:)) / opNloc
              del_ele_raw(i, e) = 2.*sum( abs(X_val(i, :)-mean_val) ) / opNloc
           end do
        end do

        ! Filtered (smoothed) element sizes
        do e=1, num_elements
           X_val=ele_val(pos, e)
           neighs=>ele_neigh(pos, e)
           num_neighs = size(neighs)

           dx_neigh_sum=0.
           do f=1, num_neighs
              dx_neigh_sum = dx_neigh_sum + del_ele_raw(:,neighs(f))
           end do

           del_fin(:,e) = (del_ele_raw(:,e) + dx_neigh_sum)/(num_neighs+1)
        end do

        deallocate(del_ele_raw)


    end subroutine amd_filter_lengths_ele_old



end module dg_les

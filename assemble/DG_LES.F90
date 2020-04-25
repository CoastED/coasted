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

#define FILTER_TYPE 1



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

  ! This scales the Van Driest effect
  real, parameter :: van_scale=1.0

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

#ifdef SHARP_LES_FILTER
                node_sum(gnode) = visc_turb
                node_visits(gnode) = 1
#else
                sgs_ele_av = sgs_ele_av + visc_turb/opNloc
#endif
            end do

#ifndef SHARP_LES_FILTER
            do ln=1, opNloc
                gnode = u_cg_ele(ln)

                node_sum(gnode) = node_sum(gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
            end do
#endif

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






    ! ========================================================================
    ! Split horizontal/vertical LES SGS viscosity
    ! ========================================================================

    subroutine calc_dg_sgs_tensor_viscosity(state, x, u)
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
        real, dimension(u%dim, u%dim) :: u_grad_node, rate_of_strain
        real :: mag_strain_horz, mag_strain_vert, mag_strain_r
        real :: sgs_horz, sgs_vert, sgs_r
        real, dimension(u%dim, u%dim) :: sgs_ele_av, visc_turb, tmp_tensor
        real :: mu, rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable,save:: node_vol_weighted_sum(:,:,:),node_neigh_total_vol(:)
        real, allocatable,save:: node_sum(:, :,:), sgs_unfiltered(:,:,:)
        integer, allocatable, save :: node_visits(:)

        real :: visc_norm2, visc_norm2_max

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0

        real, parameter :: alpha=0.5

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val
        integer :: i

        print*, "In calc_dg_sgs_tensor_viscosity()"

#if FILTER_TYPE == 1
        print*, "- using original LES filter"
#elif FILTER_TYPE == 2
        print*, "- using box LES filter"
#elif FILTER_TYPE == 3
        print*, "- using hat LES filter"
#else
#error "Unsupported Large Eddy Simulation FILTER_TYPE"
#endif

        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)

        sgs_visc_mag => extract_scalar_field(state, "TensorEddyViscosityMagnitude", stat=state_flag)

        if (state_flag /= 0) then
            FLAbort("DG_LES: TensorEddyViscosityMagnitude absent for tensor DG LES. (This should not happen)")
         end if

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")

!        ! Viscosity. Here we assume isotropic viscosity, ie. Newtonian fluid
!        ! (This will be checked for elsewhere)

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
                deallocate(sgs_unfiltered)
            end if

            allocate(node_sum(u%dim, u%dim, num_nodes))
            allocate(node_visits(num_nodes))
            allocate(node_vol_weighted_sum(u%dim, u%dim, num_nodes))
            allocate(node_neigh_total_vol(num_nodes))
            allocate(sgs_unfiltered(u%dim, u%dim, num_nodes))
        end if

        node_sum(:,:,:)=0.0
        node_visits(:)=0
        node_vol_weighted_sum(:,:,:)=0.0
        node_neigh_total_vol(:)=0.0
        
        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:,:,:)=0.0


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


#if FILTER_TYPE == 1 || FILTER_TYPE == 2
                ! Original and box filter needs the element average
                sgs_ele_av = sgs_ele_av + visc_turb
#else
                ! hat-type filter averages over values co-occupying nodes
                node_sum(:,:, gnode) = node_sum(:,:, gnode) + visc_turb
                node_visits(gnode) = node_visits(gnode) + 1
                node_vol_weighted_sum(:,:, gnode) = node_vol_weighted_sum(:,:, gnode) + ele_vol*visc_turb
                node_neigh_total_vol(gnode) = node_neigh_total_vol(gnode) + ele_vol
#endif

            end do

            sgs_ele_av = sgs_ele_av / opNloc
#if FILTER_TYPE == 1
            ! Extra calculations for original filter
            do ln=1, opNloc
                gnode = u_cg_ele(ln)

                node_sum(:,:, gnode) = node_sum(:,:, gnode) + sgs_ele_av
                node_visits(gnode) = node_visits(gnode) + 1
                node_vol_weighted_sum(:,:, gnode) = node_vol_weighted_sum(:,:, gnode) + ele_vol*sgs_ele_av
                node_neigh_total_vol(gnode) = node_neigh_total_vol(gnode) + ele_vol

            end do
#endif
        end do

        ! For box filter, we'll just do all the stuff in here, rather than loop
        ! round the global nodes
#if FILTER_TYPE == 2
        if(have_van_driest) then
            do e=1, num_elements
                u_cg_ele=ele_nodes(u_cg, e)

                do ln=1, opNloc
                    gnode = u_cg_ele(ln)
                    u_grad_node = u_grad%val(:,:, gnode)
                    y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(gnode)
                    vd_damping = (1-exp(-y_plus/A_plus))**pow_m

                    tmp_tensor = vd_damping * rho * alpha * sgs_ele_av

                    ! L2 Norm of tensor (tensor magnitude)
                    sgs_visc_val=0.0
                    do i=1, opDim
                        sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
                    end do
                    sgs_visc_val = sqrt(sgs_visc_val)

                    call set(sgs_visc, gnode,  tmp_tensor)
                    call set(sgs_visc_mag, gnode, sgs_visc_val )
                end do
            end do
        else
            do e=1, num_elements
                u_cg_ele=ele_nodes(u_cg, e)

                do ln=1, opNloc
                    gnode = u_cg_ele(ln)
                    tmp_tensor = rho * alpha * sgs_ele_av

                    ! L2 Norm of tensor (tensor magnitude)
                    sgs_visc_val=0.0
                    do i=1, opDim
                        sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
                    end do
                    sgs_visc_val = sqrt(sgs_visc_val)

                    call set(sgs_visc, gnode,  tmp_tensor)
                    call set(sgs_visc_mag, gnode, sgs_visc_val )
                end do
            end do
        end if
#endif

        ! This part only works for original and hat filters
#if FILTER_TYPE == 1 || FILTER_TYPE == 3
        ! Set final values. Two options here: one with Van Driest damping, one without.
        ! We multiply by rho here.
        if(have_van_driest) then

            do n=1, num_nodes
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m

#if FILTER_TYPE == 1
                tmp_tensor = vd_damping * rho * &
                ( alpha * sgs_unfiltered(:,:,n) + (1-alpha) * node_sum(:,:,n)/node_visits(n) )

#elif FILTER_TYPE == 3
                tmp_tensor = vd_damping * rho * node_sum(:,:,n) &
                     / node_visits(n)
#endif

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

#if FILTER_TYPE == 1
                tmp_tensor = rho * &
                ( alpha * sgs_unfiltered(:,:,n) + (1-alpha) * node_sum(:,:,n)/node_visits(n) )
#elif FILTER_TYPE == 3
                tmp_tensor = rho * node_sum(:,:,n) / node_visits(n)
#endif
!                tmp_tensor = rho * node_vol_weighted_sum(:,:,n) &
!                     / node_neigh_total_vol(n))

                ! L2 Norm of tensor (tensor magnitude)
                sgs_visc_val=0.0
                do i=1, opDim
                    sgs_visc_val = sgs_visc_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
                end do
                sgs_visc_val = sqrt(sgs_visc_val)

                call set(sgs_visc, n, tmp_tensor )
                call set(sgs_visc_mag, n, sgs_visc_val )

            end do
        end if
        ! End of code for original and hat filters
#endif

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        call halo_update(sgs_visc_mag)

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_tensor_viscosity


    ! ========================================================================
    ! Use Anisotropic Minimum Dissipation (AMD) LES filter.
    ! See Rozema et al, Computational Methods in Engineering, 2020.
    ! ========================================================================

    subroutine calc_dg_sgs_amd_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: dist_to_wall
        type(scalar_field), pointer :: sgs_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes, ln

        real :: ele_vol
        integer :: state_flag, gnode

        integer :: not_first_call
        real :: blend

        real, allocatable,save:: node_vol_weighted_sum(:,:,:), node_neigh_total_vol(:)
        real, allocatable,save:: node_sum(:, :,:)
        integer, allocatable, save :: node_visits(:)
        real, allocatable, save :: node_weight_sum(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_wall_distance, have_reference_density

        ! Reference density
        real :: rho

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val
        integer :: i, j, d

        ! AMD stuff
        real, allocatable, save :: del(:,:)
        real, dimension(u%dim, u%dim) :: S, dudx_n, del_gradu, B
        real :: BS, topbit, btmbit, Cpoin


        print*, "In calc_dg_sgs_amd_viscosity()"

        t1=mpi_wtime()

        nullify(dist_to_wall)
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

        ! Do we have a DistanceToWall field?
        have_wall_distance=.false.
        dist_to_wall=>extract_scalar_field(state, "DistanceToWall", stat=state_flag)
        if (state_flag==0) then
            print*, "DistanceToWall field: applying lengthscale limiting"
            have_wall_distance=.true.
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
                deallocate(node_visits)
                deallocate(node_vol_weighted_sum)
                deallocate(node_weight_sum)
                deallocate(node_neigh_total_vol)
                deallocate(del)
            end if

            allocate(node_sum(u%dim, u%dim, num_nodes))
            allocate(node_weight_sum(num_nodes))
            allocate(node_visits(num_nodes))
            allocate(node_vol_weighted_sum(u%dim, u%dim, num_nodes))
            allocate(node_neigh_total_vol(num_nodes))
            allocate(del(u%dim, num_nodes))

            ! We can apply filter-length limiting if we have a distance_to_wall field
            if(have_wall_distance) then
                call amd_filter_lengths(X, del, dist_to_wall)
            else
                call amd_filter_lengths(X, del)
            end if
        end if

        node_sum(:,:,:)=0.0
        node_visits(:)=0
        node_vol_weighted_sum(:,:,:)=0.0
        node_neigh_total_vol(:)=0.0

        ! Set entire SGS visc field to zero value initially
        if( not_first_call==0 ) then
            sgs_visc%val(:)=0.0
        end if
        not_first_call=1

        ! Set final values. We multiply by rho here.
        do n=1, num_nodes
            dudx_n = u_grad%val(:,:,n)

            S = 0.5 * (dudx_n + transpose(dudx_n))

            do concurrent (i=1:opDim, j=1:opDim)
                del_gradu(i,j) = del(j,n) * dudx_n(j,i)
            end do

            B = matmul(transpose(del_gradu), del_gradu)

            BS=0.0

            do i=1, opDim
                do j=1, opDim
                    BS = BS+B(i,j)*S(i,j)
                end do
            end do

            topbit = rho * Cpoin * max(-BS, 0.)

!            AS = 0.
!            do i=1, opDim
!                do j=1, opDim
!                    AS = AS + &
!                    (del(1,n)*dudx_n(i,1)+del(2,n)*dudx_n(i,2)+del(3,n)*dudx_n(i,3))&
!                   *(del(1,n)*dudx_n(j,1)+del(2,n)*dudx_n(j,2)+del(3,n)*dudx_n(j,3))&
!                   *S(i,j)
!
!! Not sure if this one isn't correct.
!!                    (del(1,n)*dudx_n(1,i)+del(2,n)*dudx_n(2,i)+del(3,n)*dudx_n(3,i))&
!!                   *(del(1,n)*dudx_n(1,j)+del(2,n)*dudx_n(2,j)+del(3,n)*dudx_n(3,j))&
!!                   *S(i,j)
!                end do
!            end do
!
!
!            topbit=max(-AS,0.)

            btmbit=0.
            do i=1, opDim
                do j=1, opDim
                    btmbit = btmbit + dudx_n(i,j)**2
                end do
            end do

            ! If the dominator is vanishing small, then set the SGS viscosity
            ! to zero
            ! We use a randomly perturbated blend between old and new values.
            !blend=0.5*rand()+0.5
            blend=0.75
            if(btmbit < 10e-10) then
                sgs_visc_val = 0. + (1-blend)* sgs_visc%val(n)
            else
                sgs_visc_val = blend * (topbit / btmbit) &
                    + (1-blend)*sgs_visc%val(n)

            end if

            call set(sgs_visc, n,  sgs_visc_val)
        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_amd_viscosity



    ! ========================================================================
    ! Horizontal and vertical lengthscales for anisotropic meshes,
    ! from Roman (2010). Only works for 3D linear tets.
    ! ========================================================================

    subroutine les_length_scales_squared_mk2(positions, ele, horzSq, vertSq)
        type(vector_field), intent(in) :: positions
        integer :: ele

        real :: horzSq, vertSq
        real, dimension(opDim, opNloc) :: X_val
        real, dimension(opDim) :: dx

        integer :: i

        X_val=ele_val(positions, ele)

        ! Calculate largest dx, dy, dz for element nodes
        do i=1, opDim
            dx(i) = maxval( X_val(i,:)) - minval (X_val(i,:))
        end do

        ! Why squares? Used in LES visc calcs, no need for expensive square roots
        horzSq = dx(1)**2 + dx(2)**2
        vertSq = dx(3)**2

    end subroutine les_length_scales_squared_mk2


    ! This one is not per-element as above, but loops over all elements
    subroutine amd_filter_lengths(pos, del_fin, distwall)
        type(vector_field), intent(in) :: pos
        real, dimension(:,:), intent(inout) :: del_fin
        type(scalar_field), optional, intent(in) :: distwall

        real, dimension(pos%dim, opNloc) :: X_val
        real, dimension(pos%dim) :: dx, del

        real, allocatable :: dx_sum(:,:), dx_ele_raw(:,:), dx_ele_filt(:,:)
        real, dimension(pos%dim) :: dx_neigh_sum, dx_neigh_average
        real :: del_max, mean_val
        integer, allocatable :: visits(:)
        integer, pointer :: neighs(:)
        integer :: local_gnodes(ele_loc(pos,1))

        integer :: e, n, i, gn, f
        integer :: num_elements, num_nodes, num_neighs

        logical :: have_wall_distance

        num_elements = ele_count(pos)
        num_nodes = node_count(pos)

        allocate(dx_sum(pos%dim, num_nodes))
        allocate(visits(num_nodes))
        allocate(dx_ele_raw(pos%dim, num_elements))
        allocate(dx_ele_filt(pos%dim, num_elements))

        ! Check for distance to wall field as argument
        have_wall_distance=present(distwall)

        ! Reset counters and sums
        visits=0
        dx_sum=0.


        ! Raw sizes per element
        do e=1, num_elements
           X_val=ele_val(pos, e)
           do i=1, opDim
              mean_val = sum(X_val(i,:)) / opNloc
              dx(i) = 2.*sum( abs(X_val(i, :)-mean_val) ) / opNloc
           end do

           dx_ele_raw(:,e) = dx(:)
        end do

        ! Filtered (smoothed) element sizes
        do e=1, num_elements
           X_val=ele_val(pos, e)
           neighs=>ele_neigh(pos, e)
           num_neighs = size(neighs)

           dx_neigh_sum=0.
           do f=1, num_neighs
              dx_neigh_sum = dx_neigh_sum + dx_ele_raw(:,neighs(f))
           end do

           dx_ele_filt(:,e) = (dx_ele_raw(:,e) + dx_neigh_sum)/(num_neighs+1)
        end do

        ! Now create sizes per-node
        do e=1, num_elements
            local_gnodes = ele_nodes(pos, e)

            X_val=ele_val(pos, e)

            do n=1, opNloc
               gn=local_gnodes(n)

               ! Add element sizes to node size sum at each corner
               dx_sum(:, gn) = dx_sum(:, gn) + dx_ele_filt(:,e)
               visits(gn) = visits(gn)+1
            end do
        end do

        ! Now go around all nodes.
        do n=1, num_nodes
            del = 2.0*dx_sum(:, n)/visits(n)
            del_max=maxval(del)

            ! Do some sensible limiting
            if(have_wall_distance) then
                do i=1, opDim
                    if(abs(del(i)/del_max) < 0.05) del(i)=0.05*del_max
                    del_fin(i, n) = (min(del(i), abs(distwall%val(n))))
                end do
            else
                do i=1, opDim
                    if(abs(del(i)/del_max) < 0.05) del(i)=0.05*del_max
                    del_fin(i, n) = del(i)
                end do
            end if
        end do

        deallocate(dx_sum)
        deallocate(dx_ele_raw)
        deallocate(dx_ele_filt)
        deallocate(visits)

    end subroutine amd_filter_lengths



end module dg_les

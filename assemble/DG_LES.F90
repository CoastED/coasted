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

#define NDIM 3
#define NLOC 4

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

  implicit none

  private

  public :: calc_dg_sgs_scalar_viscosity, calc_dg_sgs_tensor_viscosity

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
        real :: sgs_ele_av, visc_turb, mu
        real :: rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable, save :: node_sum(:), node_vol_weighted_sum(:), &
             node_neigh_total_vol(:)
        integer, allocatable, save :: node_visits(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=25.6, pow_m=2.0

        print*, "calc_dg_sgs_scalar_viscosity"

        t1=mpi_wtime()

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)


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
                deallocate(node_vol_weighted_sum)
                deallocate(node_visits)
                deallocate(node_neigh_total_vol)
            end if

            allocate(node_sum(num_nodes))
            allocate(node_vol_weighted_sum(num_nodes))
            allocate(node_visits(num_nodes))
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

            length=2.*(ele_vol**0.333333333333333333333)
            ! Factor of two included in length_scale_scalar
            Cs_length_sq = (Cs* length)**2.0

            ! This is the contribution to nu_sgs from each co-occupying node
            sgs_ele_av=0.0
            do ln=1, NLOC
                gnode = u_cg_ele(ln)
                rate_of_strain = 0.5 * (u_grad%val(:,:, gnode) + transpose(u_grad%val(:,:, gnode)))
                visc_turb = Cs_length_sq * norm2(2.0 * rate_of_strain)

                sgs_ele_av = sgs_ele_av + visc_turb/NLOC
                node_sum(gnode) = node_sum(gnode) + visc_turb
                node_visits(gnode) = node_visits(gnode) + 1
            end do

            ! This is the weighted contribution from each element
            do ln=1, NLOC
                gnode = u_cg_ele(ln)

                node_vol_weighted_sum(gnode) = node_vol_weighted_sum(gnode) + ele_vol*sgs_ele_av
                node_neigh_total_vol(gnode) = node_neigh_total_vol(gnode) + ele_vol
            end do
        end do

        ! Set final values. Two options here: one with Van Driest damping, 
        ! one without.
        if(have_van_driest) then
            do n=1, num_nodes
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node+transpose(u_grad_node)) * rho) * dist_to_wall%val(n)
                vd_damping = 1.0 - exp(-y_plus/A_plus)

                call set(sgs_visc, n, &
                    vd_damping * rho*0.5*(node_sum(n) / node_visits(n) &
                    + node_vol_weighted_sum(n) / node_neigh_total_vol(n)) )

            end do
        else
            do n=1, num_nodes
                call set(sgs_visc, n, &
                    rho*0.5*(node_sum(n) / node_visits(n) &
                    + node_vol_weighted_sum(n) / node_neigh_total_vol(n)) )
            end do
        end if

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

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes, ln
        integer :: u_cg_ele(ele_loc(u,1))

        real :: Cs_horz, Cs_length_horz_sq, Cs_vert, Cs_length_vert_sq
        real :: length_horz_sq, length_vert_sq, ele_vol
        real, dimension(u%dim, u%dim) :: u_grad_node, rate_of_strain
        real :: mag_strain_horz, mag_strain_vert
        real :: sgs_horz, sgs_vert
        real, dimension(u%dim, u%dim) :: sgs_ele_av, visc_turb
        real :: mu, rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real, allocatable,save:: node_sum(:, :,:), node_vol_weighted_sum(:,:,:),&
             node_neigh_total_vol(:)
        integer, allocatable, save :: node_visits(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=25.0, pow_m=2.0

        print*, "In calc_dg_sgs_tensor_viscosity()"

        t1=mpi_wtime()

        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)


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
                deallocate(node_vol_weighted_sum)
                deallocate(node_visits)
                deallocate(node_neigh_total_vol)
            end if

            allocate(node_sum(u%dim, u%dim, num_nodes))
            allocate(node_vol_weighted_sum(u%dim, u%dim, num_nodes))
            allocate(node_visits(num_nodes))
            allocate(node_neigh_total_vol(num_nodes))
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
            do ln=1, NLOC
                gnode = u_cg_ele(ln)
                rate_of_strain = 0.5 * (u_grad%val(:,:, gnode) + transpose(u_grad%val(:,:, gnode)))

                mag_strain_horz = sqrt(2.0* rate_of_strain(1,1)**2.0 &
                                                + 2.0* rate_of_strain(2,2)**2.0 &
                                                + 4.0* rate_of_strain(1,2)**2.0 )

                mag_strain_vert = sqrt(4.0* rate_of_strain(1,3)**2.0 &
                                                + 2.0* rate_of_strain(3,3)**2.0 &
                                                + 4.0* rate_of_strain(3,2)**2.0 )

                ! Note, this is without density. That comes later.
                sgs_horz = Cs_length_horz_sq * mag_strain_horz
                sgs_vert = Cs_length_vert_sq * mag_strain_vert

                ! As per Roman et al, 2010.
                visc_turb(1:2, 1:2) = sgs_horz
                visc_turb(3, :) = sgs_vert
                visc_turb(:, 3) = sgs_vert

                sgs_ele_av = sgs_ele_av + visc_turb/NLOC
                node_sum(:,:, gnode) = node_sum(:,:, gnode) + visc_turb
                node_visits(gnode) = node_visits(gnode) + 1
            end do

            ! This is the weighted contribution from each element
            do ln=1, NLOC
                gnode = u_cg_ele(ln)

                node_vol_weighted_sum(:,:, gnode) = node_vol_weighted_sum(:,:, gnode) + ele_vol*sgs_ele_av
                node_neigh_total_vol(gnode) = node_neigh_total_vol(gnode) + ele_vol
            end do
        end do

        ! Set final values. Two options here: one with Van Driest damping, one without.
        ! We multiply by rho here.
        if(have_van_driest) then
            do n=1, num_nodes
                u_grad_node = u_grad%val(:,:, n)
                y_plus = sqrt(norm2(u_grad_node) * rho * mu) * dist_to_wall%val(n) 
                vd_damping =(( 1- exp(-y_plus/A_plus))**pow_m)*van_scale+(1-van_scale)

                call set(sgs_visc, n, &
                    vd_damping * rho*0.5*(node_sum(:,:,n) / node_visits(n) &
                    + node_vol_weighted_sum(:,:,n) / node_neigh_total_vol(n)) )
            end do
        else
            do n=1, num_nodes
                call set(sgs_visc, n, &
                    rho*0.5*(node_sum(:,:,n) / node_visits(n) &
                    + node_vol_weighted_sum(:,:, n) / node_neigh_total_vol(n)) )
            end do
        end if

        call deallocate(u_grad)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

    end subroutine calc_dg_sgs_tensor_viscosity




    ! ========================================================================
    ! Horizontal and vertical lengthscales for anisotropic meshes,
    ! from Roman (2010). Only works for 3D linear tets.
    ! ========================================================================

    subroutine les_length_scales_squared_mk2(positions, ele, horzSq, vertSq)
        type(vector_field), intent(in) :: positions
        integer :: ele

        real :: horzSq, vertSq
        real, dimension(NDIM, NLOC) :: X_val
        real, dimension(NDIM) :: dx

        integer :: i

        X_val=ele_val(positions, ele)

        ! Calculate largest dx, dy, dz for element nodes
        do i=1, NDIM
            dx(i) = maxval( X_val(i,:)) - minval (X_val(i,:))
        end do

        ! Why squares? Used in LES visc calcs, no need for expensive square roots
        horzSq = dx(1)**2 + dx(2)**2
        vertSq = dx(3)**2

    end subroutine les_length_scales_squared_mk2
end module dg_les

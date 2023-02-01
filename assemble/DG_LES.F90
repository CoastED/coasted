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

!        call calc_dg_sgs_amd_viscosity_node(state, x, u)
        ! Calling QR for now, see how it does.
!        call calc_dg_sgs_qr_viscosity(state, x, u)
!        call calc_dg_sgs_amd_viscosity_new(state, x, u)
!        call calc_dg_sgs_amd_viscosity_new_mk2(state, x, u)

    end subroutine calc_dg_sgs_amd_viscosity




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


    end subroutine calc_dg_sgs_scalar_viscosity








    ! =========================================================================
    !  Vreman Model.
    ! =========================================================================


    subroutine calc_dg_sgs_vreman_viscosity(state, x, u)
        ! Passed parameters
        type(state_type), intent(in) :: state

        type(vector_field), intent(in) :: u, x
        type(scalar_field), pointer :: sgs_visc
        type(scalar_field), pointer :: artificial_visc

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad
        integer :: n, num_nodes, e, num_elements
        integer :: state_flag

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_reference_density, have_filter_field
        logical :: have_artificial_visc

        ! Reference density
        real :: rho, mu

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val, sgs_ele_av
        integer :: i, j, k

        ! Vreman stuff
        real, dimension(:, :), allocatable :: alpha_ij, beta_ij, dudx
        real :: alpha_sum, B_beta, Cpoin
        integer :: udim
        type(vector_field), pointer :: elelen, nodelen

        print*, "In calc_dg_sgs_vreman_viscosity()"

        t1=mpi_wtime()

        allocate( alpha_ij(opDim,opDim), beta_ij(opDim,opDim), dudx(opDim,opDim))

! I think this is the suspect call -- not needed.
        nullify(mviscosity)

        ! Length scales for filter
        elelen => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodelen => extract_vector_field(state, "NodeLengthScales", stat=state_flag)

        if(state_flag == 0) then
           print*, "**** Have ScalarElementLengthScales and ScalarNodeLengthScales fields"
        else
           FLAbort("Error: must have ElementLengthsScales and NodeLengthScales fields for Scalar Vreman DG LES. This should have been automatically created.")
        end if

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        udim = u_cg%dim

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_scalar_field(state, "ScalarEddyViscosity", &
             stat=state_flag)
        if (state_flag /= 0) then
            FLAbort("DG_LES: ScalarEddyViscosity absent for DG Vreman LES. (This should not happen)")
        end if


        ! We can use this in areas of insufficient resolution
        have_artificial_visc = .false.
        artificial_visc => extract_scalar_field(state, "ArtificialViscosity", &
             stat=state_flag)
        if(state_flag == 0) then
           print*, "ArtificialViscosity field detected."
           have_artificial_visc = .true.
        end if

        call grad(u_cg, x, u_grad)
        ! Crucially, update halos for use
        call halo_update(u_grad)

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

        num_nodes = u_cg%mesh%nodes


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

        do n=1, num_nodes

!           dudx(:,1) = u_grad%val(:,n)
!           dudx(:,2) = v_grad%val(:,n)
!           dudx(:,3) = w_grad%val(:,n)
           dudx(:,1) = u_grad%val(1,:,n)
           dudx(:,2) = u_grad%val(2,:,n)
           dudx(:,3) = u_grad%val(3,:,n)

           alpha_ij = dudx

           alpha_sum = 0.0
           beta_ij = 0.0           

           do i=1, opDim
              do j=1, opDim
                 alpha_sum = alpha_sum + alpha_ij(i,j)*alpha_ij(i,j)
                 do k=1, opDim
                    beta_ij = beta_ij + (nodelen%val(k,n)**2.)*alpha_ij(k,i)*alpha_ij(k,j)
                 end do
              end do
           end do
           
           B_beta = &
                  (beta_ij(1,1)*beta_ij(2,2) - beta_ij(1,2)**2.) &
                + (beta_ij(1,1)*beta_ij(3,3) - beta_ij(1,3)**2.) &
                + (beta_ij(2,2)*beta_ij(3,3) - beta_ij(2,3)**2.) 

           if( alpha_sum < 10e-10) then
              sgs_visc_val = 0.0
           else
              sgs_visc_val = Cpoin * sqrt( B_beta / alpha_sum )
           end if
           
           if(have_artificial_visc) then
              sgs_visc_val = sgs_visc_val + artificial_visc%val(n)
           end if

            ! Limiter
            if(sgs_visc_val > mu*10e4) then
               if(sgs_visc%val(n) < mu*10e4) then
                  sgs_visc_val=sgs_visc%val(n)
               else
                  sgs_visc_val=0.
               end if
            end if

           call set(sgs_visc, n,  sgs_visc_val)
        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)

        call deallocate(u_grad)
        deallocate( alpha_ij, beta_ij, dudx )

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
        type(scalar_field), pointer :: sgs_visc, dist_to_wall

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg, filt_len
        type(tensor_field), pointer :: mviscosity
        type(tensor_field) :: u_grad

        integer :: e, num_elements, n, num_nodes

        integer :: state_flag

        real, allocatable, save :: node_filter_lengths(:,:)
        real, allocatable :: node_vort_lengths(:)

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density, have_filter_field

        ! Reference density
        real :: rho, mu

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val, sgs_ele_av, Cs
        integer :: i, j, ln, gnode

        real :: blend

        integer, allocatable :: u_cg_ele(:)

        real, dimension(:, :), allocatable :: S, dudx_n

        print*, "In calc_dg_sgs_chauvet_viscosity()"

        t1=mpi_wtime()


        allocate( S(opDim,opDim), dudx_n(opDim,opDim) )
        allocate( u_cg_ele(opNloc) )

        ! Van Driest wall damping
        have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")

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
            // "smagorinsky_coefficient", Cs)

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        allocate(node_vort_lengths(num_nodes))
        ! We only allocate if mesh connectivity unchanged from
        ! last iteration; reuse saved arrays otherwise

        if(new_mesh_connectivity) then
            if(allocated(node_filter_lengths)) then
                deallocate(node_filter_lengths)
            end if

            allocate(node_filter_lengths(opDim, num_nodes))
            call aniso_filter_lengths(x, node_filter_lengths)

        end if


        call vorticity_filter_lengths(x, u_grad, &
                    node_filter_lengths, node_vort_lengths)

        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:)=0.0

        do n=1, num_nodes
            dudx_n = u_grad%val(:,:,n)
            S = 0.5 * (dudx_n + transpose(dudx_n))

            sgs_visc_val = ((Cs*node_vort_lengths(n))**2.) &
                        * rho *norm2(2.*S)

            call set(sgs_visc, n,  sgs_visc_val)

        end do

        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        ! call halo_update(filt_len)

        call deallocate(u_grad)
        deallocate( S, dudx_n )
        deallocate( node_vort_lengths)

        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

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
        type(scalar_field), pointer :: dist_to_wall
        type(tensor_field), pointer :: sgs_visc
        type(scalar_field), pointer :: sgs_visc_mag

        ! Velocity (CG) field, pointer to X field, and gradient
        type(vector_field), pointer :: u_cg
        type(tensor_field), pointer :: mviscosity
        type(scalar_field), pointer :: artificial_visc

        type(vector_field), pointer :: elelen, nodelen
        type(tensor_field) :: u_grad
        real, allocatable :: dx_ele_raw(:,:)

        integer :: e, ix, num_elements, n, num_nodes
        integer :: u_cg_ele(ele_loc(u,1))

        real :: minmaxlen(2)

        real :: Cs_horz, Cs_length_horz_sq, Cs_vert, Cs_length_vert_sq
        real :: length_horz_sq, length_vert_sq, ele_vol
        real, dimension(u%dim, u%dim) :: u_grad_node, rate_of_strain
        real :: mag_strain_horz, mag_strain_vert
        real :: sgs_horz, sgs_vert

        real :: sgs_visc_mag_val

        real, dimension(u%dim, u%dim) :: visc_turb, tmp_tensor
        real, dimension(u%dim) :: tmp_vector
        real :: mu, rho, y_plus, vd_damping

        integer :: state_flag, gnode

        real :: visc_norm2, visc_norm2_max

        real (kind=8) :: t1, t2
        real (kind=8), external :: mpi_wtime

        logical :: have_van_driest, have_reference_density
        logical :: have_lengths_field

        ! Constants for Van Driest damping equation
        real, parameter :: A_plus=17.8, pow_m=2.0

        real, parameter :: alpha=0.5

        ! For scalar tensor eddy visc magnitude field
        real :: sgs_visc_val
        integer :: i

        logical :: have_artificial_visc

        print*, "In calc_dg_sgs_roman_viscosity()"

        t1=mpi_wtime()

!        nullify(elelen)
!        nullify(nodelen)
        nullify(dist_to_wall)
        nullify(mviscosity)

        ! Velocity projected to continuous Galerkin
        u_cg=>extract_vector_field(state, "VelocityCG", stat=state_flag)

        ! Allocate gradient field
        call allocate(u_grad, u_cg%mesh, "VelocityCGGradient")

        sgs_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=state_flag)
        call grad(u_cg, x, u_grad)

        ! Length scales for filter
        elelen => extract_vector_field(state, "ElementLengthScales", stat=state_flag)
        nodelen => extract_vector_field(state, "NodeLengthScales", stat=state_flag)

        if(state_flag == 0) then
           print*, "**** Have ElementLengthScales and NodeLengthScales fields"
           have_lengths_field = .true.
        else
           FLAbort("Error: must have ElementLengthsScales and NodeLengthScales fields for Roman DG LES. This should have been automatically created.")
        end if
        sgs_visc_mag => extract_scalar_field(state, "TensorEddyViscosityMagnitude", stat=state_flag)

        ! We can use this in areas of insufficient resolution
        have_artificial_visc = .false.
        artificial_visc => extract_scalar_field(state, "ArtificialViscosity", &
             stat=state_flag)
        if(state_flag == 0) then
           print*, "ArtificialViscosity field detected."
           have_artificial_visc = .true.
        end if

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
            // "spatial_discretisation/discontinuous_galerkin/les_model/roman/" &
            // "smagorinsky_coefficient_vertical")) then

            call get_option(trim(u%option_path)//"/prognostic/" &
                // "spatial_discretisation/discontinuous_galerkin/les_model/roman/" &
                // "smagorinsky_coefficient_vertical", Cs_vert)
        else
            FLAbort("DG_LES: you've requested anisotropic LES, but have not specified smagorinsky_coefficient_vertical")
        end if

        num_elements = ele_count(u_cg)
        num_nodes = u_cg%mesh%nodes

        allocate(dx_ele_raw(3, num_elements))

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


        ! Set entire SGS visc field to zero value initially
        sgs_visc%val(:,:,:)=0.0

        ! calculate it at each node
        do n=1, num_nodes

            Cs_length_horz_sq = (Cs_horz**2) * (nodelen%val(1, n)**2 + nodelen%val(2, n)**2)
            Cs_length_vert_sq = (Cs_vert * nodelen%val(3, n))**2

            ! This is the contribution to nu_sgs from each co-occupying node
            rate_of_strain = 0.5 * (u_grad%val(:,:, n) + transpose(u_grad%val(:,:, n)))

            mag_strain_horz = sqrt(2.0* rate_of_strain(1,1)**2.0 &
                                 + 2.0* rate_of_strain(2,2)**2.0 &
                                 + 4.0* rate_of_strain(1,2)**2.0 )

            mag_strain_vert = sqrt(4.0* rate_of_strain(1,3)**2.0 &
                             + 2.0* rate_of_strain(3,3)**2.0 &
                             + 4.0* rate_of_strain(3,2)**2.0 )

            ! Note, this is without density. That comes later.
            sgs_horz = rho * Cs_length_horz_sq * mag_strain_horz
            sgs_vert = rho * Cs_length_vert_sq * mag_strain_vert

            ! As per Roman et al, 2010.
            visc_turb(1:2, 1:2) = sgs_horz
            visc_turb(1, 3) = sgs_vert
            visc_turb(3, 1) = sgs_vert
            visc_turb(2, 3) = sgs_vert
            visc_turb(3, 2) = sgs_vert
            visc_turb(3, 3) = sgs_vert
            
            ! visc_turb(3, 3) = sgs_horz + (-2.*sgs_vert) + 2.*sgs_r)
            ! Otherwise this is potentially negative (!)
            ! visc_turb(3, 3) = sgs_horz + 2.*sgs_vert + 2.*sgs_r

            ! According to Roman et al, 2010, only two sgs viscosities are used.
            ! visc_turb(3, 3) = sgs_vert


            ! Account for wall-damping if enabled
            if(have_van_driest) then
                y_plus = sqrt(norm2(u_grad%val(:,:, n)) * rho / mu) * dist_to_wall%val(n)
                vd_damping = (1-exp(-y_plus/A_plus))**pow_m
            else
            ! No Van Driest, no damping
                vd_damping = 1.0
            end if

            tmp_tensor = vd_damping * rho * visc_turb

            if(have_artificial_visc) then
                tmp_tensor(:,:) = tmp_tensor(:,:) + artificial_visc%val(n)
            end if

            ! L2 Norm of tensor (tensor magnitude)
            sgs_visc_mag_val=0.0
            do i=1, opDim
                sgs_visc_mag_val = sgs_visc_mag_val + dot_product(tmp_tensor(i,:),tmp_tensor(i,:))
            end do
            sgs_visc_mag_val = sqrt(sgs_visc_mag_val)

            call set(sgs_visc, n,  tmp_tensor)
            call set(sgs_visc_mag, n, sgs_visc_mag_val )

        end do



        ! Must be done to avoid discontinuities at halos
        call halo_update(sgs_visc)
        call halo_update(sgs_visc_mag)

        call deallocate(u_grad)
        deallocate(dx_ele_raw)


        t2=mpi_wtime()

        print*, "**** DG_LES_execution_time:", (t2-t1)

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
!                            stpair(1)=n
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
!                if(n /= stpair(1)) then
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

!            ! For unstructured meshes, do something slightly different.
    ! Do this for both structured and unstructured meshes, for now.
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
           del(1) = ( (6*vol) / (pi*del(3)) )**0.5
           del(2) = del(1)
        end if

!        end if

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

        ! call quick_smooth_ele_lengths_field(state, 3)

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

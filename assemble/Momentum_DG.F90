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

#include "fdebug.h"
#include "compile_opt_defs.h"


module momentum_DG
    ! This module contains the Discontinuous Galerkin form of the momentum
    ! equation.
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
    use sparsity_patterns
    use petsc_tools
    use turbine
    use diagnostic_fields
    use slope_limiters_dg
    use smoothing_module
    use fields_manipulation
    use field_derivatives
    use field_options
    use sparsity_patterns_meshes
    use colouring
    use Profiler
#ifdef _OPENMP
  use omp_lib
#endif
    use multiphase_module

    use dg_les

    implicit none

    ! Buffer for output messages.
    character(len=255), private :: message
    private
    public construct_momentum_dg, &
        momentum_DG_check_options, correct_velocity_dg, &
        assemble_poisson_rhs_dg, allocate_big_m_dg, &
        subcycle_momentum_dg

    ! Module private variables for model options. This prevents us having to
    ! do dictionary lookups for every element (or element face!)
    real :: dt, theta, theta_nl
    logical :: lump_mass, lump_abs, lump_source, subcycle

    ! Whether the advection term is only integrated by parts once.
    logical :: integrate_by_parts_once=.false.
    ! Whether the conservation term is integrated by parts or not
    logical :: integrate_conservation_term_by_parts=.false.
    ! Whether or not to integrate the surface tension term by parts
    logical :: integrate_surfacetension_by_parts

    ! Weight between conservative and non-conservative forms of the advection
    ! equation.
    ! 1 is for conservative 0 is for non-conservative.
    real :: beta

    ! Discretisation to use for viscosity term.
    integer :: viscosity_scheme
    integer, parameter :: ARBITRARY_UPWIND=1
    integer, parameter :: BASSI_REBAY=2
    integer, parameter :: CDG=3
    integer, parameter :: IP=4

    ! Method for getting h0 in IP
    integer :: edge_length_option
    integer, parameter :: USE_FACE_INTEGRALS=1
    integer, parameter :: USE_ELEMENT_CENTRES=2

    ! Parameters for interior penalty method
    real :: Interior_Penalty_Parameter, edge_length_power, h0

    ! Flag indicating whether equations are being solved in acceleration form.
    logical :: acceleration

    ! Flag indicating whether to include pressure bcs (not for cv pressure)
    logical :: l_include_pressure_bcs
  
    ! which terms do we have?
    logical :: have_mass
    logical :: have_source
    logical :: have_gravity
    logical :: on_sphere
    logical :: have_absorption
    logical :: have_vertical_stabilization
    logical :: have_implicit_buoyancy
    logical :: have_vertical_velocity_relaxation
    logical :: have_swe_bottom_drag
    ! implicit absorption is corrected by the pressure correction
    ! by combining the implicit part of absorption with the mass term of u^{n+1}
    logical :: pressure_corrected_absorption
    logical :: have_viscosity
    logical :: have_surfacetension
    logical :: have_coriolis
    logical :: have_advection
    logical :: move_mesh
    logical :: have_pressure_bc
    logical :: subtract_out_reference_profile
  
    real :: gravity_magnitude

    ! CDG stuff
    real, dimension(3) :: switch_g
    logical :: CDG_penalty
    logical :: remove_penalty_fluxes

    ! Are we running a multi-phase flow simulation?
    logical :: multiphase

    ! Declaring assembly arrays at module level
    real, allocatable, dimension(:, :) :: Coriolis_mat, rho_mat, rho_move_mat, mass_mat, inverse_mass_mat, Advection_mat, Source_mat
    real, allocatable, dimension(:, :, :) :: ele2grad_mat, Abs_mat
    real, allocatable, dimension(:, :, :, :) :: Abs_mat_sphere
    real, allocatable, dimension(:, :) :: Abs_lump
    real, allocatable, dimension(:, :, :) :: Abs_lump_sphere
    real, allocatable, dimension(:) :: source_lump
    real, allocatable, dimension(:, :) :: Q_inv
    real, allocatable, dimension(:, :, :) :: Grad_u_mat_q, Div_u_mat_q
    real, allocatable, dimension(:, :, :, :) :: Viscosity_mat
    real, allocatable, dimension(:) :: node_stress_diag, resid_stress_term

    real, allocatable, dimension(:, :) :: big_m_diag_addto, rhs_addto
    real, allocatable, dimension(:, :, :, :) :: big_m_tensor_addto
    logical, allocatable, dimension(:, :) :: diagonal_block_mask, off_diagonal_block_mask
    real, allocatable, dimension(:, :, :, :) :: subcycle_m_tensor_addto
    real, allocatable, dimension(:, :, :, :) :: sh_tensout, sh_outtens
    real, allocatable, dimension(:, :) :: sh_to_temp
    real, allocatable, dimension(:, :) :: sh_ot_temp
    real, allocatable, dimension(:, :) :: sh_dt_temp

    real, allocatable, dimension(:, :, :) :: Viscosity_ele
    real, allocatable, dimension(:, :, :) :: visc_ele_quad
    real, allocatable, dimension(:, :) :: x_val, x_val_2, u_val
    real, allocatable, dimension(:, :, :, :) :: kappa_mat
    real, allocatable, dimension(:) :: l_MassLump, l_move_masslump
    integer, allocatable, dimension(:) :: local_glno
    real, allocatable, dimension(:) :: detwei, detwei_old, detwei_new, coefficient_detwei, detwei_rhoq
    real, allocatable, dimension(:, :, :) :: du_t, dug_t, dq_t
    real, allocatable, dimension(:) :: Rho_q, Coriolis_q
    real, allocatable, dimension(:, :) :: u_nl_q
    real, allocatable, dimension(:) :: u_nl_div_q
    real, allocatable, dimension(:, :, :) :: tension
    real, allocatable, dimension(:, :) :: dtensiondj

    real, allocatable, dimension(:, :) :: absorption_gi
    real, allocatable, dimension(:, :, :) :: tensor_absorption_gi
    real, allocatable, dimension(:, :, :) :: vvr_abs
    real, allocatable, dimension(:, :) :: vvr_abs_diag
    real, allocatable, dimension(:) :: depth_at_quads
    real, allocatable, dimension(:, :, :) :: ib_abs
    real, allocatable, dimension(:, :) :: ib_abs_diag
    real, allocatable, dimension(:, :, :) :: dt_rho
    real, allocatable, dimension(:, :) :: grav_at_quads, grad_rho
    real, allocatable, dimension(:, :) :: ele_grav_val
    real, allocatable, dimension(:) :: drho_dz

    real, allocatable, dimension(:, :) :: ele_u_mesh_quad
    real, allocatable, dimension(:) :: ele_centre, neigh_centre, face_centre, face_centre_2
    real, allocatable, dimension(:) :: alpha_u_quad

    real, allocatable, dimension(:) :: face_Rho_q
    real, allocatable, dimension(:) :: face_Rho_val
    real, allocatable, dimension(:, :) :: face_normal, face_u_nl_q, face_u_f_q, face_u_f2_q, face_div_u_f_q
    real, allocatable, dimension(:, :) :: face_u_mesh_quad
    logical, allocatable, dimension(:) :: face_inflow
    real, allocatable, dimension(:) :: face_u_nl_q_dotn, face_income
    real, allocatable, dimension(:) :: face_detwei, face_detwei_work
    real, allocatable, dimension(:) :: face_inner_advection_integral, face_outer_advection_integral
    real, allocatable, dimension(:, :) :: face_nnAdvection_out
    real, allocatable, dimension(:, :) :: face_nnAdvection_in
    real, allocatable, dimension(:, :, :, :) :: face_mnCT
    real, allocatable, dimension(:, :, :) :: face_kappa_gi
    real, allocatable, dimension(:, :, :) :: face_visc_val
    real, allocatable, dimension(:, :, :) :: face_tension_q

    real, allocatable, dimension(:, :, :) :: tmp_face_tensor
    real, allocatable, dimension(:, :, :) :: face_primal_fluxes_mat
    real, allocatable, dimension(:, :) :: face_shape_shape_work
    real, allocatable, dimension(:, :, :) :: face_penalty_fluxes_mat
    real, allocatable, dimension(:, :, :) :: face_normal_mat
    real, allocatable, dimension(:, :, :) :: face_kappa_normal_mat

    real, allocatable, dimension(:) :: matmul_dut_visc


    ! =======================================================================
    ! Make assembly arrays private to each OpenMP thread
    ! =======================================================================
    !$OMP THREADPRIVATE(Coriolis_mat, rho_mat, rho_move_mat, mass_mat, inverse_mass_mat, Advection_mat, Source_mat)
    !$OMP THREADPRIVATE(ele2grad_mat, Abs_mat)
    !$OMP THREADPRIVATE(Abs_mat_sphere)
    !$OMP THREADPRIVATE(Abs_lump)
    !$OMP THREADPRIVATE(Abs_lump_sphere)
    !$OMP THREADPRIVATE(source_lump)
    !$OMP THREADPRIVATE(Q_inv)
    !$OMP THREADPRIVATE(Grad_u_mat_q, Div_u_mat_q)
    !$OMP THREADPRIVATE(Viscosity_mat)
    !$OMP THREADPRIVATE(node_stress_diag, resid_stress_term)
    !$OMP THREADPRIVATE(big_m_diag_addto, rhs_addto)

    !$OMP THREADPRIVATE(big_m_tensor_addto)
    !$OMP THREADPRIVATE(diagonal_block_mask, off_diagonal_block_mask)
    !$OMP THREADPRIVATE(subcycle_m_tensor_addto)
    !$OMP THREADPRIVATE(sh_tensout, sh_outtens)
    !$OMP THREADPRIVATE(sh_to_temp)
    !$OMP THREADPRIVATE(sh_ot_temp)
    !$OMP THREADPRIVATE(sh_dt_temp)

    !$OMP THREADPRIVATE(Viscosity_ele)
    !$OMP THREADPRIVATE(visc_ele_quad)
    !$OMP THREADPRIVATE(x_val, x_val_2, u_val)
    !$OMP THREADPRIVATE(kappa_mat)
    !$OMP THREADPRIVATE(l_MassLump, l_move_masslump)
    !$OMP THREADPRIVATE(local_glno)
    !$OMP THREADPRIVATE(detwei, detwei_old, detwei_new, coefficient_detwei, detwei_rhoq)
    !$OMP THREADPRIVATE(du_t, dug_t, dq_t)
    !$OMP THREADPRIVATE(Rho_q, Coriolis_q)
    !$OMP THREADPRIVATE(u_nl_q)
    !$OMP THREADPRIVATE(u_nl_div_q)
    !$OMP THREADPRIVATE(tension)
    !$OMP THREADPRIVATE(dtensiondj)

    !$OMP THREADPRIVATE(absorption_gi)
    !$OMP THREADPRIVATE(tensor_absorption_gi)
    !$OMP THREADPRIVATE(vvr_abs)
    !$OMP THREADPRIVATE(vvr_abs_diag)
    !$OMP THREADPRIVATE(depth_at_quads)
    !$OMP THREADPRIVATE(ib_abs)
    !$OMP THREADPRIVATE(ib_abs_diag)
    !$OMP THREADPRIVATE(dt_rho)
    !$OMP THREADPRIVATE(grav_at_quads, grad_rho)
    !$OMP THREADPRIVATE(ele_grav_val)
    !$OMP THREADPRIVATE(drho_dz)

    !$OMP THREADPRIVATE(ele_u_mesh_quad)
    !$OMP THREADPRIVATE(ele_centre, neigh_centre, face_centre, face_centre_2)
    !$OMP THREADPRIVATE(alpha_u_quad)

    !$OMP THREADPRIVATE(face_Rho_q)
    !$OMP THREADPRIVATE(face_Rho_val)
    !$OMP THREADPRIVATE(face_normal, face_u_nl_q, face_u_f_q, face_u_f2_q, face_div_u_f_q)
    !$OMP THREADPRIVATE(face_u_mesh_quad)
    !$OMP THREADPRIVATE(face_inflow)
    !$OMP THREADPRIVATE(face_u_nl_q_dotn, face_income)
    !$OMP THREADPRIVATE(face_detwei, face_detwei_work)
    !$OMP THREADPRIVATE(face_inner_advection_integral, face_outer_advection_integral)
    !$OMP THREADPRIVATE(face_nnAdvection_out)
    !$OMP THREADPRIVATE(face_nnAdvection_in)
    !$OMP THREADPRIVATE(face_mnCT)
    !$OMP THREADPRIVATE(face_kappa_gi)
    !$OMP THREADPRIVATE(face_visc_val)
    !$OMP THREADPRIVATE(face_tension_q)

    !$OMP THREADPRIVATE(tmp_face_tensor)
    !$OMP THREADPRIVATE(face_primal_fluxes_mat)
    !$OMP THREADPRIVATE(face_shape_shape_work)
    !$OMP THREADPRIVATE(face_penalty_fluxes_mat)
    !$OMP THREADPRIVATE(face_normal_mat)
    !$OMP THREADPRIVATE(face_kappa_normal_mat)

    !$OMP THREADPRIVATE(matmul_dut_visc)

    
contains


    ! Conditionally including the optimised CDG assembly code.

#include "Construct_Momentum_Element_DG_opt.F90"


    subroutine construct_momentum_dg(u, p, rho, x, &
        & big_m, rhs, state, &
        & inverse_masslump, inverse_mass, mass, &
        & acceleration_form, include_pressure_bcs,&
        & subcycle_m)
        !!< Construct the momentum equation for discontinuous elements in
        !!< acceleration form. If acceleration_form is present and false, the
        !!< matrices will be constructed for use in conventional solution form.

        !! velocity and coordinate
        type(vector_field), intent(inout) :: u, x
        !! pressure and density
        type(scalar_field), intent(inout) :: p, rho

        !! Main momentum matrix.
        type(petsc_csr_matrix), intent(inout) :: big_m
        !! Explicit subcycling matrix.
        type(block_csr_matrix), intent(inout), optional :: subcycle_m
        !! Momentum right hand side vector for each point.
        type(vector_field), intent(inout) :: rhs
        !! Collection of fields defining system state.
        type(state_type) :: state
    
        !! Inverse of the lumped mass lumping at each point.
        !! NOTE: only allocated and calculated if (lump_mass)
        type(vector_field), intent(inout), optional :: inverse_masslump
        !! Optional separate mass matrix.
        !! NOTE: if provided the mass matrix, won't be added to big_m
        !! NOTE2: this mass matrix does not include density, bcs or absorption factors
        !! NOTE3: mass is not allocated here (unlike inverse_masslump and inverse_mass)
        type(csr_matrix), intent(inout), optional :: mass
        !! Inverse mass matrix
        !! NOTE: only allocated and calculated if (.not. lump_mass)
        !! NOTE2: diagonal blocks may be different due to dirichlet bcs and/or absorption
        type(block_csr_matrix), intent(inout), optional :: inverse_mass

        !! whether to include the dirichlet pressure bc integrals to the rhs
        logical, intent(in), optional :: include_pressure_bcs

        !! Optional indicator of whether we are solving in acceleration form.
        !!
        !! If not solving in acceleration form then big_m will be formed with
        !! an effective theta of 1.0 and dt of 1.0 . In addition, only boundary
        !! terms will be inserted on the right hand side. This is sufficient to
        !! enable the full discrete equations to be written using matrix
        !! multiplies outside this routine.
        logical, intent(in), optional :: acceleration_form

        !! Position, velocity and source fields.
        type(vector_field), pointer :: U_mesh,  X_old, X_new
        type(vector_field), target, save :: U_nl
        !! Projected (non-linear) velocity field
        type(vector_field), target, save ::  pvelocity
        type(vector_field), pointer :: advecting_velocity
        !! Mesh for projected velocity.
        type(mesh_type) :: pmesh
        character(len=FIELD_NAME_LEN) :: pmesh_name

        !! Viscosity
        type(tensor_field) :: Viscosity

        !! Momentum source and absorption fields
        type(scalar_field) :: buoyancy
        type(vector_field) :: Source, gravity, Abs, Abs_wd
        !! Surface tension field
        type(tensor_field) :: surfacetension

        ! Dummy fields in case state doesn't contain the above fields
        type(scalar_field), pointer :: dummyscalar

        ! Fields for the subtract_out_reference_profile option under the Velocity field
        type(scalar_field), pointer :: hb_density, hb_pressure

        !! field over the entire surface mesh, giving bc values
        type(vector_field), save :: velocity_bc
        type(scalar_field) :: pressure_bc
        !! for each surface element, the bc type to be applied there
        !! integer value determined by ordering in call to get_entire_boundary_condition
        integer, dimension(:,:), allocatable :: velocity_bc_type
        integer, dimension(:), allocatable :: pressure_bc_type
    
        !! Sparsity for inverse mass
        type(csr_sparsity):: mass_sparsity
    
        !! Element index
        integer :: ele

        !! Status variable for field extraction.
        integer :: stat

        !! Mesh for auxiliary variable
        type(mesh_type), save :: q_mesh, turbine_conn_mesh

        ! Fields for vertical velocity relaxation
        type(scalar_field), pointer :: dtt, dtb
        type(scalar_field) :: depth
        integer :: node
        real :: vvr_sf ! A scale factor for the absorption
     
        ! Min vertical density gradient for implicit buoyancy
        real :: ib_min_grad
   
        !! Wetting and drying
        type(scalar_field), pointer :: wettingdrying_alpha
        type(scalar_field) :: alpha_u_field
        logical :: have_wd_abs
        real, dimension(u%dim) :: abs_wd_const

        !! shallow water bottom drag
        type(scalar_field) :: swe_bottom_drag, old_pressure
        type(vector_field) :: swe_u_nl

        !!
        type(integer_set), dimension(:), pointer :: colours
        integer :: len, clr, nnid
        !! Is the transform_to_physical cache we prepopulated valid
        logical :: cache_valid
        integer :: num_threads

        ! Volume fraction fields for multi-phase flow simulation
        type(scalar_field), pointer :: vfrac
        type(scalar_field) :: nvfrac ! Non-linear approximation to the PhaseVolumeFraction

        ! Free surface stabilisation
        logical :: have_free_stab
        real :: free_stab_param

        ! Partial stress - sp911
        logical :: partial_stress

        ! LES - sp911
        logical :: have_les = .false., have_isotropic_les=.false.
        logical :: have_amd_les=.false.
        logical :: have_vreman_les=.false.
        logical :: have_roman_les=.false.
        logical :: have_van_driest = .false., have_vel_cg=.false.
        real :: smagorinsky_coefficient
        type(scalar_field), pointer :: eddy_visc, prescribed_filter_width, distance_to_wall, &
            & y_plus_debug, les_filter_width_debug
        type(vector_field), pointer, save :: u_cg
        type(tensor_field), pointer :: tensor_eddy_visc

        type(tensor_field), pointer :: uGrad

        logical :: freeStab

        ! Benchmarking stuff
        real (kind=8) :: t0, t1, assemble_dt, total_dt, percent_dg
        real (kind=8) :: inner_t0, inner_t1, inner_assemble_dt, inner_percent_dg
        real (kind=8), save :: lastt
        real (kind=8), external :: mpi_wtime

        ! Shape functions
        type(element_type), pointer :: u_shape, p_shape, q_shape

        ! element lists per colour
        type ColouredListType
            integer, allocatable :: list(:)
        end type ColouredListType

        type(ColouredListType), allocatable :: coloured_ele_lists(:)

        t0 = mpi_wtime()

        ewrite(1, *) "In construct_momentum_dg"

        call profiler_tic("construct_momentum_dg")
        assert(continuity(u)<0)

        acceleration= .not. present_and_false(acceleration_form)
        ewrite(2, *) "Acceleration form? ", acceleration

        ! Free surface stabilisation
        have_free_stab=.false.
        if (have_option(trim(U%option_path)//"/type::free_surface/surface_stabilisation")) then
            have_free_stab=.true.
            call get_option(trim(U%option_path)//"/type::free_surface/surface_stabilisation/scale_factor", free_stab_param)

        end if

        if(present(include_pressure_bcs)) then
            l_include_pressure_bcs = include_pressure_bcs
        else
            l_include_pressure_bcs = .true.
        end if
    

        ! These names are based on the CGNS SIDS.
        if (.not.have_option(trim(U%option_path)//"/prognostic"//&
            &"/spatial_discretisation/discontinuous_galerkin"//&
            &"/advection_scheme/none")) then
            U_nl=extract_vector_field(state, "NonlinearVelocity")
            call incref(U_nl)
 
            if(.not.has_vector_field(state, "ProjectedNonlinearVelocity")) then
               
                call get_option(trim(U%option_path)//"/prognostic"//&
                    &"/spatial_discretisation/discontinuous_galerkin"//&
                    &"/advection_scheme/project_velocity_to_continuous"//&
                    &"/mesh/name",pmesh_name)
                pmesh = extract_mesh(state, pmesh_name)
                call allocate(pvelocity, U_nl%dim, pmesh, &
                    &"ProjectedNonlinearVelocity")
                call project_field(U, pvelocity, X)
                call insert(state, pvelocity, "ProjectedNonlinearVelocity")
                advecting_velocity => pvelocity
               
                ! Discard the additional reference.
                call deallocate(pvelocity)
            else
                pvelocity = extract_vector_field(state, &
                    &"ProjectedNonlinearVelocity")
            end if

            advecting_velocity => U_nl
            have_advection = .true.

        else
            ! Forcing a zero NonlinearVelocity will disable advection.
            call allocate(U_nl, U%dim,  U%mesh, "NonlinearVelocity", &
                FIELD_TYPE_CONSTANT)
            call zero(U_nl)
            have_advection=.false.
            advecting_velocity => U_nl
        end if
        ewrite(2, *) "Include advection? ", have_advection

        allocate(dummyscalar)
        call allocate(dummyscalar, u%mesh, "DummyScalar", field_type=FIELD_TYPE_CONSTANT)
        call zero(dummyscalar)
        dummyscalar%option_path=""

        Source=extract_vector_field(state, "VelocitySource", stat)
        have_source = (stat==0)
        if (.not.have_source) then
            call allocate(Source, U%dim,  U%mesh, "VelocitySource", FIELD_TYPE_CONSTANT)
            call zero(Source)
        else
            ! Grab an extra reference to cause the deallocate below to be safe.
            call incref(Source)
            ewrite_minmax(source)
        end if

        Abs=extract_vector_field(state, "VelocityAbsorption", stat)
        have_absorption = (stat==0)
        if (.not.have_absorption) then
            call allocate(Abs, U%dim, U%mesh, "VelocityAbsorption", FIELD_TYPE_CONSTANT)
            call zero(Abs)
        else
            ! Grab an extra reference to cause the deallocate below to be safe.
            call incref(Abs)
            ewrite_minmax(Abs)
        end if

        have_wd_abs=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption")
        ! Absorption term in dry zones for wetting and drying
        if (have_wd_abs) then
            call allocate(Abs_wd, U%dim, U%mesh, "VelocityAbsorption_WettingDrying", FIELD_TYPE_CONSTANT)
            call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption", abs_wd_const)
            call set(Abs_wd, abs_wd_const)
        ! else
        !    call zero(Abs_wd)
        end if

        ! Check if we have either implicit absorption term
        have_vertical_stabilization=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation").or. &
            have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")

        ! If we have vertical velocity relaxation set then grab the required fields
        ! sigma = n_z*g*dt*_rho_o/depth
        have_vertical_velocity_relaxation=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation")
        if (have_vertical_velocity_relaxation) then
            call get_option(trim(U%option_path)//"/prognostic/vertical_stabilization/vertical_velocity_relaxation/scale_factor", vvr_sf)
            dtt => extract_scalar_field(state, "DistanceToTop")
            dtb => extract_scalar_field(state, "DistanceToBottom")
            call allocate(depth, dtt%mesh, "Depth")
            do node=1,node_count(dtt)
                call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
            end do
        endif

        ! Implicit buoyancy (theta*g*dt*drho/dr)
        have_implicit_buoyancy=have_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy")
        call get_option(trim(U%option_path)//"/prognostic/vertical_stabilization/implicit_buoyancy/min_gradient"&
            , ib_min_grad, default=0.0)

        have_swe_bottom_drag = have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater/bottom_drag')
        if (have_swe_bottom_drag) then
            ! Note that we don't do this incref business, instead we just pass uninitialised fields if .not. have_swe_bottom_drag
            swe_bottom_drag = extract_scalar_field(state, "BottomDragCoefficient")
            assert(.not. have_vertical_stabilization)
            depth = extract_scalar_field(state, "BottomDepth") ! we reuse the field that's already passed for VVR
            old_pressure = extract_scalar_field(state, "OldPressure")
            call get_option(trim(U%option_path)//&
                &"/prognostic/temporal_discretisation/relaxation", theta_nl)
            ! because of the kludge above with advecting velocity, let's just have our own u_nl
            ! can be on whatever mesh
            swe_u_nl = extract_vector_field(state, "NonlinearVelocity")
        end if

        call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, stat)
        have_gravity = stat==0
        if (have_option(trim(u%option_path)//'/prognostic/equation::ShallowWater')) then
            ! for the swe there's no buoyancy term
            have_gravity = .false.
            assert(stat==0) ! we should have a gravity_magnitude though
        end if

        if(have_gravity) then
            buoyancy=extract_scalar_field(state, "VelocityBuoyancyDensity")
            call incref(buoyancy)
            gravity=extract_vector_field(state, "GravityDirection", stat)
            call incref(gravity)

        else
            call allocate(buoyancy, u%mesh, "VelocityBuoyancyDensity", FIELD_TYPE_CONSTANT)
            call zero(buoyancy)
            call allocate(gravity, u%dim, u%mesh, "GravityDirection", FIELD_TYPE_CONSTANT)
            call zero(gravity)
        end if
        ewrite_minmax(buoyancy)

        ! Splits up the Density and Pressure fields into a hydrostatic component (') and a perturbed component ('').
        ! The hydrostatic components, denoted p' and rho', should satisfy the balance: grad(p') = rho'*g
        ! We subtract the hydrostatic component from the density used in the buoyancy term of the momentum equation.
        if (have_option(trim(state%option_path)//'/equation_of_state/compressible/subtract_out_reference_profile')) then
            subtract_out_reference_profile = .true.
            hb_density => extract_scalar_field(state, "HydrostaticReferenceDensity")

            if(l_include_pressure_bcs) then
                hb_pressure => extract_scalar_field(state, "HydrostaticReferencePressure")
            else
                hb_pressure => dummyscalar
            end if
        else
            subtract_out_reference_profile = .false.
            hb_density => dummyscalar
            hb_pressure => dummyscalar
        end if

        Viscosity=extract_tensor_field(state, "Viscosity", stat)
        have_viscosity = (stat==0)
        if (.not.have_viscosity) then
            call allocate(Viscosity, U%mesh, "Viscosity", FIELD_TYPE_CONSTANT)
            call zero(Viscosity)
        else
            ! Grab an extra reference to cause the deallocate below to be safe.
            call incref(Viscosity)
            ewrite_minmax(viscosity)
        end if

        surfacetension = extract_tensor_field(state, "VelocitySurfaceTension", stat)
        have_surfacetension = (stat == 0)
        if(.not. have_surfacetension) then
            call allocate(surfacetension, u%mesh, "VelocitySurfaceTension", FIELD_TYPE_CONSTANT)
            call zero(surfacetension)
        else
            call incref(surfacetension)
            ewrite_minmax(surfacetension)
        end if

        ! Are we running a multi-phase simulation?
        if(option_count("/material_phase/vector_field::Velocity/prognostic") > 1) then
            multiphase = .true.

            vfrac => extract_scalar_field(state, "PhaseVolumeFraction")
            call allocate(nvfrac, vfrac%mesh, "NonlinearPhaseVolumeFraction")
            call zero(nvfrac)
            call get_nonlinear_volume_fraction(state, nvfrac)

            ewrite_minmax(nvfrac)

        else
            multiphase = .false.
            nullify(vfrac)
        end if

        have_coriolis = have_option("/physical_parameters/coriolis")
    
        q_mesh=Viscosity%mesh

        on_sphere = have_option('/geometry/spherical_earth/')

        ! Extract model parameters from options dictionary.
        if (acceleration) then
            call get_option(trim(U%option_path)//&
                &"/prognostic/temporal_discretisation/theta", theta)
            call get_option("/timestepping/timestep", dt)
        else
            theta=1.0
            dt=1.0
        end if

        have_mass = .not. have_option(trim(u%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/mass_terms/exclude_mass_terms")
        lump_mass=have_option(trim(U%option_path)//&
            &"/prognostic/spatial_discretisation"//&
            &"/discontinuous_galerkin/mass_terms/lump_mass_matrix")
        lump_abs=have_option(trim(U%option_path)//&
            &"/prognostic/vector_field::Absorption"//&
            &"/lump_absorption")
        pressure_corrected_absorption=have_option(trim(u%option_path)//&
            &"/prognostic/vector_field::Absorption"//&
            &"/include_pressure_correction") .or. (have_vertical_stabilization)
        
        if (pressure_corrected_absorption) then
            ! as we add the absorption into the mass matrix
            ! lump_abs needs to match lump_mass
            lump_abs = lump_mass
        end if
        lump_source=have_option(trim(u%option_path)//&
            &"/prognostic/vector_field::Source"//&
            &"/lump_source")
        call get_option(trim(U%option_path)//"/prognostic/spatial_discretisation"//&
            &"/conservative_advection", beta)

        ! mesh movement here only matters for the mass terms
        ! other terms are evaluated using "Coordinate" which is evaluated at t+theta*dt
        move_mesh = have_option("/mesh_adaptivity/mesh_movement") .and. &
            have_mass
        if (move_mesh) then
            X_old => extract_vector_field(state, "OldCoordinate")
            X_new => extract_vector_field(state, "IteratedCoordinate")
            U_mesh => extract_vector_field(state, "GridVelocity")
        end if
    
        ! by default we assume we're integrating by parts twice
        integrate_by_parts_once = have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/once"        )

        integrate_conservation_term_by_parts = have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_conservation_term_by_parts"        )

        ! Determine the scheme to use to discretise viscosity.
        if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/bassi_rebay"            )) then
                viscosity_scheme=BASSI_REBAY
            else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme&
         &/compact_discontinuous_galerkin"                )) then
                    !=================Compact Discontinuous Galerkin
                    viscosity_scheme=CDG
                    !Set the switch vector
                    switch_g = 0.
                    switch_g(1) = exp(sin(3.0+exp(1.0)))
                    if(mesh_dim(U)>1) switch_g(2) = (cos(exp(3.0)/sin(2.0)))**2
                    if(mesh_dim(U)>2) switch_g(3) = sin(cos(sin(cos(3.0))))
                    switch_g = switch_g/sqrt(sum(switch_g**2))

                    remove_penalty_fluxes = .true.
                    interior_penalty_parameter = 0.0
                    if(have_option(trim(U%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/viscosity_scheme"//&
                        &"/compact_discontinuous_galerkin/penalty_parameter")) then
                        remove_penalty_fluxes = .false.
                        edge_length_power = 0.0
                        call get_option(trim(U%option_path)//&
                            &"/prognostic/spatial_discretisation"//&
                            &"/discontinuous_galerkin/viscosity_scheme"//&
                            &"/compact_discontinuous_galerkin/penalty_parameter"&
                            &,Interior_Penalty_Parameter)
                    end if

                    CDG_penalty = .true.
                    edge_length_option = USE_FACE_INTEGRALS

                else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation"//&
                    &"/discontinuous_galerkin/viscosity_scheme/arbitrary_upwind")) then
                    viscosity_scheme=ARBITRARY_UPWIND
                else if (have_option(trim(U%option_path)//&
                    &"/prognostic/spatial_discretisation"//&
                    &"/discontinuous_galerkin/viscosity_scheme/interior_penalty")) then
                    remove_penalty_fluxes = .false.
                    viscosity_scheme=IP
                    CDG_penalty = .false.
                    call get_option(trim(U%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/viscosity_scheme"//&
                        &"/interior_penalty/penalty_parameter",Interior_Penalty_Parameter)
                    call get_option(trim(U%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/viscosity_scheme"//&
                        &"/interior_penalty/edge_length_power",edge_length_power)
                    edge_length_option = USE_FACE_INTEGRALS
                    if(have_option(trim(U%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/viscosity_scheme"//&
                        &"/interior_penalty/edge_length_option/use_element_centres")) then
                        edge_length_option = USE_ELEMENT_CENTRES
                    end if
                else
                    FLAbort("Unknown viscosity scheme - Options tree corrupted?")
                end if

                partial_stress = have_option(trim(u%option_path)//&
                    &"/prognostic/spatial_discretisation"//&
                    &"/discontinuous_galerkin/viscosity_scheme"//&
                    &"/partial_stress_form")
                ewrite(2,*) 'partial stress? ', partial_stress


                ! Options parsing for Large Eddy Simulation (LES).
                have_les = have_option(trim(u%option_path)//&
                    "/prognostic/spatial_discretisation"//&
                    "/discontinuous_galerkin/les_model")

                if(have_les) then
                   if(.not. partial_stress) then
                      FLAbort("Need to enable partial_stress viscosity scheme for DG LES")
                   end if
                   
                   ! Are we using the isotropic grid SGS eddy viscosity,
                   ! or AMD (anisotropic grids)?

                    have_isotropic_les = &
                        have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/isotropic")


                    have_amd_les = &
                        have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/amd")


                    have_vreman_les = &
                        have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/vreman")

                    have_roman_les = &
                        have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/roman")



                    ! Will eventually use partial_stress as an indicator
                    !                    if(partial_stress) then

                    ! Extract scalar or tensor eddy field
                    if(have_isotropic_les .or. have_amd_les &
                         .or. have_vreman_les) then
                      
                        ewrite(1,*) "*** Scalar-based DG LES"
                        ! les eddy visc field - needs to be nullified if non-existent
                        nullify(tensor_eddy_visc)
                        eddy_visc => extract_scalar_field(state, "ScalarEddyViscosity", stat=stat)
                        ! Theoretically it should always work - unless someone
                        ! deletes the field in a checkpoint FLML
                        if (stat/=0) then
                            FLAbort("Can't do scalar DG LES without a ScalarEddyViscosity field")
                        else
                            ewrite(1,*) "Found ScalarEddyViscosity field"
                         end if

                      else
                         
                        ewrite(1,*) "*** Tensor-based DG LES"
                        ! les eddy visc field - needs to be nullified if non-existent
                        nullify(eddy_visc)
                        tensor_eddy_visc => extract_tensor_field(state, "TensorEddyViscosity", stat=stat)
                        ! Theoretically it should always work - unless someone
                        ! deletes the field in a checkpoint FLML
                        if (stat/=0) then
                            FLAbort("Can't do tensor DG LES without a ScalarEddyViscosity field")
                        else
                            ewrite(1,*) "Found TensorEddyViscosity field"
                        end if
                        
                    end if

                    ! Options not specific to either isotropic or anisotropic LES
                    call get_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/smagorinsky_coefficient", &
                        smagorinsky_coefficient)

                    prescribed_filter_width => extract_scalar_field(state, "FilterWidth", stat=stat)
                    if (stat/=0) then
                        nullify(prescribed_filter_width)
                    end if

                    ! Van Driest damping requested?
                    have_van_driest = have_option(trim(u%option_path)//&
                        &"/prognostic/spatial_discretisation"//&
                        &"/discontinuous_galerkin/les_model"//&
                        &"/van_driest_damping")

                    if(have_van_driest) then
                        distance_to_wall=> extract_scalar_field(state, "DistanceToWall", stat=stat)
                        if (stat/=0) then
                            FLAbort("Van Driest damping requested, but no DistanceToWall scalar field exists")
                        end if
                    else
                        nullify(distance_to_wall)
                    end if

                    y_plus_debug => extract_scalar_field(state, "YPlus", stat=stat)
                    if (stat/=0) then
                        nullify(y_plus_debug)
                    end if

                    les_filter_width_debug => extract_scalar_field(state, "DampedFilterWidth", stat=stat)
                    if (stat/=0) then
                        nullify(les_filter_width_debug)
                    end if

                end if
                !  end of les variables

                integrate_surfacetension_by_parts = have_option(trim(u%option_path)//&
                    &"/prognostic/tensor_field::SurfaceTension"//&
                    &"/diagnostic/integrate_by_parts")

                assert(has_faces(X%mesh))
                assert(has_faces(P%mesh))

                call zero(big_m)
                subcycle=.false.
                if(present(subcycle_m)) subcycle=.true.
                if(subcycle) then
                    call zero(subcycle_m)
                end if
                call zero(RHS)
    
                if(present(inverse_masslump) .and. lump_mass) then
                    call allocate(inverse_masslump, u%dim, u%mesh, "InverseLumpedMass")
                    call zero(inverse_masslump)
                end if
                if(present(inverse_mass) .and. .not. lump_mass) then
                    assert(u%mesh%continuity<0)
                    mass_sparsity=make_sparsity_dg_mass(u%mesh)

                    if (pressure_corrected_absorption .or. has_boundary_condition(u, "dirichlet")) then
                        ! the diagonal blocks are different
                        call allocate( inverse_mass, mass_sparsity, (/ u%dim, u%dim /), &
                            diagonal=.true., name="InverseMassMatrix")
                    else
                        ! diagonal blocks are the same and all point to the same memory
                        call allocate( inverse_mass, mass_sparsity, (/ u%dim, u%dim /), &
                            diagonal=.true., equal_diagonal_blocks=.true., name="InverseMassMatrix")
                    end if
                    ! Drop the extra reference to sparsity.
                    call deallocate(mass_sparsity)
                end if
    
                ! get bc type and values on entire surface mesh
                ! numbering of types, determined by ordering here, i.e.
                ! weakdirichlet=1, free_surface=2
                allocate(velocity_bc_type(U%dim, surface_element_count(U)))
                call get_entire_boundary_condition(U, (/ &
                    "weakdirichlet       ", &
                    "free_surface        ", &
                    "no_normal_flow      ", &
                    "turbine_flux_penalty", &
                    "turbine_flux_dg     " /), velocity_bc, velocity_bc_type)

                ! the turbine connectivity mesh is only needed if one of the boundaries is a turbine.
                if (any(velocity_bc_type==4) .or. any(velocity_bc_type==5)) then
                    turbine_conn_mesh=get_periodic_mesh(state, u%mesh)
                end if

                ! same for pressure
                allocate(pressure_bc_type(surface_element_count(P)))
                call get_entire_boundary_condition(P, (/ &
                    "weakdirichlet", &
                    "dirichlet    "/), pressure_bc, pressure_bc_type)
                have_pressure_bc = any(pressure_bc_type>0)

                if (have_wd_abs) then
                    if (.not. has_scalar_field(state, "WettingDryingAlpha")) then
                        FLExit("Wetting and drying needs the diagnostic field WettingDryingAlpha activated.")
                    end if
                    ! The alpha fields lives on the pressure mesh, but we need it on the velocity, so let's remap it.
                    wettingdrying_alpha => extract_scalar_field(state, "WettingDryingAlpha")
                    call allocate(alpha_u_field, u%mesh, "alpha_u")
                    call remap_field(wettingdrying_alpha, alpha_u_field)
                end if

                call profiler_tic(u, "element_loop-omp_overhead")

#ifdef _OPENMP
    num_threads = omp_get_max_threads()
#else 
                num_threads=1
#endif

                if (have_viscosity) then
                    call get_mesh_colouring(state, u%mesh, COLOURING_DG2, colours)
                else
                    call get_mesh_colouring(state, u%mesh, COLOURING_DG0, colours)
                end if
#ifdef _OPENMP
    cache_valid = prepopulate_transform_cache(X)
    assert(cache_valid)
    if (have_coriolis) then
       call set_coriolis_parameters
    end if
#endif
                call profiler_toc(u, "element_loop-omp_overhead")
    
                call profiler_tic(u, "element_loop")

                ! Do necessary LES calculations
                if(have_les) then
                    if(have_isotropic_les) then
                       call calc_dg_sgs_scalar_viscosity(state, x, u)
                       
                    elseif(have_amd_les) then
                       call calc_dg_sgs_amd_viscosity(state, x, u)
                       
                    elseif(have_vreman_les) then
                       call calc_dg_sgs_vreman_viscosity(state, x, u)
                       
                    elseif(have_roman_les) then
                       call calc_dg_sgs_roman_viscosity(state, x, u)
                       
                    else
                        FLExit("Error: unsupported LES model (choose standard Smagorinsky or AMD")
                    end if
                end if

                ! This will eventually have a switch statement.
                ! switch(sim_type)
                ! case(optimised_assembly_case)
                !   call optimised_DG_assembly_elements()
                ! else
                !   Run unoptimised case (See loop below)

                !----------------------------------------------------------------------
                ! Establish local shape functions
                !----------------------------------------------------------------------

                if(new_mesh_connectivity .or. .not. allocated(coloured_ele_lists)) then
                    if(allocated(coloured_ele_lists)) then
                        do clr=1, size(colours)
                            deallocate(coloured_ele_lists(clr)%list)
                        end do
                        deallocate(coloured_ele_lists)
                    else
                        allocate(coloured_ele_lists(size(colours)))
                        do clr=1, size(colours)
                            allocate(coloured_ele_lists(clr)%list(key_count(colours(clr))))
                        end do
                    end if

                    do clr=1, size(colours)
                        len = key_count(colours(clr))

                        do nnid=1, len
                            coloured_ele_lists(clr)%list(nnid)=fetch(colours(clr), nnid)
                        end do
                    end do
                end if


                if(U%dim==opDim .and. P%mesh%shape%degree==opPresDeg &
                    .and. have_viscosity .and. viscosity_scheme==CDG &
                    .and. ele_loc(U,1)==opNloc .and. ele_ngi(U,1)==opNgi ) then
                    print*, "Optimised DG assembly: Compact DG"

                    inner_t0 = mpi_wtime()

                    !$OMP PARALLEL DEFAULT(SHARED) &
                    !$OMP PRIVATE(clr, nnid, ele, len)

                    call allocateDGAssemblyArrays()

                    colour_loop: do clr = 1, size(colours)
                      len = key_count(colours(clr))

                     !$OMP DO SCHEDULE(STATIC)
                      element_loop: do nnid = 1, len
                       ele = fetch(colours(clr), nnid)

                        call construct_momentum_elements_dg_opt( &
                             ele, big_m, rhs, &
                             X, U, advecting_velocity, U_mesh, X_old, X_new, &
                             Source, Buoyancy, hb_density, hb_pressure, gravity, Abs, Viscosity, &
                             swe_bottom_drag, swe_u_nl, &
                             P, old_pressure, Rho, surfacetension, q_mesh, &
                             velocity_bc, velocity_bc_type, &
                             pressure_bc, pressure_bc_type, &
                             turbine_conn_mesh, on_sphere, depth, have_wd_abs, &
                             alpha_u_field, Abs_wd, vvr_sf, ib_min_grad, nvfrac, &
                             inverse_mass=inverse_mass, &
                             inverse_masslump=inverse_masslump, &
                             mass=mass, subcycle_m=subcycle_m, partial_stress=partial_stress, &
                             have_les=have_les, have_isotropic_les=have_isotropic_les, &
                             have_amd_les=have_amd_les, &
                             have_vreman_les=have_vreman_les, &
                             have_roman_les=have_roman_les, &
                             smagorinsky_coefficient=smagorinsky_coefficient, &
                             eddy_visc=eddy_visc, tensor_eddy_visc=tensor_eddy_visc, &
                             prescribed_filter_width=prescribed_filter_width, &
                             distance_to_wall=distance_to_wall, y_plus_debug=y_plus_debug, &
                             les_filter_width_debug=les_filter_width_debug, &
                             have_free_stab=have_free_stab, free_stab_param=free_stab_param )

                      end do element_loop
                      !$OMP END DO

                    end do colour_loop

                    call deallocateDGAssemblyArrays()

                    !$OMP END PARALLEL
                    
                    inner_t1 = mpi_wtime()

                 else
                    print*, "udim, opdim:", U%dim, opDim
                    print*, "ele_loc, opnNloc",  ele_loc(U,1), opNloc
                    print*, "ele_ngi, opNgi",  ele_ngi(U,1), opNgi
                    
                    FLExit("Non-optimised DG assembly no longer supported")

                end if


                call profiler_toc(u, "element_loop")

                if (have_wd_abs) then
                    ! the remapped field is not needed anymore.
                    call deallocate(alpha_u_field)
                    !  deallocate(alpha_u_field)
                    call deallocate(Abs_wd)
                end if

                if (present(inverse_masslump) .and. lump_mass) then
                    call apply_dirichlet_conditions_inverse_mass(inverse_masslump, u)
                    ewrite_minmax(inverse_masslump)
                end if
                if (present(inverse_mass) .and. .not. lump_mass) then
                    call apply_dirichlet_conditions_inverse_mass(inverse_mass, u)
                    ewrite_minmax(inverse_mass)
                end if
                ewrite_minmax(rhs)



                ! Drop the reference to the fields we may have made.

!                call deallocate(uGrad)
!                deallocate(uGrad)
                !    print*,"after deallocs()"

                call deallocate(Viscosity)
                call deallocate(Abs)
                call deallocate(Source)
                call deallocate(U_nl)
                call deallocate(velocity_bc)
                call deallocate(pressure_bc)
                deallocate(velocity_bc_type)
                deallocate(pressure_bc_type)
                call deallocate(surfacetension)
                call deallocate(buoyancy)
                call deallocate(gravity)
                if(multiphase) then
                    call deallocate(nvfrac)
                end if
                call deallocate(dummyscalar)
                deallocate(dummyscalar)
    
                ewrite(1, *) "Exiting construct_momentum_dg"

                call profiler_toc("construct_momentum_dg")

                t1 = mpi_wtime()

                assemble_dt = t1-t0
                inner_assemble_dt = inner_t1-inner_t0
                if (lastt > 1e-10) then
                    total_dt = t1 - lastt
                    percent_dg =  (assemble_dt/total_dt)*100.0
                    inner_percent_dg =  (inner_assemble_dt/total_dt)*100.0
                else
                    percent_dg = 0.0
                    inner_percent_dg = 0.0
                end if
                lastt = t1

                print*, "**** DG_time_spent_in_assemble:", assemble_dt
                print*, "**** DG_%_in_assemble:", percent_dg
                print*, "**** DG_inner_loop_time_spent_in_assemble:", inner_assemble_dt
                print*, "**** DG_inner_loop_%_in_assemble:", inner_percent_dg
                print*, "**** DG_time_since_last_call:", total_dt

    
            end subroutine construct_momentum_dg




subroutine subcycle_momentum_dg(u, mom_rhs, subcycle_m, inverse_mass, state)
    type(vector_field), intent(inout) :: u
    type(vector_field), intent(inout):: mom_rhs
    type(block_csr_matrix), intent(in):: subcycle_m, inverse_mass
    type(state_type), intent(inout):: state

    type(vector_field) :: u_sub, m_delta_u, delta_u
    type(scalar_field), pointer :: courant_number_field
    type(scalar_field) :: u_cpt
    real :: max_courant_number
    integer :: d, i, subcycles
    logical :: limit_slope

    ! Arrays used in limit_vb_opt
    type(vector_field) :: T_max, T_min !, T_limit
    type(mesh_type), pointer :: vertex_mesh

    ! Benchmarking stuff
    real (kind=8) :: t0, t1, subcycle_dt, total_dt, percent_subcycle
    real (kind=8), save :: lastt
    real (kind=8), external :: mpi_wtime

    logical :: run_optimal

    ! Benchmarking
    t0=mpi_wtime()

#ifdef USE_CTO
    ! Only works for 3D
    run_optimal=(u%dim==opDim)
#else
    run_optimal=.false.
#endif

    ! Thought for today -- can I make this second order?
    ! Kuzmin, J. Comp. Appl. Math., 2010 seem to think so.
    ! - Angus, 4th June 2019.

    ewrite(1,*) 'Inside subcycle_momentum_dg'

    !Always limit slope using VB limiter if subcycling
    !If we get suitable alternative limiter options we shall use them
    limit_slope = .true.

    call get_option(trim(u%option_path)//&
        &"/prognostic/temporal_discretisation"//&
        &"/discontinuous_galerkin/maximum_courant_number_per_subcycle",&
        &max_courant_number)
    courant_number_field => &
        extract_scalar_field(state, "DG_CourantNumber")
    call calculate_diagnostic_variable(state, &
        "DG_CourantNumber", &
        & courant_number_field)
    subcycles = ceiling( maxval(courant_number_field%val)&
        &/max_courant_number)
    call allmax(subcycles)
    ewrite(2,*) 'Number of subcycles: ', subcycles
    if (subcycles==0) return

    call allocate(u_sub, u%dim, u%mesh, "SubcycleU")
    u_sub%option_path = trim(u%option_path)
    call set(u_sub, u)

    ! aux. field to store increment between subcycles
    call allocate(delta_u, u%dim, u%mesh, "SubcycleDeltaU")
    ! aux. field that incrementally computes M (u^sub-u^n)/dt
    call allocate(m_delta_u, u%dim, u%mesh, "SubcycleMDeltaU")
    call zero(m_delta_u)

    ! Allocate fields for limit_vb_opt

    ! returns linear version of T%mesh (if T%mesh is periodic, so is vertex_mesh)


    if(run_optimal) then
        print*, "subcycle_momentum_dg: using optimised code"

        call find_linear_parent_mesh(state, u_sub%mesh, vertex_mesh)
        call allocate(T_max, u%dim, vertex_mesh, trim(u_sub%name)//"LimitMax")
        call allocate(T_min, u%dim, vertex_mesh, trim(u_sub%name)//"LimitMin")
    else
        print*, "subcycle_momentum_dg: unoptimised method"
    end if

    if(subcycles > 5) subcycles=5

    do i=1, subcycles
        if (limit_slope) then

            if(run_optimal) then
                ! filter wiggles from u
                call limit_vb_opt(u_sub)
            else
               ! filter wiggles from u
                do d =1, u%dim
                    u_cpt = extract_scalar_field_from_vector_field(u_sub,d)
                    call limit_vb(state, u_cpt)
                end do
            end if

        end if


        ! du = advection * u
        call mult(delta_u, subcycle_m, u_sub)
        ! M*du/dt = M*du/dt - advection * u
        call addto(m_delta_u, delta_u, scale=-1.0/subcycles)

        ! we're only interested in m_delta_u, so we may leave early:
        if (i==subcycles) exit

        ! du = m^(-1) du
        call dg_apply_mass(inverse_mass, delta_u)

        ! u = u - dt/s * du
        call addto(u_sub, delta_u, scale=-dt/subcycles)
        call halo_update(u_sub)


     ! strictly speaking we should have another halo_update here, but
      ! we can assume that the limiting inside halo 1 elements can be
      ! performed locally

    end do

    ewrite_minmax(delta_u)

    !update RHS of momentum equation

    ! here is the low-down:
    !
    ! This is what we get from construct_momentum_dg:
    !   big_m = M + dt*theta*K, where K are any terms not included in subcycling (viscosity, coriolis etc.)
    !   mom_rhs = f - K u^n
    ! This is what we want to solve:
    !   M (u^sub - u^n)/dt + A u^n = 0, assuming one subcycle here
    !   M (u^n+1 - u^sub)/dt + K u^n+theta = f
    ! The last eqn can be rewritten:
    !   M (u^n+1 - u^n)/dt - M (u^sub - u^n)/dt + K u^n + dt*theta*K (u^n+1-u^n)/dt = f
    ! i.o.w.:
    !   big_m (u^n+1 - u^n)/dt = f - K u^n + M (u^sub - u^n)/dt
    ! This means mom_rhs needs to have M (u^sub - u^n)/dt added in
    ! and the implicit big_m solve computes a du/dt starting from u^n and not u^sub!
    ! Therefor this sub doesn't actually change u,  but only adds in the explicit advection
    ! to the rhs of the mom eqn.

    call addto(mom_rhs, m_delta_u)

    call deallocate(m_delta_u)
    call deallocate(u_sub)
    call deallocate(delta_u)

    if(run_optimal) then
        call deallocate(T_max)
        call deallocate(T_min)
    end if

    t1 = mpi_wtime()

    subcycle_dt = t1-t0
    if (lastt > 1e-10) then
        total_dt = t1 - lastt
        percent_subcycle=  (subcycle_dt/total_dt)*100.0
    else
        percent_subcycle = 0.0
    end if
    lastt = t1

    print*, "**** subcycle: DG_time_spent_in:", subcycle_dt
    print*, "**** subcycle: DG_time_since_last_call:", total_dt
    print*, "**** subcycle: DG_%_in:", percent_subcycle

contains
    ! There is an optimised version of this in Momentum_DG.F90
    ! limit_vb is only ever called from there

#if USE_CTO
    subroutine limit_vb_opt(T)
        !Vertex-based (not Victoria Bitter) limiter from
        !Kuzmin, J. Comp. Appl. Math., 2010
        ! doi:10.1016/j.cam.2009.05.028
        type(vector_field), intent(inout) :: T
        !

        ! counters
        integer :: ele, node, i
        ! local numbers
        ! integer, dimension(:), pointer :: T_ele
        ! gradient scaling factor
        real, dimension(opDim) :: alpha
        ! Global node lists for elements on DG and vertex meshes
        integer, dimension(:), pointer :: T_ele, V_ele
        ! local field values. P1DG, so always NLOC in size
        real, dimension(opDim, opNloc) :: T_val, T_val_slope, T_val_min,T_val_max, T_val_minus_bar
        real, dimension(opDim) :: Tbar

        if (.not. element_degree(T%mesh, 1)==1 .or. continuity(T%mesh)>=0) then
            FLExit("The vertex based slope limiter only works for P1DG fields.")
        end if
 
        call set(T_max, (/-huge(0.0),-huge(0.0),-huge(0.0)/))
        call set(T_min, (/huge(0.0),huge(0.0),huge(0.0)/))

        ! For each vertex in the mesh store the min and max values of the P1DG
        ! nodes directly surrounding it
        ! The subtle bit is that T is on a DG mesh, whereas T_min and T_max
        ! are on a linear CG mesh. This means the same nodes on T_min &
        ! T_max will be visited several times. (Necessary for min/max values)

        do ele = 1, ele_count(T)
            T_ele => T%mesh%ndglno(opNloc*(ele-1)+1:opNloc*ele)
            T_val = T%val(:,T_ele)
            do concurrent (i=1:opDim)
                Tbar(i) = sum(T_val(i, :))/opNloc
            end do

            ! do maxes and mins
            V_ele => vertex_mesh%ndglno(opNloc*(ele-1)+1:opNloc*ele)
            T_val_max = T_max%val(:,V_ele)
            T_val_min = T_min%val(:,V_ele)

            do node = 1,opNloc
                do concurrent (i=1:opDim)
                    T_val_min(i, node) = min(T_val_min(i, node), Tbar(i))
                    T_val_max(i, node) = max(T_val_max(i, node), Tbar(i))
                end do
            end do

            T_min%val(:, V_ele) = T_val_min
            T_max%val(:, V_ele) = T_val_max
        end do

        ! now for each P1DG node make sure the field value is between the recorded vertex min and max
        ! this is done without changing the element average (Tbar)
        do ele = 1, ele_count(T)
            !Set slope factor to 1
            alpha = 1.

            ! Element node numbers (T mesh and vertex mesh), with Tbar etc.
            T_ele => T%mesh%ndglno(opNloc*(ele-1)+1:opNloc*ele)
            T_val = T%val(:,T_ele)
            do concurrent (i=1:opDim)
                Tbar(i) = sum(T_val(i, :))/opNloc
            end do

            do concurrent (node=1:opNloc)
                T_val_slope(:,node) = T_val(:,node) - Tbar
                T_val_minus_bar(:,node) = T_val(:,node) - Tbar
            end do

            V_ele => vertex_mesh%ndglno(opNloc*(ele-1)+1:opNloc*ele)
            T_val_max = T_max%val(:,V_ele)
            T_val_min = T_min%val(:,V_ele)

            !loop over nodes, adjust alpha
            do node = 1, opNloc
                do concurrent (i=1:opDim)
                    !check whether to use max or min, and avoid floating point algebra errors due to round-off and underflow
                    if(T_val(i,node)>Tbar(i)*(1.0+sign(1.0e-12,Tbar(i))) .and. T_val_minus_bar(i,node) > tiny(0.0)*1e10) then
                        alpha(i) = min(alpha(i),(T_val_max(i,node)-Tbar(i))/T_val_minus_bar(i,node) )
                    else if(T_val(i,node)<Tbar(i)*(1.0-sign(1.0e-12,Tbar(i))) .and. T_val_minus_bar(i,node)  < -tiny(0.0)*1e10) then
                        alpha(i) = min(alpha(i),(T_val_min(i,node)-Tbar(i))/T_val_minus_bar(i,node) )
                    end if
                end do

            !call set(T, T_ele, Tbar + alpha*T_val_slope)
            T%val(:,T_ele(node))=Tbar+alpha*T_val_slope(:,node)

            end do
        end do

    end subroutine limit_vb_opt
#endif

end subroutine subcycle_momentum_dg

    
! The Coordinate and Solution fields of a turbine simulation live on a non-periodic mesh (that is with option remove-periodicity).
! This function takes such a field's mesh and returns the periodic mesh from which it is derived.
recursive function get_periodic_mesh(state, mesh) result(periodic_mesh)
    type(state_type), intent(in) :: state
    type(mesh_type), intent(in) :: mesh
    type(mesh_type) :: periodic_mesh
    character(len=OPTION_PATH_LEN) :: option_path
    character(len=4096) :: derived_meshname
    integer :: stat

    option_path=mesh%option_path
    if (have_option(trim(mesh%option_path) // '/from_mesh')) then
        call get_option(trim(mesh%option_path) // '/from_mesh/mesh/name', derived_meshname, stat)
        assert(stat==0)
        if (have_option(trim(mesh%option_path) // '/from_mesh/periodic_boundary_conditions/remove_periodicity')) then
            periodic_mesh=extract_mesh(state, derived_meshname, stat)
        else
            periodic_mesh=get_periodic_mesh(state, extract_mesh(state, derived_meshname, stat))
        end if
        assert(stat==0)
    else
        FLExit("A periodic mesh with remove_periodicity has to be used in combination with the turbine model.")
    end if
end function get_periodic_mesh

subroutine allocate_big_m_dg(state, big_m, u)
    !!< This routine allocates big_m as a petsc_csr_matrix without explicitly
    !!< constructing a sparsity, but only working the number of local and non-local
    !!< nonzero entries per row. As this should be a reasonably cheap operation this
    !!< is done every non-linear iteration.
    !!< Assumptions:
    !!< - contiguous numbering of owned nodes and elements
    !!< - number of nodes per element is the same
    !!< - both test and trial space are discontinuous
    type(state_type) :: state
    type(petsc_csr_matrix), intent(out):: big_m
    type(vector_field), intent(in):: u

    !! NOTE: use_element_blocks only works if all element have the same number of nodes
    logical:: use_element_blocks
      
    character(len=FIELD_NAME_LEN):: pc
    type(halo_type), pointer:: halo
    integer, dimension(:), pointer:: neighbours, neighbours2, nodes
    integer, dimension(:), allocatable:: dnnz, onnz
    logical:: compact_stencil, have_viscosity, have_coriolis, have_advection, have_turbine, partial_stress
    integer:: rows_per_dim, rows, nonods, elements
    integer:: owned_neighbours, foreign_neighbours, coupled_components, coupled_components_ele
    integer:: i, j, dim, ele, nloc
    type(mesh_type) :: neigh_mesh
      
    assert( continuity(u)<0 )
    
    compact_stencil = have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation"//&
        &"/discontinuous_galerkin/viscosity_scheme"//&
        &"/interior_penalty") .or. &
        &have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation"//&
        &"/discontinuous_galerkin/viscosity_scheme"//&
        &"/compact_discontinuous_galerkin")
                
    ! NOTE: this only sets the local have_viscosity, have_advection, have_coriolis and partial stress
    have_viscosity = have_option(trim(u%option_path)//&
        &"/prognostic/tensor_field::Viscosity")
    have_advection = .not. have_option(trim(u%option_path)//"/prognostic"//&
        &"/spatial_discretisation/discontinuous_galerkin"//&
        &"/advection_scheme/none")
    have_coriolis = have_option("/physical_parameters/coriolis")
    partial_stress = have_option(trim(u%option_path)//&
        &"/prognostic/spatial_discretisation"//&
        &"/discontinuous_galerkin/viscosity_scheme"//&
        &"/partial_stress_form")

    ! It would be enough to set this variable to true only if there is a flux turbine.
    ! However, for performance reasons, this is done whenever a turbine model is in use.
    have_turbine = have_option("/turbine_model")
    
    ! some preconditioners do not support petsc block matrix
    call get_option(trim(u%option_path)// &
        &"/prognostic/solver/preconditioner/name", pc)
    use_element_blocks = .not. (pc=="eisenstat" .or. pc=="mg" &
        .or. compact_stencil)

    if (have_turbine) then
        neigh_mesh=get_periodic_mesh(state, u%mesh)
    else
        neigh_mesh=u%mesh
    end if
    if (associated(u%mesh%halos)) then
        halo => u%mesh%halos(1)
        rows_per_dim=halo_nowned_nodes(halo)
    else
        nullify(halo)
        rows_per_dim=node_count(u)
    end if
    if (use_element_blocks) rows_per_dim=rows_per_dim/ele_loc(u,1)
    
    rows=rows_per_dim*u%dim
    allocate( dnnz(1:rows), onnz(1:rows) )
    
    coupled_components = 0
    coupled_components_ele = 0
    if (partial_stress) then
        coupled_components = u%dim - 1
    else if (have_coriolis) then
        coupled_components_ele = u%dim -1
    end if
    
    ! we first work everything out for rows corresponding to the first component
    do ele=1, element_count(u)
        ! we only have to provide nnz for owned rows. The owner
        ! therefore needs to specify the correct nnzs including
        ! contributions from others.
        ! NOTE: that the allocate interface assumes a contiguous
        ! numbering of owned nodes and elements
        if (.not. element_owned(u, ele)) cycle
      
        ! for each element work out the number of neighbours it talks to
      
        ! this is for zeroth order (i.e. without advection and viscosity)
        owned_neighbours = 0
        foreign_neighbours = 0
      
        if (have_viscosity .or. have_advection) then
            ! start with first order
            neighbours => ele_neigh(neigh_mesh, ele)
            do i=1, size(neighbours)
                ! skip boundaries
                if (neighbours(i)<=0) cycle
                if (element_owned(u, neighbours(i))) then
                    owned_neighbours = owned_neighbours+1
                else
                    foreign_neighbours = foreign_neighbours+1
                end if
            end do
        end if
      
        ! Added brackes around (.not. compact_stencil), check this
        if (have_viscosity .and. (.not. compact_stencil)) then
            ! traverse the second order neighbours
            do i=1, size(neighbours)
                ! skip boundaries
                if (neighbours(i)<=0) cycle
          
                neighbours2 => ele_neigh(neigh_mesh, neighbours(i))
                do j=1, size(neighbours2)
                    ! skip boundaries:
                    if (neighbours2(j)<=0) cycle
                    ! prevent double counting:
                    if (neighbours2(j)==ele .or. any(neighbours==neighbours2(j))) cycle

                    if (element_owned(u, neighbours2(j))) then
                        owned_neighbours = owned_neighbours + 1
                    else
                        foreign_neighbours = foreign_neighbours + 1
                    end if
                end do
            end do
        end if

        if (.not. use_element_blocks) then
            nodes => ele_nodes(u, ele)
            ! NOTE: there is an assumption here that n/o nodes of the neighbours
            ! is equal to that of ele (so in fact is the same for all elements)
            ! We need to do something more complicated if this is no longer true
            nloc = size(nodes)
            do i=1, nloc
                ! this break down as follows:
                ! 1                       for node-node coupling of the same component within the element
                ! owned_neighbours        for node-node coupling of the same component with 1st or 2nd order neighbours
                ! coupled components_ele  for node-node coupling with different components only within the element
                ! note: no coupling with different components of neighbouring elements as long as we're in tensor form
                ! coupled components      for node-node coupling with different components
                dnnz( nodes(i) ) = ( (1+owned_neighbours)*(coupled_components+1) + coupled_components_ele) * nloc
                ! this breaks down as follows:
                ! foreign_neighbours  for node-node coupling of the same component with neighbours that are owned by an other process
                ! note: coriolis only couples within the element and is therefore always completely local
                onnz( nodes(i) ) = foreign_neighbours*(coupled_components+1) * nloc
            end do
        else
            ! see above for reasoning
            dnnz(ele)=(1+owned_neighbours)*(coupled_components+1) + coupled_components_ele
            onnz(ele)=foreign_neighbours*(coupled_components+1)
        end if
    end do
      
    ! then copy to rows of other components
    do dim=2, u%dim
        dnnz( (dim-1)*rows_per_dim+1:dim*rows_per_dim ) = dnnz(1:rows_per_dim)
        onnz( (dim-1)*rows_per_dim+1:dim*rows_per_dim ) = onnz(1:rows_per_dim)
    end do
      
    if (use_element_blocks) then
        ! local owned and non-elements
        elements=element_count(u)
        call allocate(big_m, elements, elements, &
            dnnz, onnz, (/ u%dim, u%dim /), "BIG_m", halo=halo, &
            element_size=ele_loc(u,1))
    else
        ! local owned and non-owned nodes
        nonods=node_count(u)
        call allocate(big_m, nonods, nonods, &
            dnnz, onnz, (/ u%dim, u%dim /), "BIG_m", halo=halo)
    end if
      
end subroutine allocate_big_m_dg

subroutine correct_velocity_dg(U, inverse_mass, CT, delta_P)
    !!< Given the pressure correction delta_P, correct the velocity.
    !!<
    !!< U_new = U_old + M^{-1} * C * delta_P
    type(vector_field), intent(inout) :: U
    type(block_csr_matrix), intent(in):: inverse_mass
    type(block_csr_matrix), intent(in) :: CT
    type(scalar_field), intent(in) :: delta_P
    
    ! Correction to U one dimension at a time.
    type(scalar_field) :: delta_U1, delta_U2
    
    integer :: dim

    ewrite(1,*) 'correct_velocity_dg'

    call allocate(delta_U1, U%mesh, "Delta_U1")
    call allocate(delta_U2, U%mesh, "Delta_U2")
    
    do dim=1,U%dim

        call mult_T(delta_U1, block(CT,1,dim), delta_P)
        call mult(delta_U2, block(inverse_mass,dim, dim), delta_U1)

        call addto(U, dim, delta_U2)
      
    end do

    call halo_update(u)
    ewrite_minmax(u)

    call deallocate(delta_U1)
    call deallocate(delta_U2)

end subroutine correct_velocity_dg
    
subroutine assemble_poisson_rhs_dg(poisson_rhs, ctp_m, inverse_mass, &
    mom_rhs, ct_rhs, velocity, dt, theta_pg)

    type(scalar_field), intent(inout) :: poisson_rhs
    type(block_csr_matrix), intent(in) :: ctp_m
    type(block_csr_matrix), intent(in) :: inverse_mass
    type(vector_field), intent(inout) :: mom_rhs
    type(scalar_field), intent(inout) :: ct_rhs
    type(vector_field), intent(inout) :: velocity
    real, intent(in) :: dt, theta_pg

    type(vector_field) :: l_mom_rhs, minv_mom_rhs
    type(halo_type), pointer :: halo

    ewrite(1,*) 'Entering assemble_poisson_rhs_dg'
    
    ! poisson_rhs = ct_rhs/dt - C^T ( M^-1 mom_rhs + velocity/dt )

    if (IsParallel()) then

        call allocate(l_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssemblePoissonMomRHS")
        call set(l_mom_rhs, mom_rhs)
      
        ! we need to still add up the non-owned contributions from the global assembly of the mom_rhs
        ! this is done via a slight hack: assemble it as a petsc vector where petsc will add up the local
        ! contributions, and copy it back again
        halo => mom_rhs%mesh%halos(1)
        call addup_global_assembly(l_mom_rhs, halo)
      
    else
    
        l_mom_rhs =  mom_rhs

    end if
    
    ! compute M^-1 mom_rhs
    call allocate(minv_mom_rhs, mom_rhs%dim, mom_rhs%mesh, name="AssembleMinvPoissonMomRHS")
    call mult(minv_mom_rhs, inverse_mass, l_mom_rhs)
    call halo_update(minv_mom_rhs)
      
    call addto(minv_mom_rhs, velocity, scale=1.0/dt/theta_pg)
    call mult(poisson_rhs, ctp_m, minv_mom_rhs)

    call scale(poisson_rhs, -1.0)
    
    call addto(poisson_rhs, ct_rhs, scale=1.0/dt/theta_pg)

    call deallocate(minv_mom_rhs)
    if (IsParallel()) then
        call deallocate(l_mom_rhs)
    end if
    
    ewrite_minmax(poisson_rhs%val(1:nowned_nodes(poisson_rhs)))

end subroutine assemble_poisson_rhs_dg

subroutine momentum_DG_check_options
    
    character(len=OPTION_PATH_LEN) :: phase_path, velocity_path, dg_path
    integer :: i
    integer :: nstates ! number of states

    nstates=option_count("/material_phase")
    
    state_loop: do i=0, nstates-1

        phase_path="/material_phase["//int2str(i)//"]"
        velocity_path=trim(phase_path)//"/vector_field::Velocity/prognostic"
        dg_path=trim(velocity_path)//"/spatial_discretisation/discontinuous_galerkin"
       
        if (have_option(dg_path)) then
            if (have_option(trim(velocity_path)//"/solver/iterative_method::cg") &
                &.and. &
                &(  (.not. have_option(trim(dg_path)//"/advection_scheme/none")) &
                &    .or. have_option("/physical_parameters/coriolis"))) then
            
                ewrite(0,*) "Warning: You have selected conjugate gradient &
                &as a solver for"                
                ewrite(0,*) "    "//trim(phase_path)//&
                    &"/vector_field::Velocity"
                ewrite(0,*) "which is probably an asymmetric matrix"
            end if
        end if

        if (((have_option(trim(velocity_path)//"vertical_stabilization/vertical_velocity_relaxation") .or. &
            have_option(trim(velocity_path)//"vertical_stabilization/implicit_buoyancy")).and. &
            have_option(trim(velocity_path)//"vector_field::Absorption")) .and. &
            (.not. have_option(trim(velocity_path)//"vector_field::Absorption/include_pressure_correction"))) then
            ewrite(0,*) "Warning: You have selected a vertical stabilization but have not set"
            ewrite(0,*) "include_pressure_correction under your absorption field."
            ewrite(0,*) "This option will now be turned on by default."
        end if

        if (have_option(trim(dg_path)//"/viscosity_scheme/partial_stress_form") .and. .not. &
            (have_option(trim(dg_path)//"/viscosity_scheme/bassi_rebay") .or. &
            have_option(trim(dg_path)//"/viscosity_scheme/compact_discontinuous_galerkin") )) then
            FLAbort("partial stress form is only implemented for the bassi-rebay and CDG viscosity schemes in DG")
        end if

    end do state_loop

end subroutine momentum_DG_check_options

! Generated by script
subroutine allocateDGAssemblyArrays()
    allocate( Coriolis_mat(opNloc, opNloc) )
    allocate( rho_mat(opNloc, opNloc) )
    allocate( rho_move_mat(opNloc, opNloc) )
    allocate( mass_mat(opNloc, opNloc) )
    allocate( inverse_mass_mat(opNloc, opNloc) )
    allocate( Advection_mat(opNloc, opNloc) )
    allocate( Source_mat(opNloc, opNloc) )
    allocate( ele2grad_mat(opDim, opNloc, opNloc) )
    allocate( Abs_mat(opDim, opNloc, opNloc) )
    allocate( Abs_mat_sphere(opDim, opDim, opNloc, opNloc) )
    allocate( Abs_lump(opDim, opNloc) )
    allocate( Abs_lump_sphere(opDim, opDim, opNloc) )
    allocate( source_lump(opNloc) )
    allocate( Q_inv(opNloc, opNloc) )
    allocate( Grad_u_mat_q(opDim, opNloc, opEFloc) )
    allocate( Div_u_mat_q(opDim, opNloc, opEFloc) )
    allocate( Viscosity_mat(opDim, opDim, opEFloc, opEFloc) )
    allocate( node_stress_diag(opDim) )
    allocate( resid_stress_term(opDim) )

    allocate( big_m_diag_addto(opDim, opEFloc) )
    allocate( rhs_addto(opDim, opEFloc) )
    allocate( big_m_tensor_addto(opDim, opDim, opEFloc, opEFloc) )
    allocate( diagonal_block_mask(opDim, opDim) )
    allocate( off_diagonal_block_mask(opDim, opDim) )
    allocate( subcycle_m_tensor_addto(opDim, opDim, opEFloc, opEFloc) )
    allocate( sh_tensout(opDim, opDim, opNloc, opNloc) )
    allocate( sh_outtens(opDim, opDim, opNloc, opNloc) )
    allocate( sh_to_temp(opDim, opNloc) )
    allocate( sh_ot_temp(opNloc, opDim) )
    allocate( sh_dt_temp(opNloc, opNloc) )

    allocate( Viscosity_ele(opDim, opDim, opNloc) )
    allocate( visc_ele_quad(opDim, opDim, opNgi) )
    allocate( x_val(opDim, opNloc) )
    allocate( x_val_2(opDim, opNloc) )
    allocate( u_val(opDim, opNloc) )
    allocate( kappa_mat(opDim, opDim, opNloc, opNloc) )
    allocate( l_MassLump(opNloc) )
    allocate( l_move_masslump(opNloc) )
    allocate( local_glno(opEFloc) )
    allocate( detwei(opNgi) )
    allocate( detwei_old(opNgi) )
    allocate( detwei_new(opNgi) )
    allocate( coefficient_detwei(opNgi) )
    allocate( detwei_rhoq(opNgi) )
    allocate( du_t(opNloc, opNgi, opDim) )
    allocate( dug_t(opNloc, opNgi, opDim) )
    allocate( dq_t(opNloc, opNgi, opDim) )
    allocate( Rho_q(opNgi) )
    allocate( Coriolis_q(opNgi) )
    allocate( u_nl_q(opDim, opNgi) )
    allocate( u_nl_div_q(opNgi) )
    allocate( tension(opDim, opDim, opNgi) )
    allocate( dtensiondj(opDim, opNgi) )

    allocate( absorption_gi(opDim, opNgi) )
    allocate( tensor_absorption_gi(opDim, opDim, opNgi) )
    allocate( vvr_abs(opDim, opDim, opNgi) )
    allocate( vvr_abs_diag(opDim, opNgi) )
    allocate( depth_at_quads(opNgi) )
    allocate( ib_abs(opDim, opDim, opNgi) )
    allocate( ib_abs_diag(opDim, opNgi) )
    allocate( dt_rho(opNloc, opNgi, opDim) )
    allocate( grav_at_quads(opDim, opNgi) )
    allocate( grad_rho(opDim, opNgi) )
    allocate( ele_grav_val(opDim, opNloc) )
    allocate( drho_dz(opNgi) )

    allocate( ele_u_mesh_quad(opDim, opNgi) )
    allocate( ele_centre(opDim) )
    allocate( neigh_centre(opDim) )
    allocate( face_centre(opDim) )
    allocate( face_centre_2(opDim) )
    allocate( alpha_u_quad(opNgi) )

    allocate( face_Rho_q(opFngi) )
    allocate( face_Rho_val(opFloc) )
    allocate( face_normal(opDim, opFngi) )
    allocate( face_u_nl_q(opDim, opFngi) )
    allocate( face_u_f_q(opDim, opFngi) )
    allocate( face_u_f2_q(opDim, opFngi) )
    allocate( face_div_u_f_q(opDim, opFngi) )
    allocate( face_u_mesh_quad(opDim, opFngi) )
    allocate( face_inflow(opFngi) )
    allocate( face_u_nl_q_dotn(opFngi) )
    allocate( face_income(opFngi) )
    allocate( face_detwei(opFngi) )
    allocate( face_detwei_work(opFngi) )
    allocate( face_inner_advection_integral(opFngi) )
    allocate( face_outer_advection_integral(opFngi) )
    allocate( face_nnAdvection_out(opFloc, opFloc) )
    allocate( face_nnAdvection_in(opFloc, opFloc) )
    allocate( face_mnCT(1, opDim, opPFloc, opFloc) )
    allocate( face_kappa_gi(opDim, opDim, opFngi) )
    allocate( face_visc_val(opDim, opDim, opFloc) )
    allocate( face_tension_q(opDim, opDim, opFngi) )

    allocate( tmp_face_tensor(opDim, opDim, opFloc) )
    allocate( face_primal_fluxes_mat(2, opFloc, opNloc) )
    allocate( face_shape_shape_work(opFloc, opFloc) )
    allocate( face_penalty_fluxes_mat(2, opFloc, opFloc) )
    allocate( face_normal_mat(opDim, opFloc, opFloc) )
    allocate( face_kappa_normal_mat(opDim, opFloc, opFloc) )

    allocate( matmul_dut_visc(opDim) )

end subroutine allocateDGAssemblyArrays

! Generated by script
subroutine deallocateDGAssemblyArrays()
    deallocate( Coriolis_mat )
    deallocate( rho_mat )
    deallocate( rho_move_mat )
    deallocate( mass_mat )
    deallocate( inverse_mass_mat )
    deallocate( Advection_mat )
    deallocate( Source_mat )
    deallocate( ele2grad_mat )
    deallocate( Abs_mat )
    deallocate( Abs_mat_sphere )
    deallocate( Abs_lump )
    deallocate( Abs_lump_sphere )
    deallocate( source_lump )
    deallocate( Q_inv )
    deallocate( Grad_u_mat_q )
    deallocate( Div_u_mat_q )
    deallocate( Viscosity_mat )
    deallocate( node_stress_diag )
    deallocate( resid_stress_term )

    deallocate( big_m_diag_addto )
    deallocate( rhs_addto )
    deallocate( big_m_tensor_addto )
    deallocate( diagonal_block_mask )
    deallocate( off_diagonal_block_mask )
    deallocate( subcycle_m_tensor_addto )
    deallocate( sh_tensout )
    deallocate( sh_outtens )
    deallocate( sh_to_temp )
    deallocate( sh_ot_temp )
    deallocate( sh_dt_temp )

    deallocate( Viscosity_ele )
    deallocate( visc_ele_quad )
    deallocate( x_val )
    deallocate( x_val_2 )
    deallocate( u_val )
    deallocate( kappa_mat )
    deallocate( l_MassLump )
    deallocate( l_move_masslump )
    deallocate( local_glno )
    deallocate( detwei )
    deallocate( detwei_old )
    deallocate( detwei_new )
    deallocate( coefficient_detwei )
    deallocate( detwei_rhoq )
    deallocate( du_t )
    deallocate( dug_t )
    deallocate( dq_t )
    deallocate( Rho_q )
    deallocate( Coriolis_q )
    deallocate( u_nl_q )
    deallocate( u_nl_div_q )
    deallocate( tension )
    deallocate( dtensiondj )

    deallocate( absorption_gi )
    deallocate( tensor_absorption_gi )
    deallocate( vvr_abs )
    deallocate( vvr_abs_diag )
    deallocate( depth_at_quads )
    deallocate( ib_abs )
    deallocate( ib_abs_diag )
    deallocate( dt_rho )
    deallocate( grav_at_quads )
    deallocate( grad_rho )
    deallocate( ele_grav_val )
    deallocate( drho_dz )

    deallocate( ele_u_mesh_quad )
    deallocate( ele_centre )
    deallocate( neigh_centre )
    deallocate( face_centre )
    deallocate( face_centre_2 )
    deallocate( alpha_u_quad )

    deallocate( face_Rho_q )
    deallocate( face_Rho_val )
    deallocate( face_normal )
    deallocate( face_u_nl_q )
    deallocate( face_u_f_q )
    deallocate( face_u_f2_q )
    deallocate( face_div_u_f_q )
    deallocate( face_u_mesh_quad )
    deallocate( face_inflow )
    deallocate( face_u_nl_q_dotn )
    deallocate( face_income )
    deallocate( face_detwei )
    deallocate( face_detwei_work )
    deallocate( face_inner_advection_integral )
    deallocate( face_outer_advection_integral )
    deallocate( face_nnAdvection_out )
    deallocate( face_nnAdvection_in )
    deallocate( face_mnCT )
    deallocate( face_kappa_gi )
    deallocate( face_visc_val )
    deallocate( face_tension_q )

    deallocate( tmp_face_tensor )
    deallocate( face_primal_fluxes_mat )
    deallocate( face_shape_shape_work )
    deallocate( face_penalty_fluxes_mat )
    deallocate( face_normal_mat )
    deallocate( face_kappa_normal_mat )

    deallocate( matmul_dut_visc )
    
  end subroutine deallocateDGAssemblyArrays


end module momentum_DG

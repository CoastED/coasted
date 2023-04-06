!    Copyright (C) 2023 Dr. Angus Creech and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
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

#define MESH_DIM 3
#define NLOC 4

module momentum_FV
  ! This module contains the Finite Volume form of the momentum  equation. 
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
       COLOURING_DG0
  use coriolis_module
  use halos
  use sparsity_patterns
  use petsc_tools
  use turbine
  use diagnostic_fields
  use slope_limiters_dg
  use smoothing_module
  use fields_manipulation
  use field_options
  use sparsity_patterns_meshes
  use colouring
  use Profiler
#ifdef _OPENMP
  use omp_lib
#endif
  use multiphase_module

  use les_module

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

contains

  subroutine construct_momentum_fv(u, p, rho, x, &
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
    type(vector_field), pointer :: U_mesh, X_old, X_new
    type(vector_field), target :: U_nl
    !! Projected velocity field for them as needs it. 
    type(vector_field), target :: pvelocity
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
    type(vector_field) :: velocity_bc
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

    ! Partial stress - sp911
    logical :: partial_stress

    ! LES - sp911
    logical :: have_les = .false., have_van_driest = .false.
    real :: smagorinsky_coefficient
    type(scalar_field), pointer :: eddy_visc, prescribed_filter_width, distance_to_wall, &
         & y_plus_debug, les_filter_width_debug
    type(tensor_field), pointer :: tensor_eddy_visc


    ! Benchmarking stuff
    real :: t0, t1, assemble_dt, total_dt, percent_dg
    real, save :: lastt
    real, external :: mpi_wtime


    t0 = mpi_wtime()

    ewrite(1, *) "In construct_momentum_fv"

    call profiler_tic("construct_momentum_fv")
    assert(continuity(u)<0)

    acceleration= .not. present_and_false(acceleration_form)
    ewrite(2, *) "Acceleration form? ", acceleration


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

       if(have_option(trim(U%option_path)//"/prognostic"//&
            &"/spatial_discretisation/discontinuous_galerkin"//&
            &"/advection_scheme/project_velocity_to_continuous")) then
          ewrite(3,*) 'CREATING PROJECTEDNONLINEARVELOCITY, cjc'
          if(.not.has_scalar_field(state, "ProjectedNonlinearVelocity")) then
          
             call get_option(trim(U%option_path)//"/prognostic"//&
                  &"/spatial_discretisation/discontinuous_galerkin"//&
                  &"/advection_scheme/project_velocity_to_continuous"//&
                  &"/mesh/name",pmesh_name)
             pmesh = extract_mesh(state, pmesh_name)
             call allocate(pvelocity, U_nl%dim, pmesh, &
                  &"ProjectedNonlinearVelocity")
             call project_field(U_nl, pvelocity, X)
             call insert(state, pvelocity, "ProjectedNonlinearVelocity")
             advecting_velocity => pvelocity

             ! Discard the additional reference.
             call deallocate(pvelocity)
          else
             pvelocity = extract_vector_field(state, &
                  &"ProjectedNonlinearVelocity")

             advecting_velocity => pvelocity
          end if
       else
          advecting_velocity => U_nl
       end if
       have_advection = .true.
    else
       ! Forcing a zero NonlinearVelocity will disable advection.
       call allocate(U_nl, MESH_DIM,  U%mesh, "NonlinearVelocity", &
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
       call allocate(Source, MESH_DIM,  U%mesh, "VelocitySource", FIELD_TYPE_CONSTANT)
       call zero(Source)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Source)
       ewrite_minmax(source)
    end if

    Abs=extract_vector_field(state, "VelocityAbsorption", stat)   
    have_absorption = (stat==0)
    if (.not.have_absorption) then
       call allocate(Abs, MESH_DIM, U%mesh, "VelocityAbsorption", FIELD_TYPE_CONSTANT)
       call zero(Abs)
    else
       ! Grab an extra reference to cause the deallocate below to be safe.
       call incref(Abs)
       ewrite_minmax(Abs)
    end if

    have_wd_abs=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/dry_absorption")
    ! Absorption term in dry zones for wetting and drying
    if (have_wd_abs) then
       call allocate(Abs_wd, MESH_DIM, U%mesh, "VelocityAbsorption_WettingDrying", FIELD_TYPE_CONSTANT)
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
         &discontinuous_galerkin/advection_scheme/integrate_advection_by_parts/once")

    integrate_conservation_term_by_parts = have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/advection_scheme/integrate_conservation_term_by_parts")

    ! Determine the scheme to use to discretise viscosity.
    if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme/bassi_rebay")) then
       viscosity_scheme=BASSI_REBAY
    else if (have_option(trim(U%option_path)//"/prognostic/spatial_discretisation/&
         &discontinuous_galerkin/viscosity_scheme&
         &/compact_discontinuous_galerkin")) then
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

    have_van_driest = .false.
    have_les = have_option(trim(u%option_path)//&
            "/prognostic/spatial_discretisation"//&
            "/discontinuous_galerkin/les_model")

    if( have_les) then
          if(partial_stress) then
              call get_option(trim(u%option_path)//&
                 &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/les_model"//&
                &"/smagorinsky_coefficient", &
                smagorinsky_coefficient)

              have_van_driest = have_option(trim(u%option_path)//&
                 &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/les_model"//&
                &"/van_driest_damping")

            ! les variables - need to be nullified if non-existent
            eddy_visc => extract_scalar_field(state, "DGLESScalarEddyViscosity", stat=stat)
            if (stat/=0) then
                nullify(eddy_visc)
            end if
            prescribed_filter_width => extract_scalar_field(state, "FilterWidth", stat=stat)
            if (stat/=0) then
                nullify(prescribed_filter_width)
            end if
        else
            ewrite(1,*) "*** Tensor-based DG LES (experimental)"

              call get_option(trim(u%option_path)//&
                 &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/les_model"//&
                &"/smagorinsky_coefficient", &
                smagorinsky_coefficient)

              have_van_driest = have_option(trim(u%option_path)//&
                 &"/prognostic/spatial_discretisation"//&
                &"/discontinuous_galerkin/les_model"//&
                &"/van_driest_damping")

            ! les variables - need to be nullified if non-existent
               tensor_eddy_visc => extract_tensor_field(state, "DGLESTensorEddyViscosity", stat=stat)
                if (stat/=0) then
                    ewrite(1,*) "*** Found DGLESTensorEddyViscosity field"
                    nullify(tensor_eddy_visc)
                end if

            if (stat/=0) then
                nullify(prescribed_filter_width)
            end if
        end if
    end if


    distance_to_wall => extract_scalar_field(state, "DistanceToWall", stat=stat)
    if (stat/=0) then
      if(have_van_driest) then
        FLAbort("You must create a prescribed scalar field DistanceToWall to use Van Driest damping with LES")
      end if
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
    allocate(velocity_bc_type(MESH_DIM, surface_element_count(U)))
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

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(clr, nnid, ele, len)

    colour_loop: do clr = 1, size(colours) 
      len = key_count(colours(clr))

      !$OMP DO SCHEDULE(STATIC)
      element_loop: do nnid = 1, len
       ele = fetch(colours(clr), nnid)
       call construct_momentum_element_dg( ele, big_m, rhs, &
            & X, U, advecting_velocity, U_mesh, X_old, X_new, &
            & Source, Buoyancy, hb_density, hb_pressure, gravity, Abs, Viscosity, &
            & swe_bottom_drag, swe_u_nl, &
            & P, old_pressure, Rho, surfacetension, q_mesh, &
            & velocity_bc, velocity_bc_type, &
            & pressure_bc, pressure_bc_type, &
            & turbine_conn_mesh, on_sphere, depth, have_wd_abs, &
            & alpha_u_field, Abs_wd, vvr_sf, ib_min_grad, nvfrac, &
            & inverse_mass=inverse_mass, &
            & inverse_masslump=inverse_masslump, &
            & mass=mass, subcycle_m=subcycle_m, partial_stress=partial_stress, &
            have_les=have_les, smagorinsky_coefficient=smagorinsky_coefficient, &
            eddy_visc=eddy_visc, tensor_eddy_visc=tensor_eddy_visc, &
            prescribed_filter_width=prescribed_filter_width, &
            distance_to_wall=distance_to_wall, y_plus_debug=y_plus_debug, &
            les_filter_width_debug=les_filter_width_debug )

      end do element_loop
      !$OMP END DO

    end do colour_loop
    !$OMP END PARALLEL

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

    if (associated(eddy_visc)) then
      ! eddy visc is calculated in momentum_dg element loop. we need to do a halo_update
      call halo_update(eddy_visc)
    end if

    if (associated(tensor_eddy_visc)) then
      ! tensor eddy visc is calculated in momentum_dg element loop. we need to do a halo_update
      call halo_update(tensor_eddy_visc)
    end if


    ! Drop the reference to the fields we may have made.
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
    
    ewrite(1, *) "Exiting construct_momentum_fv"

    call profiler_toc("construct_momentum_fv")

    t1 = mpi_wtime()

    assemble_dt = t1-t0
    if (lastt > 1e-10) then
        total_dt = t1 - lastt
        percent_dg =  (assemble_dt/total_dt)*100.0
    else
        percent_dg = 0.0
    end if
    lastt = t1

    print*, "**** FV - time spent in assemble:", assemble_dt
    print*, "**** FV - time since last call:", total_dt
    print*, "**** FV - % in assemble:", percent_dg

    
  end subroutine construct_momentum_fv

  subroutine momentum_fv_check_options
    
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
            have_option(trim(dg_path)//"/viscosity_scheme/bassi_rebay")) then
         FLAbort("partial stress form is only implemented for the bassi-rebay viscosity scheme in DG")
       end if

    end do state_loop

  end subroutine momentum_fv_check_options



end module momentum_fv

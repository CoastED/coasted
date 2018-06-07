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

! This module contains subroutines for creating useful option defaults,
! eg. projected CG velocity (VelocityCG) and DG_CourantNumber fields.
! Note, these defaults can always be overridden by the user creating the
! fields themselves in Diamond - as was the case before.

module DG_prep

    use global_parameters
    use futils
    use spud
    use state_module

    implicit none

    private
    public :: check_and_add_DG_defaults

contains

    ! ========================================================================
    ! The main routine and entry point
    ! ========================================================================

    subroutine check_and_add_DG_defaults(state)
        type(state_type), dimension(:), pointer :: state

        integer :: ph
        character(len=OPTION_PATH_LEN) :: phase_path, dg_path

        ! Loop over phases. If any DG, optionally set up default fields and
        ! associated options so the user doesn't have to.
        do ph = 1, size(state)

            phase_path="/material_phase["//int2str(ph-1)//"]/"

            dg_path=trim(phase_path)&
                    //"vector_field::Velocity/prognostic/spatial_discretisation/discontinuous_galerkin/"

            ! If we have a Discontinuous Galerkin field
            if(have_option(trim(dg_path))) then

                ! Create Velocity CG vector field options in FLML tree using DG velocity
                ! and create DG_CourantNumber options if necessary
                call create_velocity_cg_options(ph)

                ! If we have adaptive timestepping, create a DG_CourantNumber
                ! field if we haven't already
                if(have_option("/timestepping/adaptive_timestep")) then
                    call create_dg_courant_field_options(ph)
                end if

                if(have_option(trim(dg_path)//"les_model")) then
                    call create_dg_les_field_options(phase_path, dg_path)
                end if

            end if


        end do
    end subroutine check_and_add_DG_defaults

    ! ========================================================================
    ! Create VelocityCG field options for given phase
    ! ========================================================================

    subroutine create_velocity_cg_options(ph)
        integer :: ph

        character(len=OPTION_PATH_LEN) :: opt_path
        integer :: stat

        opt_path="/material_phase["//int2str(ph-1)//"]/vector_field::VelocityCG/"

        ! Only create if we don't have this field already. The user can
        ! can override by creating their own field

        if(.not. have_option(trim(opt_path))) then

            ewrite(1,*) "Creating VelocityCG field"

            call add_option(trim(opt_path), stat)
            call set_option_attribute(trim(opt_path)//"rank", "1", stat)
            call add_option(trim(opt_path)//"diagnostic", stat)

            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/name", &
                "vector_galerkin_projection", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/material_phase_support", &
                "single", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/source_field_type", &
                "vector", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/source_field_name", &
                "Velocity", stat)

            ! Creating generic solver options

            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/solver/iterative_method/name", &
                "cg", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/solver/preconditioner/name", &
                "sor", stat)
            call set_option(trim(opt_path)//"diagnostic/algorithm/solver/relative_error", &
                1.0e-7, stat)
            call set_option(trim(opt_path)//"diagnostic/algorithm/solver/max_iterations", &
                30, stat)
            call add_option(trim(opt_path)//"diagnostic/algorithm/solver/never_ignore_solver_failures", stat)
            call add_option(trim(opt_path)//"diagnostic/algorithm/solver/diagnostics", stat)
            call add_option(trim(opt_path)//"diagnostic/algorithm/solver/diagnostics/monitors", stat)

            ! Diagnostic outputs, etc.

            call add_option(trim(opt_path)//"diagnostic/mesh/name", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/mesh/name", &
                "CoordinateMesh", stat)

            call add_option(trim(opt_path)//"diagnostic/output", stat)
            call add_option(trim(opt_path)//"diagnostic/stat", stat)
            call add_option(trim(opt_path)//"diagnostic/stat/include_in_stat", stat)
            call add_option(trim(opt_path)//"diagnostic/convergence", stat)
            call add_option(trim(opt_path)//"diagnostic/convergence/exclude_from_convergence", stat)
            call add_option(trim(opt_path)//"diagnostic/detectors", stat)
            call add_option(trim(opt_path)//"diagnostic/detectors/include_in_detectors", stat)
            call add_option(trim(opt_path)//"diagnostic/steady_state", stat)
            call add_option(trim(opt_path)//"diagnostic/steady_state/include_in_steady_state", stat)
            call add_option(trim(opt_path)//"diagnostic/consistent_interpolation", stat)
        end if

    end subroutine create_velocity_cg_options


    ! ========================================================================
    ! Create DG_CourantNumber field for a given phase
    ! ========================================================================

    subroutine create_dg_courant_field_options(ph)
        integer :: ph

        character(len=OPTION_PATH_LEN) :: opt_path
        integer :: stat

        opt_path="/material_phase["//int2str(ph-1)//"]/scalar_field::DG_CourantNumber/"

        ! Only create if we don't have this field already. The user can
        ! can override by creating their own field
        if(.not. have_option(trim(opt_path))) then

            ewrite(1,*) "Creating DG_CourantNumber field"

            call add_option(trim(opt_path), stat)
            call set_option_attribute(trim(opt_path)//"rank", "0", stat)
            call add_option(trim(opt_path)//"diagnostic", stat)

            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/name", &
                "Internal", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/algorithm/material_phase_support", &
                "multiple", stat)

            call add_option(trim(opt_path)//"diagnostic/mesh/name", stat)
            call set_option_attribute(trim(opt_path)//"diagnostic/mesh/name", &
                "VelocityMesh", stat)

            call add_option(trim(opt_path)//"diagnostic/output", stat)
            call add_option(trim(opt_path)//"diagnostic/output/exclude_from_vtu", stat)
            call add_option(trim(opt_path)//"diagnostic/stat", stat)
            call add_option(trim(opt_path)//"diagnostic/convergence", stat)
            call add_option(trim(opt_path)//"diagnostic/convergence/include_in_convergence", stat)
            call add_option(trim(opt_path)//"diagnostic/detectors", stat)
            call add_option(trim(opt_path)//"diagnostic/detectors/include_in_detectors", stat)
            call add_option(trim(opt_path)//"diagnostic/steady_state", stat)
            call add_option(trim(opt_path)//"diagnostic/steady_state/include_in_steady_state", stat)
        end if

    end subroutine create_dg_courant_field_options

    ! ========================================================================
    !
    ! ========================================================================
    subroutine create_dg_les_field_options(phase_path, dg_path)
        character(len=OPTION_PATH_LEN) :: phase_path, dg_path, &
            tensor_eddy_visc_path, scalar_eddy_visc_path

        logical :: have_les_option, have_les_visc_field, have_isotropic_les, have_partial_stress
        integer :: stat

        scalar_eddy_visc_path = trim(phase_path)//"scalar_field::ScalarEddyViscosity/"
        tensor_eddy_visc_path = trim(phase_path)//"tensor_field::TensorEddyViscosity/"

        have_isotropic_les = have_option(trim(dg_path)//"les_model/isotropic")
        have_partial_stress = have_option(trim(dg_path)//"viscosity_scheme/partial_stress_form")

        have_les_option = have_option(trim(dg_path)//"les_model")
        have_les_visc_field = (have_option(trim(scalar_eddy_visc_path)) &
            .or. have_option(trim(tensor_eddy_visc_path)))

        ! Partial stress needs to be on for any kind of turbulence modelling.
!        if(.not. have_partial_stress) then
!            FLExit("Error. Partial stress must be selected for Large Eddy Simulation")
!        end if


        if(have_les_option .and. .not. have_les_visc_field) then
            if(have_isotropic_les) then
                ! Create SGS Eddy Viscosity scalar field
                ewrite(1,*) "Creating ScalarEddyViscosity field"

                call add_option(trim(scalar_eddy_visc_path), stat)
                call set_option_attribute(trim(scalar_eddy_visc_path)//"rank", "0", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic", stat)

                call set_option_attribute(trim(scalar_eddy_visc_path)//"diagnostic/algorithm/name", &
                    "Internal", stat)
                call set_option_attribute(trim(scalar_eddy_visc_path)//"diagnostic/algorithm/material_phase_support", &
                    "single", stat)

                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/mesh/name", stat)
#ifdef LES_USES_DG_VEL
                call set_option_attribute(trim(scalar_eddy_visc_path)//"diagnostic/mesh/name", &
                    "VelocityMesh", stat)
#else
                call set_option_attribute(trim(scalar_eddy_visc_path)//"diagnostic/mesh/name", &
                    "CoordinateMesh", stat)
#endif
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/output", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/stat", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/convergence", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/convergence/exclude_from_convergence", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/detectors", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/detectors/include_in_detectors", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/steady_state", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/steady_state/include_in_steady_state", stat)
                call add_option(trim(scalar_eddy_visc_path)//"diagnostic/consistent_interpolation", stat)

            else
                ! Create tensor Eddy Viscosity field
                ewrite(1,*) "Creating TensorEddyViscosity field"

                call add_option(trim(tensor_eddy_visc_path), stat)
                call set_option_attribute(trim(tensor_eddy_visc_path)//"rank", "2", stat)
                call add_option(trim(tensor_eddy_visc_path)//"diagnostic", stat)

                call set_option_attribute(trim(tensor_eddy_visc_path)//"diagnostic/algorithm/name", &
                    "Internal", stat)
                call set_option_attribute(trim(tensor_eddy_visc_path)//"diagnostic/algorithm/material_phase_support", &
                    "single", stat)

                call add_option(trim(tensor_eddy_visc_path)//"diagnostic/mesh/name", stat)
                call set_option_attribute(trim(tensor_eddy_visc_path)//"diagnostic/mesh/name", &
                    "CoordinateMesh", stat)
                call add_option(trim(tensor_eddy_visc_path)//"diagnostic/output", stat)
                call add_option(trim(tensor_eddy_visc_path)//"diagnostic/stat", stat)
                call add_option(trim(tensor_eddy_visc_path)//"diagnostic/stat/include_in_stat", stat)
                call add_option(trim(tensor_eddy_visc_path)//"diagnostic/consistent_interpolation", stat)
            end if
        end if

    end subroutine

end module DG_prep

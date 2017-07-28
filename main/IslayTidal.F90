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

module islay_tidal

  use AuxilaryOptions
  use MeshDiagnostics
  use signal_vars
  use spud
  use equation_of_state
  use timers
  use adapt_state_module
  use adapt_state_prescribed_module
  use FLDebug
  use sparse_tools
  use elements
  use fields
  use boundary_conditions_from_options
  use populate_state_module
  use populate_sub_state_module
  use reserve_state_module
  use vtk_interfaces
  use Diagnostic_variables
  use diagnostic_fields_new, only : &
    & calculate_diagnostic_variables_new => calculate_diagnostic_variables, &
    & check_diagnostic_dependencies
  use diagnostic_fields_wrapper
  use diagnostic_children
  use advection_diffusion_cg
  use advection_diffusion_DG
  use advection_diffusion_FV
  use field_equations_cv, only: solve_field_eqn_cv, initialise_advection_convergence, coupled_cv_field_eqn
  use vertical_extrapolation_module
  use qmesh_module
  use checkpoint
  use write_state_module
  use synthetic_bc
  use goals
  use adaptive_timestepping
  use conformity_measurement
  ! Use Solid-fluid coupling and ALE - Julian- 18-09-06
  use ale_module
  use adjacency_lists
  use multimaterial_module
  use parallel_tools
  use SolidConfiguration

  use biology
  use momentum_equation
  use timeloop_utilities
  use field_priority_lists
  use boundary_conditions
  use spontaneous_potentials, only: calculate_electrical_potential
  use saturation_distribution_search_hookejeeves
  use discrete_properties_module

  use iceshelf_meltrate_surf_normal
  use halos
  use memory_diagnostics
  use free_surface_module
  use global_parameters, only: current_time, dt, timestep, OPTION_PATH_LEN, &
                               simulation_start_time, &
                               simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name
  use eventcounter

  use multiphase_module
  use detector_parallel, only: sync_detector_coordinates, deallocate_detector_list_array
  use momentum_diagnostic_fields, only: calculate_densities
  use sediment_diagnostics, only: calculate_sediment_flux

  implicit none

  private
        real, parameter :: minx=0, miny=0
        real, parameter :: maxx=21453.7278466, maxy=36651.1955983
        real, parameter :: sizex=21453.7278466, sizey=36651.1955983

        real, parameter :: maxlat=56.0230, minlat=55.6886
        real, parameter :: minlon=-6.2398, maxlon=-5.8951
        real, parameter :: sizelat = (maxlat-minlat), sizelon=(maxlon-minlon)

        integer, parameter :: ndim = 3

        real, parameter :: thickness = 1000

        type BoundaryPoint
            real :: x, y
            real :: m2amp, m2phase, s2amp, s2phase, n2amp, n2phase
            real :: k2amp, k2phase, k1amp, k1phase, o1amp, o1phase
        end type


  public set_islay_boundary_absorption, set_islay_boundary_nonhydrostatic_pressure


  contains
        subroutine set_islay_boundary_absorption(state)
        type(state_type), dimension(:), pointer :: state
        type(Vector_Field) :: velAbs, position
        type(Scalar_Field) :: distTop, distBottom

        integer :: i, j, ncnodes, nvnodes
        logical :: doNotRecalculate
        real :: relaxtime

        real :: absorb, wt, wt2, x, y, depth
        integer, save :: alreadyran
        real, allocatable, save :: absarray(:)
        real, parameter :: relaxhours = 0.1
!        real, parameter :: relaxhours = 0.5
!        real, parameter :: shoreHours=0.25
!        real, parameter :: minShoreDepth = 5, maxShoreDepth=10
!        real :: shoreAbs, shoreScale


        if(.not. has_vector_field(state(1), "VelocityAbsorption")) then
            FLAbort("Error: you must have a diagnostic velocity absorption field (internal algorithm) on the Coordinate mesh to use set_islay_boundary_absorption()")
        end if

        velAbs     = extract_vector_field(state(1), "VelocityAbsorption")
        position   = extract_vector_field(state(1), "Coordinate")
        distTop    = extract_scalar_field(state(1), "DistanceToTop")
        distBottom = extract_scalar_field(state(1), "DistanceToBottom")


        print*, "*** Setting sponge layer absorption for Islay"

        relaxtime = 60.*60.*relaxhours
        ncnodes = position%mesh%nodes
        nvnodes = velAbs%mesh%nodes

        if(ncnodes .ne. nvnodes) then
            FLAbort("set_islay_boundary_absorption(): VelocityAbsorption must be on Coordinate mesh")
        end if

        ! If we don't want to recalculate, then don't.
        if(have_option(trim(velAbs%option_path)//"/diagnostic/do_not_recalculate") &
            .and. allocated(absarray)) then

           velAbs%val(1, :) = absarray(:)
           velAbs%val(2, :) = absarray(:)
           velAbs%val(3, :) = absarray(:)

           return
        end if

        if( .not. allocated(absarray) ) then
            allocate(absarray(ncnodes))
        end if

        do i=1, ncnodes
            depth = distBottom%val(i)+distTop%val(i)
!             if(depth <= minShoreDepth) then
!               shoreAbs = 1.0/(shoreHours*60.0*60.0)
!             elseif(depth > maxShoreDepth) then
!                shoreAbs=0.0
!             else
!                shoreScale=1.0-(depth-minShoreDepth)/(maxShoreDepth-minShoreDepth)
!                shoreAbs = shoreScale/(shoreHours*60.0*60.0)
!            end if

            absorb = 0.

            wt = 0.
            wt2 = 0.

            x = position%val(1, i)
            y = position%val(2, i)

            if(x < minx+thickness) then
                wt=1.0-(x-minx)/thickness

            elseif(x > maxx-thickness) then
                wt=1.0-(maxx-x)/thickness

            end if

            if(y < miny+thickness) then
                wt2=1.0-(y-miny)/thickness
            elseif(y > maxy-thickness) then
                wt2=1.0-(maxy-y)/thickness
            end if

            if(wt2>wt) wt=wt2

            absorb=wt/relaxtime

 !           if(shoreAbs > absorb) absorb=shoreAbs

            velAbs%val(1, i) = absorb
            velAbs%val(2, i) = absorb
            velAbs%val(3, i) = absorb
            absarray(i) = absorb
        end do

        alreadyran = 1

    end subroutine set_islay_boundary_absorption


    subroutine init_pressure_boundary_data(northBC, westBC, southBC, eastBC)
        type(BoundaryPoint), allocatable, dimension(:) :: northBC, westBC, southBC, eastBC

        integer :: northct, westct, southct, eastct
        integer, parameter :: fd =1103

        integer :: err, i, nlines
        character(len=1000) :: lineText
        real :: lon, lat, x, y, scalex, scaley
        real :: m2amp, m2phase, s2amp, s2phase, n2amp, n2phase, k2amp, k2phase
        real :: k1amp, k1phase, o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

        open(fd, file="tidalconst/islay_tidal_const.csv", access="sequential", iostat=err)

        ! How many lines in the file?
        nlines=0

        northct=0
        westct=0
        southct=0
        eastct=0

        ! Ignore header row
        read(fd, '(A)', end=999) lineText
        do
            read(fd, *, end=999) lon, lat, m2amp, m2phase, s2amp, s2phase, &
                n2amp, n2phase, k2amp, k2phase, k1amp, k1phase, &
                o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

            ! Northern boundary
            if(abs(lat-maxlat) < 0.01) northct=northct+1
            ! Western boundary
            if(abs(lon-minlon) < 0.01) westct=westct+1
            ! Southern boundary
            if(abs(lat-minlat) < 0.01) southct=southct+1
            ! Eastern boundary
            if(abs(lon-maxlon) < 0.01) eastct=eastct+1

        end do


999  continue


        ! now we know how to allocate
        allocate( northBC(northct), southBC(southct), westBC(westct), eastBC(eastct) )

        rewind(fd)

        ! Read first line in again.
        read(fd, '(A)', end=999) lineText

        scalex=sizelon/sizex
        scaley=sizelat/sizey

        northct=0
        westct=0
        southct=0
        eastct=0

        do i=1, nlines-1
            read(fd, *) lon, lat, m2amp, m2phase, s2amp, s2phase, &
                n2amp, n2phase, k2amp, k2phase, k1amp, k1phase, &
                o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

            x = (lon-minlon) / scalex
            y = (lat-minlat) / scaley

            ! Northern boundary
            if(abs(lat-maxlat) < 0.01) then
                northBC(northct+1)=BoundaryPoint(x, y, m2amp, m2phase, &
                    s2amp, s2phase, n2amp, n2phase, k2amp, k2phase, &
                    k1amp, k1phase, o1amp, o1phase)
                northct=northct+1
            end if

            ! Western boundary
            if(abs(lon-minlon) < 0.01) then
                westBC(westct+1)=BoundaryPoint(x, y, m2amp, m2phase, &
                    s2amp, s2phase, n2amp, n2phase, k2amp, k2phase, &
                    k1amp, k1phase, o1amp, o1phase)
                westct=westct+1
            end if

            ! Southern boundary
            if(abs(lat-minlat) < 0.01) then
                southBC(southct+1)=BoundaryPoint(x, y, m2amp, m2phase, &
                    s2amp, s2phase, n2amp, n2phase, k2amp, k2phase, &
                    k1amp, k1phase, o1amp, o1phase)
                southct=southct+1
            end if

            ! Eastern boundary
            if(abs(lon-maxlon) < 0.01) then
                eastBC(eastct+1)=BoundaryPoint(x, y, m2amp, m2phase, &
                    s2amp, s2phase, n2amp, n2phase, k2amp, k2phase, &
                    k1amp, k1phase, o1amp, o1phase)
                eastct=eastct+1
            end if
        end do

        close(fd)



    end subroutine init_pressure_boundary_data

    subroutine  set_islay_boundary_nonhydrostatic_pressure(state)
        type(state_type), dimension(:), pointer :: state
        type(Scalar_Field) :: pressure
        type(Vector_Field) :: pos, remap
        type(BoundaryPoint), allocatable, dimension(:), save :: northBC, westBC, southBC, eastBC

        integer :: nnodes, i
        integer, save :: readConstituents

        real :: lon, lat, x(3)
        real :: m2amp, m2phase, s2amp, s2phase, n2amp, n2phase, k2amp, k2phase
        real :: k1amp, k1phase, o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

        print*, "*** set_islay_boundary_nonhydrostatic_pressure()"

        if(readConstituents==0) then
            call init_pressure_boundary_data(northBC, westBC, southBC, eastBC)
            readConstituents=1
        end if

        pressure = extract_scalar_field(state(1), "Pressure")
        pos = extract_vector_field(state(1), "Coordinate")

        ! Positions remapped to Pressure space
        call allocate(remap, 3, pressure%mesh, name="PressureCoordinate")
        call remap_field(pos, remap)

        nnodes = remap%mesh%nodes

        do i=1, nnodes
            x = remap%val(:,i)

        end do



    end subroutine  set_islay_boundary_nonhydrostatic_pressure

  end module islay_tidal

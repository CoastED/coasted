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

  use fetools
  use fields
  use field_options
  use field_derivatives
  use smoothing_module
  use vector_tools
  use state_module
  use state_fields_module

  use global_parameters, only: current_time, dt, timestep, OPTION_PATH_LEN, &
                               simulation_start_time, &
                               simulation_start_cpu_time, &
                               simulation_start_wall_time, &
                               topology_mesh_name, new_mesh_geometry
  use physics_from_options

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
        type(Vector_Field), pointer :: velAbs, position
        type(Scalar_Field), pointer :: distTop, distBottom

        integer :: i, j, ncnodes, nvnodes
        logical :: recalcAbs
        real :: relaxtime

        real :: absorb, wt, wt2, x, y, depth
        real, allocatable, save :: absarray(:)
        real, parameter :: relaxhours = 0.1
!        real, parameter :: relaxhours = 0.5
!        real, parameter :: shoreHours=0.25
!        real, parameter :: minShoreDepth = 5, maxShoreDepth=10
!        real :: shoreAbs, shoreScale

        if(.not. has_vector_field(state(1), "VelocityAbsorption")) then
            FLAbort("Error: you must have a diagnostic velocity absorption field (internal algorithm) on the Coordinate mesh to use set_islay_boundary_absorption()")
        end if

        velAbs     => extract_vector_field(state(1), "VelocityAbsorption")
        position   => extract_vector_field(state(1), "Coordinate")
!        distTop    => extract_scalar_field(state(1), "DistanceToTop")
!        distBottom => extract_scalar_field(state(1), "DistanceToBottom")


        print*, "*** Setting sponge layer absorption for Islay"

        relaxtime = 60.*60.*relaxhours
        ncnodes = position%mesh%nodes
        nvnodes = velAbs%mesh%nodes

        if(ncnodes .ne. nvnodes) then
            FLAbort("set_islay_boundary_absorption(): VelocityAbsorption must be on Coordinate mesh")
        end if

        recalcAbs=.false.
        if(new_mesh_geometry) then
            if( .not. allocated(absarray) ) then
                print*, "*** allocating absarray"
                allocate(absarray(ncnodes))
            else
                print*, "*** deallocating/reallocating absarray"
                deallocate(absarray)
                allocate(absarray(ncnodes))
            end if

            recalcAbs=.true.
        end if

        if(recalcAbs) then
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

                absarray(i) = absorb
            end do
        end if

        ! Copy to velocity absorption array.
        velAbs%val(1, :) = absarray(:)
        velAbs%val(2, :) = absarray(:)
        velAbs%val(3, :) = absarray(:)

    end subroutine set_islay_boundary_absorption


    subroutine init_pressure_boundary_data(tidalBC)
        type(BoundaryPoint), allocatable :: tidalBC(:)

        integer :: tidalct
        integer, parameter :: fd =1103

        integer :: err, i, nlines
        character(len=1000) :: lineText
        real :: lon, lat, x, y, scalex, scaley
        real :: m2amp, m2phase, s2amp, s2phase, n2amp, n2phase, k2amp, k2phase
        real :: k1amp, k1phase, o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

        print*, "init_pressure_boundary_data()"

        open(fd, file="tidalconst/boundaries.csv", access="sequential", iostat=err)

        ! How many lines in the file?
        tidalct=0

        ! Ignore header row
        read(fd, '(A)', end=999) lineText
        do
            read(fd, *, end=999) lon, lat, m2amp, m2phase, s2amp, s2phase, &
                n2amp, n2phase, k2amp, k2phase, k1amp, k1phase, &
                o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

            tidalct=tidalct+1
        end do

999  continue
        nlines=tidalct


        ! now we know how to allocate
        allocate( tidalBC(tidalct) )

        rewind(fd)

        ! Read first line in again.
        read(fd, '(A)', end=999) lineText

        scalex=sizelon/sizex
        scaley=sizelat/sizey

        tidalct=0

        do i=1, nlines-1
            read(fd, *) lon, lat, m2amp, m2phase, s2amp, s2phase, &
                n2amp, n2phase, k2amp, k2phase, k1amp, k1phase, &
                o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

            x = (lon-minlon) / scalex
            y = (lat-minlat) / scaley

            tidalBC(tidalct+1)=BoundaryPoint(x, y, m2amp, m2phase, &
                s2amp, s2phase, n2amp, n2phase, k2amp, k2phase, &
                k1amp, k1phase, o1amp, o1phase)

            tidalct=tidalct+1
        end do

        close(fd)

    end subroutine init_pressure_boundary_data



    subroutine  set_islay_boundary_nonhydrostatic_pressure(state, &
            surface_field, bc_position, bc_type_path, field_name)

        type(state_type) :: state
        type(scalar_field), intent(inout):: surface_field
        type(vector_field), intent(in):: bc_position
        character(len=*), intent(in):: bc_type_path, field_name

        type(Scalar_Field), pointer :: pressure
        type(Vector_Field) :: pos, remap
        type(BoundaryPoint), allocatable, dimension(:), save :: tidalpt
        type(BoundaryPoint) :: thisbp

        integer :: nnodes, i, j, nspecpoints
        integer, save :: readConstituents

        real :: lon, lat, x(3)
        real :: m2amp, m2phase, s2amp, s2phase, n2amp, n2phase, k2amp, k2phase
        ! real :: k1amp, k1phase, o1amp, o1phase, p1amp, p1phase, q1amp, q1phase

        real :: m2ang, s2ang, n2ang, k2ang
        real :: dist, wt, sumwt, rampTime, rampval, rho0, pval

        real :: t
        real, parameter :: pi=3.14159265359, piConv=2.0*pi/360.0, grav=9.81

        call get_fs_reference_density_from_options(rho0, state%option_path)

        call set(surface_field, 0.0)

        m2ang = piConv / 12.4206012
        s2ang = piConv / 12.0
        n2ang = piConv / 12.65834751
        k2ang = piConv / 11.96723606

        t = current_time

        ramptime=48*60*60

        print*, "*** set_islay_boundary_nonhydrostatic_pressure()"

        if(readConstituents==0) then
            call init_pressure_boundary_data(tidalpt)
            readConstituents=1
        end if

        pressure => extract_scalar_field(state, "Pressure")
        pos = extract_vector_field(state, "Coordinate")

        ! Positions remapped to Pressure space
        call allocate(remap, 3, pressure%mesh, name="PressureCoordinate")
        call remap_field(pos, remap)

        nnodes = node_count(bc_position)
        nspecpoints = size(tidalpt)

        do i=1, nnodes

            ! Only the four main tidal constituents for now
            m2amp=0
            m2phase=0
            s2amp=0
            s2phase=0
            n2amp=0
            n2phase=0
            k2amp=0
            k2phase=0

            x = node_val(bc_position, i)

            sumwt=0

            ! Calculate contributions to mesh points from each tidal boundary
            ! constituent point

            do j=1, nspecpoints
                dist = sqrt( (x(1)-tidalpt(j)%x)**2.0 + (x(2)-tidalpt(j)%y)**2.0 )
                wt=1/dist
                sumwt=sumwt+wt

                m2amp = m2amp + wt*tidalpt(j)%m2amp
                m2phase = m2phase+ wt*tidalpt(j)%m2phase * piConv
                s2amp = s2amp + wt*tidalpt(j)%s2amp
                s2phase = s2phase + wt*tidalpt(j)%s2phase * piConv

                n2amp = n2amp + wt*tidalpt(j)%n2amp
                n2phase = n2phase + wt*tidalpt(j)%n2phase * piConv
                k2amp = k2amp + wt*tidalpt(j)%k2amp
                k2phase = k2phase + wt*tidalpt(j)%k2phase * piConv
            end do

            ! now set pressure
            if(t<ramptime) then
                rampval=t/ramptime
            else
                rampval=1
            end if

            pval = rampval * sumwt * rho0 * grav &
                *  ( m2amp * cos(m2ang*t - m2phase) &
                + s2amp * cos(s2ang*t - s2phase) &
                + n2amp * cos(n2ang*t - n2phase) &
                + k2amp * cos(k2ang*t - k2phase) )

            call addto(surface_field, i, pval)

        end do


        call deallocate(remap)
    end subroutine  set_islay_boundary_nonhydrostatic_pressure

  end module islay_tidal

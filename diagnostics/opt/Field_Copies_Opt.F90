    subroutine gp_discontinuous_NLOC_NGI
      integer, parameter :: firstEle=1

      integer, dimension(:), pointer :: nodes
      type(element_type), pointer :: shape

      real, dimension(NLOC) :: little_rhs, ele_val
      real, dimension(NLOC, NLOC) :: little_mass
      !real, dimension(ele_ngi(s_field, firstEle)) :: detwei
      real, dimension(NGI) :: detwei

      integer :: ele, iloc, jloc, gi

      shape => ele_shape(s_field, ele)

      do ele = 1, ele_count(s_field)

        call transform_to_physical(positions, ele, detwei = detwei)

        !little_mass = shape_shape(shape, shape, detwei)

        ! little_mass shape_shape replacement
        little_mass=0.
        do jloc = 1, NLOC
            do iloc = 1, NLOC
                ! Main mass matrix.
                little_mass(iloc,jloc)=shape%n(iloc,1)*shape%n(jloc,1) * detwei(1)

                do gi = 2, NLOC
                    little_mass(iloc,jloc)=little_mass(iloc,jloc)+shape%n(iloc,gi)*shape%n(jloc,gi) * detwei(gi)
                end do
            end do
        end do

        little_rhs = shape_rhs(shape, detwei * ele_val_at_quad(source_field, ele))
        call solve(little_mass, little_rhs)

        nodes => ele_nodes(s_field, ele)

        call set(s_field, nodes, little_rhs)
      end do

    end subroutine gp_discontinuous_NLOC_NGI

#undef gp_discontinuous_NLOC_NGI

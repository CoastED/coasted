#include "confdefs.h"

! Wrapper for the ll2r3_rotate subroutine in femtools/Coordinates.F90
subroutine rotate_ll2cart(longitude, latitude, u, v, r3u, r3v, r3w) bind(c)
  use, intrinsic :: iso_c_binding
  use Coordinates
  implicit none
#ifdef DOUBLEP
  real(c_double), intent(in):: longitude, latitude, u, v
  real(c_double), intent(out):: r3u, r3v, r3w
#else
  real(c_float), intent(in):: longitude, latitude, u, v
  real(c_float), intent(out):: r3u, r3v, r3w
#endif

  call ll2r3_rotate(longitude, latitude, u, v, r3u, r3v, r3w)

end subroutine


/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* We do support the IEC 559 math functionality, real and complex.  */

/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */

/* We do not support C11 <threads.h>.  */

!WRF:DRIVER_LAYER:UTIL
!

SUBROUTINE wrf_message( str )
  IMPLICIT NONE

  CHARACTER*(*) str
  print *,trim(str)

END SUBROUTINE wrf_message

SUBROUTINE wrf_error_fatal( str )
  IMPLICIT NONE
  CHARACTER*(*) str
  CALL wrf_message( '-------------- FATAL CALLED ---------------' )
  CALL wrf_message( str )
  CALL wrf_message( '-------------------------------------------' )

  CALL wrf_abort
END SUBROUTINE wrf_error_fatal

SUBROUTINE wrf_abort
      STOP 'wrf_abort'
END SUBROUTINE wrf_abort

SUBROUTINE wrf_debug( level , str ) 
  IMPLICIT NONE 
  CHARACTER*(*) str 
  INTEGER , INTENT (IN) :: level 
  CALL wrf_message( str ) 
  RETURN 
END SUBROUTINE wrf_debug 

subroutine wrf_dm_bcast_real(rval, ival)
  implicit none
  real,    intent(in) :: rval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_real

subroutine wrf_dm_bcast_integer(rval, ival)
  implicit none
  integer, intent(in) :: rval
  integer, intent(in) :: ival
end subroutine wrf_dm_bcast_integer

subroutine wrf_dm_bcast_string(rval, ival)
  implicit none
  character(len=*), intent(in) :: rval
  integer,          intent(in) :: ival
end subroutine wrf_dm_bcast_string

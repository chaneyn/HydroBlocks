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

module kwm_string_utilities

contains

  subroutine strrep(string, rep1, rep2)
    ! In string 1, replace expression rep1 with expression rep2.
    ! if rep2 is shorter than rep1, fill the end of the string
    ! with blanks

    implicit none
    character(len=*) :: string, rep1, rep2
    integer :: idx, inlen, len1, len2

    do
       inlen = len(string)
       len1 = len(rep1)
       len2 = len(rep2)
       idx = index(string, rep1)-1

       if (idx == -1) then
          return
       else
          string = string(1:idx)// rep2 // &
               string((idx+len1+1):inlen) // &
               "                                                    "
       endif
    enddo

  end subroutine strrep


  character(len=1024) function unblank(string) result(return_string)
    ! Remove blanks (and tabs [char(9)] from string
    implicit none
    character(len=*), intent(in) :: string
    integer :: i, j,  lenstr

    return_string = " "

    if ( verify(string," "//char(9)) == 0 ) then
       stop 'String is all blanks.'
    endif

    j = 0
    do i = 1, len(string)
       if ((string(i:i).ne.' ').and.(string(i:i).ne.char(9))) then
          j = j + 1
          return_string(j:j) = string(i:i)
       endif
    enddo
    
  end function unblank

!KWM  character function upcase(h)
!KWM    implicit none
!KWM    character :: h
!KWM
!KWM    if ((ichar(h).ge.96) .and. (ichar(h).le.123)) then
!KWM       upcase = char(ichar(h)-32)
!KWM    else
!KWM       upcase = h
!KWM    endif
!KWM
!KWM  end function upcase

  character(len=256) function upcase(h) result(return_string)
    implicit none
    character(len=*), intent(in) :: h
    integer :: i
    
    return_string = " "

    do i = 1, len_trim(h)

       if ((ichar(h(i:i)).ge.96) .and. (ichar(h(i:i)).le.123)) then
          return_string(i:i) = char(ichar(h(i:i))-32)
       else
          return_string(i:i) = h(i:i)
       endif
    enddo

  end function upcase

!KWM  character function downcase(h)
!KWM    implicit none
!KWM    character h
!KWM
!KWM    if ((ichar(h).ge.65) .and. (ichar(h).le.90)) then
!KWM       downcase = char(ichar(h)+32)
!KWM    else
!KWM       downcase = h
!KWM    endif
!KWM  end function downcase

  character(len=256) function downcase(h) result(return_string)
    implicit none
    character(len=*), intent(in) :: h
    integer :: i

    return_string = " "

    do i = 1, len_trim(h)

       if ((ichar(h(i:i)).ge.65) .and. (ichar(h(i:i)).le.90)) then
          return_string(i:i) = char(ichar(h(i:i))+32)
       else
          return_string(i:i) = h(i:i)
       endif
    enddo

  end function downcase


end module kwm_string_utilities

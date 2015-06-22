subroutine initialize(wcolumns,wrowindex,nhsu,nvalues)

 implicit none
 integer*4,intent(in) :: nhsu,nvalues
 integer*4,intent(inout),dimension(nvalues) :: wcolumns
 integer*4,intent(inout),dimension(nhsu+1) :: wrowindex
 call initialize_dss(wcolumns,wrowindex,nhsu,nvalues)

end subroutine

subroutine finalize()

 implicit none
 call finalize_dss()

end subroutine

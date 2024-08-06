subroutine time_average(series_in,series_out,nt_out,nt_in)

 integer,intent(in) :: nt_out,nt_in
 real*8,intent(inout) :: series_in(nt_in)
 real*8,intent(out) :: series_out(nt_out)
 integer :: i
 integer :: dt
 dt = nt_in/nt_out
 do i = 1,nt_out
  series_out(i) = sum(series_in((1+(i-1)*dt):(i*dt)))/dt
 enddo

endsubroutine

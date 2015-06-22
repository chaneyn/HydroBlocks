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

subroutine update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,&
                  area,dx,dt,celerity,celerity1,storage_mask,wvalues,wcolumns,&
                  wrowindex,nthreads,maxntt,isw,nhsu,nvalues)

 implicit none
 integer*4 :: nhsu,nvalues
 real*8,intent(in) :: dt
 real*8,intent(in) :: isw
 integer*8,intent(in),dimension(nhsu) :: storage_mask
 real*8,intent(in),dimension(nhsu) :: dx,area,celerity,celerity1
 real*8,intent(inout),dimension(nhsu) :: recharge,storage,recharge1,storage1,qout,qin,qout1,qin1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*8,intent(in),dimension(nvalues) :: wvalues
 integer*4,intent(in) :: nthreads,maxntt
 real*8,dimension(nhsu) :: scarea
 real*8 :: dt_minimum,dtt
 integer*8 :: ntt,itime

 !Define the specific catchment area
 scarea = area/dx

 !Determine the required time step
 dt_minimum = minval(dx/celerity)
 ntt = int(ceiling(dt/dt_minimum))
 if (ntt .gt. maxntt) ntt = maxntt
 if (ntt .le. 0) ntt = 1
 dtt = dt/ntt

 !Process the smaller time steps
 do itime = 1,ntt
 
   !Solve the kinematic wave for this time step
   call solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                             recharge1,scarea,dx,dtt,celerity,celerity1,&
                             wvalues,wcolumns,wrowindex,isw,storage_mask)

 enddo

end subroutine update

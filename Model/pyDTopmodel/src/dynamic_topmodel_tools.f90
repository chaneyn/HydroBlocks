subroutine update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,&
                  area,dx,dt,celerity,celerity1,wvalues,wcolumns,wrowindex,&
                  nthreads,nhsu,nvalues)

 use omp_lib
 implicit none
 integer :: nhsu,nvalues
 real*8,intent(in) :: dt
 real*4,intent(in),dimension(nhsu) :: dx,area,celerity,celerity1
 real*4,intent(inout),dimension(nhsu) :: recharge,storage,recharge1,storage1,qout,qin,qout1,qin1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer,intent(in) :: nthreads
 real*4,dimension(nhsu) :: scarea
 real*4 :: dt_minimum,dtt
 integer :: ntt,itime

 !Define the specific catchment area
 scarea = area/dx

 !Determine the required time step
 dt_minimum = minval(dx/celerity)
 ntt = int(ceiling(dt/dt_minimum))
 if (ntt .le. 0) ntt = 1
 dtt = dt/ntt

 !Set the number of threads
 call mkl_set_num_threads(nthreads)
          
 !Process the smaller time steps
 do itime = 1,ntt
 
   !Solve the kinematic wave for this time step
   call solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                             recharge1,scarea,dx,dtt,celerity,celerity1,&
                             wvalues,wcolumns,wrowindex)

 enddo

end subroutine update

subroutine solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                                recharge1,scarea,dx,dtt,celerity,celerity1,&
                                wvalues,wcolumns,wrowindex)

 implicit none
 integer,intent(in) :: nhsu,nvalues
 real*4,intent(in) :: dtt
 real*4,intent(inout),dimension(nhsu) :: storage,qout,qin,recharge,storage1,qout1,qin1,recharge1
 real*4,intent(in),dimension(nhsu) :: scarea,dx,celerity,celerity1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 real*4,dimension(nvalues) :: A
 real*4,dimension(nhsu) :: denominator,numerator1,numerator2,part1,part2
 real*4 :: w
 integer*4 :: i

 !Define implicit scheme parameters
 w = 0.5 

 !Calculate as many constants as possible before going into the iterations
 denominator = (1 + dtt*w*celerity/dx)
 numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
 numerator2 = dtt*w*celerity/dx
 part1 = numerator1/denominator
 part2 = numerator2/denominator

 !Element wise multiplication on the rows
 do i = 1,nhsu
  A(wrowindex(i)+1:wrowindex(i+1)) = wvalues(wrowindex(i)+1:wrowindex(i+1))*scarea(i)*&
                                     part2(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)/&
                                     scarea(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)
 enddo

 !Solve for this time step
 !qout[:] = scipy.sparse.linalg.spsolve((I-A).T,b)

 !Set all negative fluxes to 0
 where (qout < 0.0) qout = 0.0

 !Calculate qin
 call mkl_cspblas_scsrgemv('T',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 qin = qin/scarea

 !Adjust the storages
 storage = storage + dtt*((qin - qout)/dx + recharge)

 !Set the next time step's info
 qout1 = qout
 qin1 = qin

end subroutine




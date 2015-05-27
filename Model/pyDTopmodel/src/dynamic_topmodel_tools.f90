subroutine update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,&
                  area,dx,dt,nhsu,celerity,celerity1,wvalues,wcolumns,wrowindex,nvalues,nthreads)

 use omp_lib
 implicit none
 integer :: nhsu,nvalues
 real*8,intent(in) :: dt
 !real*4,intent(in),dimension(nhsu,nhsu) :: weights
 real*4,intent(in),dimension(nhsu) :: dx,area,celerity,celerity1
 real*4,intent(inout),dimension(nhsu) :: recharge,storage,recharge1,storage1,qout,qin,qout1,qin1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer,intent(in) :: nthreads
 real*4,dimension(nhsu) :: storage_,qout_,storage1_,qout1_,qin1_
 real*8 :: dt_minimum,dtt
 integer :: ntt,itime

 !Determine the required time step
 dt_minimum = minval(dx/celerity)
 ntt = 2*(int(ceiling(dt/dt_minimum)) + 1)
 dtt = dt/ntt

 !Initialize variables to average the sub time steps
 storage_(:) = 0.0
 qout_(:) = 0.0
 storage1_(:) = storage(:)
 qout1_(:) = qout(:) 
 qin1_(:) = qin(:)

 !Set the number of threads
 !call omp_set_num_threads(nthreads)
 call mkl_set_num_threads(nthreads)
 !call openblas_set_num_threads(32)
          
 !Process the smaller time steps
 do itime = 1,ntt
  
   !Solve the kinematic wave for this time step
   call solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                             recharge1,area,dx,dtt,celerity,celerity1,&
                             wvalues,wcolumns,wrowindex)!,weights)

   !Add for the average at the end
   qout_ = qout_ + qout

 enddo

 !Set the output variables
 qout(:) = qout_/ntt
 qout1(:) = qout1_
 storage1(:) = storage1_
 qin1(:) = qin1_

end subroutine update

subroutine solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                                recharge1,area,dx,dtt,celerity,celerity1,&
                                wvalues,wcolumns,wrowindex)!,weights)

 implicit none
 integer,intent(in) :: nhsu,nvalues
 real*8,intent(in) :: dtt
 !real*4,intent(in),dimension(nhsu,nhsu) :: weights
 real*4,intent(inout),dimension(nhsu) :: storage,qout,qin,recharge,storage1,qout1,qin1,recharge1
 real*4,intent(in),dimension(nhsu) :: area,dx,celerity,celerity1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer :: i,j,iter,niter
 real*4 :: w,eps
 real*4,dimension(nhsu) :: denominator,numerator1,numerator2,qout_,ds,part1,part2

 !Define implicit scheme parameters
 w = 0.5 
 niter = 15!0
 eps = 10.0**-6.0

 !Calculate as many constants as possible before going into the iterations
 denominator = (1 + dtt*w*celerity/dx)
 numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
 numerator2 = dtt*w*celerity/dx
 part1 = numerator1/denominator
 part2 = numerator2/denominator

 !Update qout using the implicit scheme until it converges (upstream to downstream)
 do iter = 1,niter
 
  !Calculate qin
  call mkl_cspblas_scsrgemv('N',nhsu,wvalues,wrowindex,wcolumns,area*qout,qin)
  !call sgemv('N',nhsu,nhsu,1.0,weights,nhsu,qout*area,1,0.0,qin,1)
  !qout = qout*area
  !qin = 0
  !!$OMP PARALLEL DO SCHEDULE(STATIC)
  !do i = 1,nhsu
  !  do j = wrowindex(i)+1,wrowindex(i+1)
  !   qin(i) = qin(i) + qout(wcolumns(j)+1)*wvalues(j)
  !  enddo
  !enddo
  !!OMP PARALLEL END DO 
  !Divide by the area to get the flux
  qin = qin/area
 
  !Calculate qout
  qout = part1 + part2*qin

  !Set values below zero to 0
  where (qout .lt. 0.0) qout = 0.0

  !Determine if the tolerance has been achieved
  if (maxval(abs(qout_ - qout)) .le. eps) exit

  !Save values for comparison
  qout_ = qout

 enddo
 
 !Set all negative fluxes to 0
 where (qout < 0.0) qout = 0.0
 !Calculate the change in storage
 ds = dtt*((qin - qout)/dx + recharge)
 !ds = dtt*((qin - qout) + recharge)

 !Adjust the storages
 storage = storage + ds

 !Set the next time step's info
 qout1 = qout
 qin1 = qin

end subroutine




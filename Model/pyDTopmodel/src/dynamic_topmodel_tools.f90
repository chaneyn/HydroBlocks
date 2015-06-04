subroutine update(recharge,storage,qout,qin,recharge1,storage1,qout1,qin1,&
                  area,dx,dt,celerity,celerity1,wvalues,wcolumns,wrowindex,&
                  qin_outlet,area_outlet,nthreads,nhsu,nvalues,nhru_outlet)

 use omp_lib
 implicit none
 integer :: nhsu,nvalues,nhru_outlet
 real*8,intent(in) :: dt
 !real*4,intent(in),dimension(nhsu,nhsu) :: weights
 real*4,intent(in),dimension(nhsu) :: dx,area,celerity,celerity1
 real*4,intent(inout),dimension(nhsu) :: recharge,storage,recharge1,storage1,qout,qin,qout1,qin1
 real*4,intent(inout),dimension(nhru_outlet) :: qin_outlet
 real*4,intent(in),dimension(nhru_outlet) :: area_outlet
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer,intent(in) :: nthreads
 real*4,dimension(nhsu) :: storage_,qout_,storage1_,qout1_,qin1_,qin_outlet_
 real*4 :: dt_minimum,dtt
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
 qin_outlet_(:) = 0.0

 !Set the number of threads
 !call omp_set_num_threads(1)!nthreads)
 call mkl_set_num_threads(nthreads)
 !call openblas_set_num_threads(32)
          
 !Process the smaller time steps
 do itime = 1,ntt
  
   !Solve the kinematic wave for this time step
   call solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                             recharge1,area,dx,dtt,celerity,celerity1,&
                             wvalues,wcolumns,wrowindex,&
                             area_outlet,qin_outlet,nhru_outlet)!,weights)

   !Add for the average at the end
   qout_ = qout_ + qout
   qin_outlet_ = qin_outlet_ + qin_outlet

 enddo

 !Set the output variables (WHAT ABOUT qin???)
 qout(:) = qout_/ntt
 qout1(:) = qout1_
 storage1(:) = storage1_
 qin1(:) = qin1_
 qin_outlet(:) = qin_outlet_/ntt

end subroutine update

subroutine solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                                recharge1,area,dx,dtt,celerity,celerity1,&
                                wvalues,wcolumns,wrowindex,&
                                area_outlet,qin_outlet,nhru_outlet)!,weights)

 implicit none
 integer,intent(in) :: nhsu,nvalues,nhru_outlet
 real*4,intent(in) :: dtt
 !real*4,intent(in),dimension(nhsu,nhsu) :: weights
 real*4,intent(inout),dimension(nhsu) :: storage,qout,qin,recharge,storage1,qout1,qin1,recharge1
 real*4,intent(in),dimension(nhsu) :: area,dx,celerity,celerity1
 real*4,intent(inout),dimension(nhru_outlet) :: qin_outlet
 real*4,intent(in),dimension(nhru_outlet) :: area_outlet
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*4,intent(in),dimension(nvalues) :: wvalues
 integer :: i,j,iter,niter
 real*4 :: w,eps
 real*4,dimension(nhsu) :: denominator,numerator1,numerator2,qout_,ds,part1,part2
 real*4,dimension(nhsu+nhru_outlet) :: qin_all

 !Define implicit scheme parameters
 w = 0.5 
 niter = 15!0
 eps = 10.0**-6.0

 !Initialize dummy variables
 qin_all = 0.0

 !Calculate as many constants as possible before going into the iterations
 denominator = (1 + dtt*w*celerity/dx)
 numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
 numerator2 = dtt*w*celerity/dx
 part1 = numerator1/denominator
 part2 = numerator2/denominator

 !Update qout using the implicit scheme until it converges (upstream to downstream)
 do iter = 1,niter
 
  qin_all = 0.0
  !Calculate qin
  qout = area*qout
  call mkl_cspblas_scsrgemv('T',nhsu,wvalues,wrowindex,wcolumns,qout,qin_all)
  !print*,qin_all
  !Split up the qin among hrus and outlet hrus
 ! call scsrsgemv('N',nhsu,nhsu,1.0,weights,nhsu,qout*area,1,0.0,qin,1)
  !call cspblas_scsrgemv('T',nhsu,wvalues,wrowindex,wcolumns,qout,qin_all)
  !qout = qout*area
  !qin = 0
  !print*,wrowindex
  !print*,wcolumns
  !print*,size(wvalues),size(wrowindex)
  !!!$OMP PARALLEL DO SCHEDULE(STATIC)
  !do j = 1,size(wrowindex)-1
  !  do i = wrowindex(j)+1,wrowindex(j+1)
  !   !qin_all(i) = qin_all(i) + qout(wcolumns(j)+1)*wvalues(j)
  !   qin_all(wcolumns(i)+1) = qin_all(wcolumns(i)+1) + qout(i)*wvalues(j)
  !  enddo
  !enddo
  !qin_all(:) = 0.0001
  !!!OMP PARALLEL END DO 
  !Divide by the area to get the flux (hrus and outlets)
  qin = qin_all(1:nhsu)/area
  qin_outlet = qin_all(nhsu+1:nhsu+nhru_outlet)/area_outlet
 
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




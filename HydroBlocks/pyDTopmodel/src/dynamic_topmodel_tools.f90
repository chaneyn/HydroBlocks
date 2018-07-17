INCLUDE 'mkl_dss.f90' ! Include the standard DSS "header file."

module dss_variables

 use mkl_dss
 implicit none
 type(MKL_DSS_HANDLE) :: dss_handle ! Allocate storage for the solver handle.
 public :: dss_handle

end module

subroutine initialize_dss(wcolumns,wrowindex,nhsu,nvalues)

 use dss_variables
 implicit none
 integer*4,intent(in) :: nhsu,nvalues
 integer*4,intent(inout),dimension(nvalues) :: wcolumns
 integer*4,intent(inout),dimension(nhsu+1) :: wrowindex
 integer*4 :: opt,perm(nhsu)
 integer*8 :: error

 ! Initialize the solver 
 opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING
 error = dss_create(dss_handle,opt)

 ! Define the non-zero structure of the matrix.
 opt = MKL_DSS_NON_SYMMETRIC
 error = DSS_DEFINE_STRUCTURE(dss_handle,opt,wrowindex,nhsu,nhsu,wcolumns,nvalues)

 ! Reorder the matrix.
 opt = MKL_DSS_AUTO_ORDER
 error = DSS_REORDER(dss_handle,opt,perm )

end subroutine

subroutine finalize_dss()

 use dss_variables
 implicit none
 integer :: opt,error

 ! Delete the solver
 opt = MKL_DSS_DEFAULTS
 error = dss_delete(dss_handle,opt)

end subroutine

subroutine solve_kinematic_wave(nhsu,nvalues,storage,qout,qin,recharge,storage1,qout1,qin1,&
                                recharge1,scarea,dx,dtt,celerity,celerity1,&
                                wvalues,wcolumns,wrowindex,w,storage_mask)
                                !wvalues,wcolumns,wrowindex,w,dss_handle,storage_mask)
 !use mkl_dss
 use dss_variables
 implicit none
 integer*4,intent(in) :: nhsu,nvalues
 real*8,intent(in) :: dtt,w
 real*8,intent(inout),dimension(nhsu) :: storage,qout,qin,recharge,storage1,qout1,qin1,recharge1
 integer*8,intent(in),dimension(nhsu) :: storage_mask
 real*8,intent(in),dimension(nhsu) :: scarea,dx,celerity,celerity1
 integer*4,intent(in),dimension(nvalues) :: wcolumns
 integer*4,intent(in),dimension(nhsu+1) :: wrowindex
 real*8,intent(in),dimension(nvalues) :: wvalues
 !type(MKL_DSS_HANDLE),intent(inout) :: dss_handle ! Allocate storage for the solver handle.
 real*8,dimension(nvalues) :: A
 real*8,dimension(nhsu) :: denominator,numerator1,numerator2,part1,part2
 integer*8 :: i,j,error
 integer*4 :: opt,perm(nhsu)
 real*8 :: t0,t1

 !Calculate as many constants as possible before going into the iterations
 denominator = (1.0 + dtt*w*celerity/dx)
 numerator1 = qout1 + dtt*w*celerity*recharge + dtt*(1.0 - w)*celerity1*((qin1 - qout1)/dx + recharge1)
 numerator2 = dtt*w*celerity/dx
 part1 = numerator1/denominator
 part2 = numerator2/denominator

 !Prepare the A matrix (I - A) where A = ((F*(I*P2/SCAREA)).T*(I*SCAREA)).T
 do i = 1,nhsu
  A(wrowindex(i)+1:wrowindex(i+1)) = wvalues(wrowindex(i)+1:wrowindex(i+1))*&
                                     scarea(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)*&
                                     part2(i)/scarea(i)
  !A(wrowindex(i)+1:wrowindex(i+1)) = wvalues(wrowindex(i)+1:wrowindex(i+1))*&
  !                                   scarea(i)*&
  !                                   part2(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)/&
  !                                   scarea(wcolumns(wrowindex(i)+1:wrowindex(i+1))+1)
 enddo

 !I - A
 i = 0
 do j = 1,nvalues
  if ((wcolumns(j) .eq. i) .and. (j .ge. wrowindex(i+1) + 1))then
   A(j) = 1 - A(j)
   i = i + 1
  else
   A(j) = - A(j)
  endif
 enddo
   
 ! Factor the matrix.
 error = DSS_FACTOR_REAL(dss_handle,MKL_DSS_DEFAULTS,A)

 ! Solve the system of linear equations
 opt = MKL_DSS_DEFAULTS
 error = DSS_SOLVE_REAL(dss_handle,opt,part1,1,qout)

 !Set all negative fluxes to 0
 !where (qout < 0.0) qout = 0.0
 
 !Set outflows to 0 where the storage is below a threshold
 !where (storage_mask .eq. 0) qout = 0.0

 !Calculate qin
 !call mkl_cspblas_scsrgemv('N',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 call mkl_cspblas_dcsrgemv('N',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 !call mkl_cspblas_scsrgemv('T',nhsu,wvalues,wrowindex,wcolumns,scarea*qout,qin)
 qin = qin/scarea

 !Adjust the storages
 storage = storage + dtt*((qin - qout)/dx + recharge)

 !Set the next time step's info
 qout1 = qout
 qin1 = qin

end subroutine


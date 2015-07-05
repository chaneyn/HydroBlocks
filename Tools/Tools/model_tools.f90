module preprocessor

contains

subroutine calculate_connections(HSUs,dinfangle,transition_probabilities,nhsu,n,m)

 implicit none
 integer :: i,j,k,n,m
 integer,intent(in) :: nhsu
 real*8,intent(in),dimension(n,m) :: HSUs,dinfangle 
 real*8,intent(out),dimension(nhsu,nhsu) :: transition_probabilities
 integer :: is(8),js(8)
 real :: hsu_org,hsu_dst
 real :: angles(10),pi,angle,a,pct
 pi = 3.14159265359
 is = (/0,-1,-1,-1,0,1,1,1/)
 js = (/1,1,0,-1,-1,-1,0,1/)
 angles = (/-pi/4,0.,pi/4,pi/2,3*pi/4,pi,5*pi/4,6*pi/4,7*pi/4,2*pi/)

 !Iterate through the cells and determine the connections
 do i = 1,n
  do j = 1,m
   if (isnan(HSUs(i,j))) cycle
   angle = dinfangle(i,j)
   do k = 2,9
    pct = 0.0
    a = angle
    !Determine if the flow goes to that cell
    if (a > angles(k-1) .and. a < angles(k+1)) then
     if (k .eq. 1 .and. a .gt. pi) a = a - 2.0*pi
     if (a > angles(k)) then
      pct = (angles(k+1)-a)/(angles(k+1)-angles(k))
     else 
      pct = (a - angles(k-1))/(angles(k)-angles(k-1))
     endif
     if (pct < 1e-5) cycle
     !Add up the pct
     hsu_dst = HSUs(i+is(k-1),j+js(k-1)) + 1
     hsu_org = HSUs(i,j) + 1
     if (isnan(hsu_dst)) cycle
     transition_probabilities(int(hsu_org),int(hsu_dst)) = transition_probabilities(int(hsu_org),int(hsu_dst)) + pct
    endif
   enddo
  enddo
 enddo 

 !Compute the transition probabilities
 do i = 1,nhsu
  if (sum(transition_probabilities(i,:)) .eq. 0) cycle
  transition_probabilities(i,:) = transition_probabilities(i,:)/sum(transition_probabilities(i,:))
 end do

 end subroutine calculate_connections

 subroutine calculate_connections_d8(HSUs,d8dir,area,hrus_dst,hrus_org,nhsu,&
  max_nhsu,outlet_icoord,outlet_jcoord,outlet_hru,outlet_d8,n,m)

 implicit none
 integer :: i,j,k,n,m,ipos
 integer,intent(in) :: nhsu,max_nhsu
 real*8,intent(in),dimension(n,m) :: HSUs,d8dir,area
 !real*8,intent(out),dimension(nhsu,nhsu) :: transition_probabilities
 integer,intent(out),dimension(max_nhsu) :: hrus_dst,hrus_org
 integer,intent(out),dimension(max_nhsu) :: outlet_icoord,outlet_jcoord
 integer,intent(out),dimension(max_nhsu) :: outlet_hru,outlet_d8
 real*8,dimension(nhsu,nhsu) :: transition_probabilities
 real :: hsu_org,hsu_dst,minimum_area,maxarea
 hrus_dst(:) = -9999
 hrus_org(:) = -9999
 outlet_icoord(:) = -9999
 outlet_jcoord(:) = -9999
 outlet_hru(:) = -9999
 outlet_d8(:) = -9999
 ipos = 1
 minimum_area = 10**5
 maxarea = maxval(area)

 !Iterate through the cells and determine the connections
 do i = 1,n
  do j = 1,m
   if (isnan(HSUs(i,j))) cycle
     !Determine the HSU it flows into
     !0=N, 1=NE, 2=E, ... 7=NW"
     !if (d8dir(i,j) == 0) hsu_dst = HSUs(i+1,j) !N
     if (d8dir(i,j) == 0) hsu_dst = HSUs(i-1,j) !N
     !if (d8dir(i,j) == 1) hsu_dst = HSUs(i+1,j+1) !NE
     if (d8dir(i,j) == 1) hsu_dst = HSUs(i-1,j+1) !NE
     if (d8dir(i,j) == 2) hsu_dst = HSUs(i,j+1) !E
     !if (d8dir(i,j) == 3) hsu_dst = HSUs(i-1,j+1) !SE
     if (d8dir(i,j) == 3) hsu_dst = HSUs(i+1,j+1) !SE
     !if (d8dir(i,j) == 4) hsu_dst = HSUs(i-1,j) !S
     if (d8dir(i,j) == 4) hsu_dst = HSUs(i+1,j) !S
     !if (d8dir(i,j) == 5) hsu_dst = HSUs(i-1,j-1) !SW
     if (d8dir(i,j) == 5) hsu_dst = HSUs(i+1,j-1) !SW
     if (d8dir(i,j) == 6) hsu_dst = HSUs(i,j-1) !W
     !if (d8dir(i,j) == 7) hsu_dst = HSUs(i+1,j-1) !NW
     if (d8dir(i,j) == 7) hsu_dst = HSUs(i-1,j-1) !NW
     !If the grid cell is an outlet...
     if (isnan(hsu_dst) .and. (area(i,j) .ne. maxarea)) then
      !if (area(i,j) .lt. minimum_area) then
      hsu_dst = HSUs(i,j) + 1
      hsu_org = HSUs(i,j) + 1
      !hrus_dst(ipos) = hsu_dst
      hrus_org(ipos) = hsu_org
      !Send all this directly to the "outlet"
      hrus_dst(ipos) = 0
      !else
      ! print*,area(i,j)
      ! outlet_icoord(ipos) = i
      ! outlet_jcoord(ipos) = j
      ! outlet_hru(ipos) = HSUs(i,j) + 1
      ! outlet_d8(ipos) = d8dir(i,j)
      !endif
     else if (area(i,j) .eq. maxarea) then
      outlet_icoord(ipos) = i
      outlet_jcoord(ipos) = j
      outlet_hru(ipos) = HSUs(i,j) + 1
      outlet_d8(ipos) = d8dir(i,j)
     else
      hsu_dst = hsu_dst + 1
      !Determine the HSU it comes from
      hsu_org = HSUs(i,j) + 1
      !Add the count to the transition probabilities
      !transition_probabilities(int(hsu_org),int(hsu_dst)) = transition_probabilities(int(hsu_org),int(hsu_dst)) + 1
      hrus_dst(ipos) = hsu_dst
      hrus_org(ipos) = hsu_org
     endif
     ipos = ipos + 1
  enddo
 enddo
 !Compute the transition probabilities
 !do i = 1,nhsu
 ! if (sum(transition_probabilities(i,:)) .eq. 0) cycle
 ! transition_probabilities(i,:) = transition_probabilities(i,:)/sum(transition_probabilities(i,:))
 !end do

 end subroutine

 subroutine calculate_hsus()

  !Read in 3d array of covariates
  !Read in 2d array of threshholds

 end subroutine calculate_hsus

end module preprocessor




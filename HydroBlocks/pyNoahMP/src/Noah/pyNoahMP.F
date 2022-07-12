!module model
!
! implicit none
!
! !Declare variables
!
! !Miscellanous parameters
! real :: badval = 0.0!-1.E36 
! real :: co2 = 395.e-06
! real :: o2 = 0.209
!  
! !General parameters
! integer :: nsoil ! Number of soil layers
! integer :: ncells ! Number of grid cells
! integer :: nsnow = 3 ! Number of snow layers
! real :: dt ! Timestep in seconds
! real :: dx ! Spatial resolution (m)
! real :: julian ! Julian day
! integer :: yearlen
! logical :: iz0tlnd
! integer :: itime
! character(len=256) :: llanduse
! character(len=256) :: lsoil
! !character(:), allocatable :: llanduse,lsoil
! character(len=256) :: vegparm_file,soilparm_file,genparm_file,mptable_file
! character(len=19) :: nowdate
!
! !Model options
! integer :: idveg ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1      
! integer :: iopt_crs ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
! integer :: iopt_btr ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
! integer :: iopt_run ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
! integer :: iopt_sfc ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
! integer :: iopt_frz ! supercooled liquid water (1-> NY06; 2->Koren99)
! integer :: iopt_inf ! frozen soil permeability (1-> NY06; 2->Koren99)
! integer :: iopt_rad ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
! integer :: iopt_alb ! snow surface albedo (1->BATS; 2->CLASS)
! integer :: iopt_snf ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
! integer :: iopt_tbot ! lower boundary of soil temperature (1->zero-flux; 2->Noah) 
! integer :: iopt_stc ! snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)
!
! !Single layer
! integer,dimension(:),allocatable :: isnow ! actual number of snow layers
! real,dimension(:),allocatable :: z_ml ! model height [m]
! real,dimension(:),allocatable :: lwdn       ! Downward longwave radiation flux at surface [W/m2]
! real,dimension(:),allocatable :: swdn      ! Downward shortwave radiation flux at surface [W/m2]
! real,dimension(:),allocatable :: p_ml    !  Surface pressure defined at intermediate level [Pa]
! real,dimension(:),allocatable :: psfc    ! Surface pressure [Pa]
! real,dimension(:),allocatable :: prcp       ! Precipitation rate (kg m-2 s-1)
! real,dimension(:),allocatable :: t_ml     ! Air temperature (K) [Forcing]
! real,dimension(:),allocatable :: q_ml         ! Surface specific humidity (kg kg-1)
! real,dimension(:),allocatable :: u_ml       ! West-to-east component of the surface [m/s]
! real,dimension(:),allocatable :: v_ml       ! North-to-south component of the surface [m/s]
! real,dimension(:),allocatable :: fsh ! total sensible heat (w/m2)
! real,dimension(:),allocatable :: ssoil ! soil heat (w/m2)
! real,dimension(:),allocatable :: salb ! surface albedo
! real,dimension(:),allocatable :: fsno ! snow cover fraction
! real,dimension(:),allocatable :: swe ! snow water equivalent (mm)
! real,dimension(:),allocatable :: sndpth ! snow depth (m)
! real,dimension(:),allocatable :: emissi ! net surface emissivity
! real,dimension(:),allocatable :: qsfc1d ! bulk surface specific humidity (kg/kg)
! real,dimension(:),allocatable :: tv ! vegetation canopy temperature
! real,dimension(:),allocatable :: tg ! ground surface temperature
! real,dimension(:),allocatable :: canice ! canopy-intercepted ice (mm)
! real,dimension(:),allocatable :: canliq ! canopy-intercepted liquid water (mm)
! real,dimension(:),allocatable :: eah ! canopy air vapor pressure (pa)
! real,dimension(:),allocatable :: tah ! canopy air temperature (K)
! real,dimension(:),allocatable :: cm ! momentum drag coefficient
! real,dimension(:),allocatable :: ch ! sensible heat exchange coefficient
! real,dimension(:),allocatable :: fwet ! weted or snowed fraction of the canopy
! real,dimension(:),allocatable :: sneqvo ! snow mass at last time step (mm h2o)
! real,dimension(:),allocatable :: albold ! snow albedo at last time step
! real,dimension(:),allocatable :: qsnow ! snowfall on the ground (mm/s)
! real,dimension(:),allocatable :: wslake ! lake water storage (mm)
! real,dimension(:),allocatable :: zwt ! water table depth (m)
! real,dimension(:),allocatable :: dzwt ! change in water table depth (m)
! real,dimension(:),allocatable :: wa ! water in the aquifer (mm)
! real,dimension(:),allocatable :: wt ! groundwater storage (mm)
! real,dimension(:),allocatable :: smcwtd ! soil moisture content in the transmission zone
! real,dimension(:),allocatable :: deeprech ! recharge to the water table when deep
! real,dimension(:),allocatable :: rech ! recharge to the water table (diagnostic)
! real,dimension(:),allocatable :: lfmass ! leaf mass (g/m2)
! real,dimension(:),allocatable :: rtmass ! mass of fine roots (g/m2)
! real,dimension(:),allocatable :: stmass ! stem mass (g/m2)
! real,dimension(:),allocatable :: wood ! mass of wood (incl. woody roots) (g/m2)
! real,dimension(:),allocatable :: stblcp ! stable carbon in deep soil (g/m2)
! real,dimension(:),allocatable :: fastcp ! short lived carbon, shallow soil (g/m2)
! real,dimension(:),allocatable :: plai ! leaf area index
! real,dimension(:),allocatable :: psai ! stem area index
! real,dimension(:),allocatable :: tauss ! non-dimensional snow age
! real,dimension(:),allocatable :: t2mv ! 2m temperature of vegetation part
! real,dimension(:),allocatable :: t2mb ! 2m temperature of bare ground part
! real,dimension(:),allocatable :: q2mv ! 2m mixing ratio of vegetation part
! real,dimension(:),allocatable :: q2mb ! 2m mixing ratio of bare ground part
! real,dimension(:),allocatable :: trad ! surface radiative temperature (k)
! real,dimension(:),allocatable :: nee ! net ecosys exchange (g/m2/s CO2)
! real,dimension(:),allocatable :: gpp ! gross primary assimilation [g/m2/s C]
! real,dimension(:),allocatable :: npp ! net primary productivity [g/m2/s C]
! real,dimension(:),allocatable :: fvegmp ! greenness vegetation fraction [-]
! real,dimension(:),allocatable :: runsf ! surface runoff [mm/s]
! real,dimension(:),allocatable :: runsb ! subsurface runoff [mm/s]
! real,dimension(:),allocatable :: ecan ! evaporation of intercepted water (mm/s)
! real,dimension(:),allocatable :: etran ! transpiration rate (mm/s)
! real,dimension(:),allocatable :: esoil ! soil surface evaporation rate (mm/s]
! real,dimension(:),allocatable :: fsa ! total absorbed solar radiation (w/m2)
! real,dimension(:),allocatable :: fira ! total net longwave rad (w/m2) [+ to atm]
! real,dimension(:),allocatable :: apar ! photosyn active energy by canopy (w/m2)
! real,dimension(:),allocatable :: psn ! total photosynthesis (umol co2/m2/s) [+]
! real,dimension(:),allocatable :: sav ! solar rad absorbed by veg. (w/m2)
! real,dimension(:),allocatable :: sag ! solar rad absorbed by ground (w/m2)
! real,dimension(:),allocatable :: rssun ! sunlit leaf stomatal resistance (s/m)
! real,dimension(:),allocatable :: rssha ! shaded leaf stomatal resistance (s/m)
! real,dimension(:),allocatable :: bgap ! between gap fraction
! real,dimension(:),allocatable :: wgap ! within gap fraction
! real,dimension(:),allocatable :: tgv ! under canopy ground temperature [K]
! real,dimension(:),allocatable :: tgb ! bare ground temperature [K]
! real,dimension(:),allocatable :: chv ! sensible heat exchange coefficient vegetated
! real,dimension(:),allocatable :: chb ! sensible heat exchange coefficient bare-ground
! real,dimension(:),allocatable :: irc ! canopy net LW rad. [w/m2] [+ to atm]
! real,dimension(:),allocatable :: irg ! veg ground net LW rad. [w/m2] [+ to atm]
! real,dimension(:),allocatable :: shc ! canopy sen. heat [w/m2]   [+ to atm]
! real,dimension(:),allocatable :: shg ! veg ground sen. heat [w/m2]   [+ to atm]
! real,dimension(:),allocatable :: evg ! veg ground evap. heat [w/m2]  [+ to atm]
! real,dimension(:),allocatable :: ghv ! veg ground heat flux [w/m2]  [+ to soil]
! real,dimension(:),allocatable :: irb ! bare net longwave rad. [w/m2] [+ to atm]
! real,dimension(:),allocatable :: shb ! bare sensible heat [w/m2]     [+ to atm]
! real,dimension(:),allocatable :: evb ! bare evaporation heat [w/m2]  [+ to atm]
! real,dimension(:),allocatable :: ghb ! bare ground heat flux [w/m2] [+ to soil]
! real,dimension(:),allocatable :: tr ! transpiration [w/m2]  [+ to atm]
! real,dimension(:),allocatable :: evc ! canopy evaporation heat [w/m2]  [+ to atm]
! real,dimension(:),allocatable :: chleaf ! leaf exchange coefficient
! real,dimension(:),allocatable :: chuc ! under canopy exchange coefficient
! real,dimension(:),allocatable :: chv2 ! veg 2m exchange coefficient
! real,dimension(:),allocatable :: chb2 ! bare 2m exchange coefficient 
! real,dimension(:),allocatable :: cosz ! cosine zenith angle
! real,dimension(:),allocatable :: lat  ! latitude [rad]
! real,dimension(:),allocatable :: lon  ! longitude [rad]
! real,dimension(:),allocatable :: fveg ! vegetation fraction
! real,dimension(:),allocatable :: fvgmax ! annual max vegetation
! real,dimension(:),allocatable :: fpice ! snow fraction of precip
! real,dimension(:),allocatable :: fcev ! canopy evaporation heat (w/m2) [+ to atm]
! real,dimension(:),allocatable :: fgev ! ground evaporation heat (w/m2) [+ to atm]
! real,dimension(:),allocatable :: fctr ! transpiration heat flux (w/m2) [+ to atm]
! real,dimension(:),allocatable :: qsnbot ! snowmelt out bottom of pack [mm/s]
! real,dimension(:),allocatable :: ponding ! snowmelt with no pack [mm]
! real,dimension(:),allocatable :: ponding1 ! snowmelt with no pack [mm]
! real,dimension(:),allocatable :: ponding2 ! snowmelt with no pack [mm]
! real,dimension(:),allocatable :: fsr ! total reflected solar radiation (w/m2)
! real,dimension(:),allocatable :: co2pp ! CO2 partial pressure [Pa]
! real,dimension(:),allocatable :: o2pp ! O2 partial pressure [Pa]
! real,dimension(:),allocatable :: foln ! nitrogen saturation [%
! real,dimension(:),allocatable :: tbot ! deep soil temperature [K]
! integer,dimension(:),allocatable :: isurban ! if cell is urban
! integer,dimension(:),allocatable :: slopetyp
! integer,dimension(:),allocatable :: soiltyp ! soil type
! integer,dimension(:),allocatable :: vegtyp ! vegetation type
! integer,dimension(:),allocatable :: ice ! glacier cell
! integer,dimension(:),allocatable :: isc ! soil color index
! integer,dimension(:),allocatable :: ist ! surface type 1-soil; 1-lake
! real,dimension(:),allocatable :: smcmax ! saturated soil moisture
! real,dimension(:),allocatable :: smcdry ! residual soil moisture
! real,dimension(:),allocatable :: smcref ! field capacity soil moisture
! real,dimension(:),allocatable :: errwat ! water balance error
! real,dimension(:),allocatable :: si0 ! initial soil moisture deficit
! real,dimension(:),allocatable :: si1 ! final soil moisture deficit
! real,dimension(:),allocatable :: zwt0 ! initial water table depth
! real,dimension(:),allocatable :: minzwt ! minimum water table depth
! real,dimension(:),allocatable :: co2air
! real,dimension(:),allocatable :: o2air
! real,dimension(:),allocatable :: bb0
! real,dimension(:),allocatable :: drysmc0
! real,dimension(:),allocatable :: f110
! real,dimension(:),allocatable :: maxsmc0
! real,dimension(:),allocatable :: refsmc0
! real,dimension(:),allocatable :: satpsi0
! real,dimension(:),allocatable :: satdk0
! real,dimension(:),allocatable :: satdw0
! real,dimension(:),allocatable :: wltsmc0
! real,dimension(:),allocatable :: qtz0
!
! ! Noemi
! integer,dimension(:),allocatable :: root_depth
!
! !Multi layer
! real,dimension(:,:),allocatable :: stc ! snow/soil temperatures
! real,dimension(:,:),allocatable :: sh2o ! vol. soil liquid water (m3/m3)
! real,dimension(:,:),allocatable :: smc ! vol. soil moisture (m3/m3)
! real,dimension(:,:),allocatable :: smceq ! eq vol. soil moisture (m3/m3)
! real,dimension(:,:),allocatable :: zsnso ! snow layer depth (m)
! real,dimension(:,:),allocatable :: snice ! snow layer ice (mm)
! real,dimension(:,:),allocatable :: snliq ! snow layer liquid water (mm)
! real,dimension(:,:),allocatable :: ficeold ! snow layer ice fraction 
! real,dimension(:,:),allocatable :: zsoil ! depth to soil interfaces
! real,dimension(:,:),allocatable :: sldpth ! soil layer thickness
! real,dimension(:,:),allocatable :: hdiv ! horizontal divergence (mm/s) [+]
!
 subroutine initialize_parameters(llanduse,lsoil,vegparm_file,genparm_file,&
         IOPT_CRS,IOPT_BTR,IOPT_RUN,IOPT_SFC,IOPT_FRZ,IOPT_INF,IOPT_RAD,IOPT_ALB,&
         IOPT_SNF,IOPT_TBOT,IOPT_STC,IDVEG,&
         mptable_file,BB0,DRYSMC0,F110,MAXSMC0,REFSMC0,SATPSI0,SATDK0,SATDW0,WLTSMC0,QTZ0,ncells)

  USE module_sf_noahmplsm,only : read_mp_veg_parameters,noahmp_options
  implicit none
  integer :: i
  integer,intent(in) :: ncells,IOPT_CRS,IOPT_BTR,IOPT_RUN,IOPT_SFC,IOPT_FRZ,IOPT_INF,IOPT_RAD,IOPT_ALB
  integer,intent(in) :: IOPT_SNF,IOPT_TBOT,IOPT_STC,IDVEG
  character(len=256),intent(in) :: llanduse,lsoil,vegparm_file,genparm_file,mptable_file
  real,dimension(ncells),intent(in) :: BB0,DRYSMC0,F110,MAXSMC0,REFSMC0,SATPSI0,SATDK0,SATDW0,WLTSMC0,QTZ0
  character(len=256) :: soilparm_file

  ! Read our lookup tables and parameter tables:  VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
  call soil_veg_gen_parm(llanduse,lsoil,vegparm_file,soilparm_file,genparm_file,&
         BB0,DRYSMC0,F110,MAXSMC0,REFSMC0,SATPSI0,SATDK0,SATDW0,WLTSMC0,QTZ0,ncells)

  ! Read the Noah-MP table
  call read_mp_veg_parameters(llanduse,mptable_file)

  ! Define the NOAH-MP options
  CALL NOAHMP_OPTIONS(IDVEG  ,IOPT_CRS  ,IOPT_BTR  ,IOPT_RUN,IOPT_SFC,IOPT_FRZ , &
                     IOPT_INF  ,IOPT_RAD  ,IOPT_ALB  ,IOPT_SNF,IOPT_TBOT,IOPT_STC )

 end subroutine

 subroutine run_model(ncores,yearlen,idveg,ioptcrs,iopt_btr,iopt_run,iopt_sfc,iopt_frz,iopt_inf,iopt_rad,iopt_alb,& !0D,integer
                      iopt_snf,iopt_tbot,iopt_stc,itime,iz0tlnd,& !0D,integer
                      dt,dx,julian,& !0D,real
                      isnow,isurban,slopetyp,soiltyp,vegtyp,ice,isc,ist,root_depth,& !1D,integer
                      nowdate,& !1D,character
                      z_ml,lwdn,swdn,p_ml,psfc,prcp,t_ml,q_ml,u_ml,v_ml,fsh,ssoil,salb,fsno,swe,sndpth,emissi,& !1D,real
                      qsfc1d,tv,tg,canice,canliq,eah,tah,cm,ch,fwet,sneqvo,albold,qsnow,wslake,zwt,dzwt,wa,wt,smcwtd,& !1D,real
                      deeprech,rech,lfmass,rtmass,stmass,wood,stblcp,fastcp,plai,psai,tauss,t2mv,t2mb,q2mv,q2mb,& !1D,real
                      trad,nee,gpp,npp,fvegmp,runsf,runsb,ecan,etran,esoil,fsa,fira,apar,psn,sav,sag,rssun,rssha,& !1D,real
                      bgap,wgap,tgv,tgb,chv,chb,irc,irg,shc,shg,evg,ghv,irb,shb,evb,ghb,tr,evc,chleaf,chuc,chv2,& !1D,real
                      chb2,cosz,lat,lon,fveg,fvgmax,fpice,fcev,fgev,fctr,qsnbot,ponding,ponding1,ponding2,fsr,& !1D,real
                      co2pp,o2pp,foln,tbot,smcmax,smcdry,smcref,errwat,si0,si1,zwt0,minzwt,co2air,o2air,bb0,& !1D,real
                      drysmc0,f110,maxsmc0,refsmc0,satpsi0,satdk0,satdw0,wltsmc0,qtz0,& !1D,real
                      stc,sh2o,smc,smceq,zsnso,snice,snliq,ficeold,zsoil,sldpth,hdiv,&!2D,real
                      ncells,nsoil,nsnow)
                      

  implicit none
  integer :: i,ncores,omp_get_thread_num
  real*8 :: ostart,oend,ostart1,oend1,omp_get_wtime,tmp
  integer,intent(inout) :: ncells,nsoil,nsnow,yearlen,idveg,ioptcrs,iopt_btr,iopt_run,iopt_sfc,iopt_frz,iopt_inf,iopt_rad,iopt_alb,&
                        iopt_snf,iopt_tbot,iopt_stc,itime,iz0tlnd
  real,intent(inout) :: dt,dx,julian
  character(len=19),intent(inout) :: nowdate
  integer,dimension(ncells),intent(inout) :: isnow,isurban,slopetyp,soiltyp,vegtyp,ice,isc,ist,root_depth
  real,dimension(ncells),intent(inout) :: z_ml,lwdn,swdn,p_ml,psfc,prcp,t_ml,q_ml,u_ml,v_ml,fsh,ssoil,salb,fsno,swe,sndpth,emissi,&
                                          qsfc1d,tv,tg,canice,canliq,eah,tah,cm,ch,fwet,sneqvo,albold,qsnow,wslake,zwt,dzwt,wa,wt,smcwtd,&
                                          deeprech,rech,lfmass,rtmass,stmass,wood,stblcp,fastcp,plai,psai,tauss,t2mv,t2mb,q2mv,q2mb,&
                                          trad,nee,gpp,npp,fvegmp,runsf,runsb,ecan,etran,esoil,fsa,fira,apar,psn,sav,sag,rssun,rssha,&
                                          bgap,wgap,tgv,tgb,chv,chb,irc,irg,shc,shg,evg,ghv,irb,shb,evb,ghb,tr,evc,chleaf,chuc,chv2,&
                                          chb2,cosz,lat,lon,fveg,fvgmax,fpice,fcev,fgev,fctr,qsnbot,ponding,ponding1,ponding2,fsr,&
                                          co2pp,o2pp,foln,tbot,smcmax,smcdry,smcref,errwat,si0,si1,zwt0,minzwt,co2air,o2air,bb0,&
                                          drysmc0,f110,maxsmc0,refsmc0,satpsi0,satdk0,satdw0,wltsmc0,qtz0 
  real,intent(inout),dimension(ncells,-nsnow+1:nsoil) :: zsnso,stc
  real,intent(inout),dimension(ncells,-nsnow+1:0) :: snice,snliq,ficeold
  real,dimension(ncells,nsoil),intent(inout) :: sh2o,smc,smceq,zsoil,sldpth,hdiv
  call OMP_SET_NUM_THREADS(ncores)
  !$OMP PARALLEL
  !write(*,*)'Thread rank: ', OMP_GET_THREAD_NUM()
  !$OMP DO

  do i = 1,ncells

    call run_model_cell(&
            LAT(i)     , YEARLEN , JULIAN  , COSZ(i)    , & ! IN : Time/Space-related
            DT      , DX      , NSOIL   , ZSOIL(i,:)   , NSNOW,  & ! IN : Model configuration 
            FVEG(i)    , FVGMAX(i)  , VEGTYP(i)  , ISURBAN(i) , ICE(i)     , IST(i)     , & ! IN : Vegetation/Soil characteristics
            ISC(i)     , SMCEQ(i,:)   ,                                         & ! IN : Vegetation/Soil characteristics
            IZ0TLND ,                                                   & ! IN : User options
            T_ML(i)    , P_ML(i)    , PSFC(i)    , U_ML(i)    , V_ML(i)    , Q_ML(i)    , & ! IN : Forcing
            SWDN(i)    , LWDN(i)    , PRCP(i)    , TBOT(i)    , CO2PP(i)   , & ! IN : Forcing
            O2PP(i)    , FOLN(i)    , FICEOLD(i,:) , Z_ML(i)    ,           & ! IN : Forcing
            HDIV(i,:), & ! IN :
            ALBOLD(i)  , SNEQVO(i)  ,                                         & ! IN/OUT : 
            STC(i,:)     , SH2O(i,:)   , SMC(i,:)     , TAH(i)     , EAH(i)     , FWET(i)    , & ! IN/OUT : 
            CANLIQ(i)  , CANICE(i)  , TV(i)      , TG(i)      , QSFC1D(i)  , QSNOW(i)   , & ! IN/OUT : 
            ISNOW(i)   , ZSNSO(i,:)   , SNDPTH(i)  , SWE(i)     , SNICE(i,:)   , SNLIQ(i,:)   , & ! IN/OUT : 
            ZWT(i)     , WA(i)      , WT(i)      , WSLAKE(i)  , LFMASS(i)  , RTMASS(i)  , & ! IN/OUT : 
            STMASS(i)  , WOOD(i)    , STBLCP(i)  , FASTCP(i)  , PLAI(i)    , PSAI(i)    , & ! IN/OUT : 
            CM(i)      , CH(i)      , TAUSS(i)   ,                               & ! IN/OUT : 
            SMCWTD(i)  ,DEEPRECH(i) , RECH(i)    ,                               & ! IN/OUT :
            FSA(i)     , FSR(i)     , FIRA(i)    , FSH(i)     , SSOIL(i)   , FCEV(i)    , & ! OUT : 
            FGEV(i)    , FCTR(i)    , ECAN(i)    , ETRAN(i)   , ESOIL(i)   , TRAD(i)    , & ! OUT : 
            TGB(i)     , TGV(i)     , T2MV(i)    , T2MB(i)    , Q2MV(i)    , Q2MB(i)    , & ! OUT : 
            RUNSF(i)   , RUNSB(i)   , APAR(i)    , PSN(i)    , SAV(i)     , SAG(i)     , & ! OUT : 
            FSNO(i)    , NEE(i)     , GPP(i)     , NPP(i)     , FVEGMP(i)  , SALB(i)    , & ! OUT : 
            QSNBOT(i)  , PONDING(i) , PONDING1(i), PONDING2(i), RSSUN(i)   , RSSHA(i)   , & ! OUT : 
            BGAP(i)    , WGAP(i)    , CHV(i)     , CHB(i)     , EMISSI(i)  ,           & ! OUT : 
            SHG(i)     , SHC(i)     , SHB(i)     , EVG(i)     , EVB(i)     , GHV(i)     , & ! OUT :
            GHB(i)     , IRG(i)     , IRC(i)     , IRB(i)     , TR(i)      , EVC(i)     , & ! OUT :
            CHLEAF(i)  , CHUC(i)    , CHV2(i)    , CHB2(i)    , FPICE(i), &
            SOILTYP(i), SLOPETYP(i),LON(i),NOWDATE,ITIME,SLDPTH(i,:),DZWT(i),ERRWAT(i),si0(i),si1(i), &
            zwt0(i),minzwt(i),ROOT_DEPTH(i)) ! OTHER
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

 end subroutine

!end module model

subroutine run_model_cell(&
            LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
            DT      , DX      , NSOIL   , ZSOIL   , NSNOW,  & ! IN : Model configuration 
            FVEG    , FVGMAX  , VEGTYP  , ISURBAN , ICE     , IST     , & ! IN : Vegetation/Soil characteristics
            ISC     , SMCEQ   ,                                         & ! IN : Vegetation/Soil characteristics
            IZ0TLND ,                                                   & ! IN : User options
            T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
            SWDN    , LWDN    , PRCP    , TBOT    , CO2PP   , & ! IN : Forcing
            O2PP    , FOLN    , FICEOLD , Z_ML    ,           & ! IN : Forcing
            HDIV    , & ! IN : Horizontal
            ALBOLD  , SNEQVO  ,                                         & ! IN/OUT : 
            STC     , SMH2O   , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT : 
            CANLIQ  , CANICE  , TV      , TG      , QSFC1D  , QSNOW   , & ! IN/OUT : 
            ISNOW   , ZSNSO   , SNDPTH  , SWE     , SNICE   , SNLIQ   , & ! IN/OUT : 
            ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT : 
            STMASS  , WOOD    , STBLCP  , FASTCP  , PLAI    , PSAI    , & ! IN/OUT : 
            CM      , CH      , TAUSS   ,                               & ! IN/OUT : 
            SMCWTD  ,DEEPRECH , RECH    ,                               & ! IN/OUT :
            FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT : 
            FGEV    , FCTR    , ECAN    , ETRAN   , ESOIL   , TRAD    , & ! OUT : 
            TGB     , TGV     , T2MV    , T2MB    , Q2MV    , Q2MB    , & ! OUT : 
            RUNSF   , RUNSB   , APAR    , PSN     , SAV     , SAG     , & ! OUT : 
            FSNO    , NEE     , GPP     , NPP     , FVEGMP  , SALB    , & ! OUT : 
            QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT : 
            BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT : 
            SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
	    GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
	    CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE, &
            SOILTYP, SLOPETYP,LON,NOWDATE,ITIME,SLDPTH,DZWT,ERRWAT,si0,si1,zwt0,minzwt,ROOT_DEPTH) ! OTHER

  USE module_sf_noahmplsm
  USE module_sf_noahmp_glacier

  implicit none

    !OTHER
    real*8 :: omp_get_wtime,ostart,oend
    REAL :: JULIAN,DT,DX,LON
    CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE
    INTEGER :: NSOIL,ISURBAN,IZ0TLND,NSNOW
    
  
    !SFLX
    REAL                                :: COSZ         ! cosine zenith angle
    REAL                                :: LAT          ! latitude [rad]
    REAL                                :: Z_ML         ! model height [m]
    INTEGER                             :: VEGTYP       ! vegetation type
    INTEGER                             :: SOILTYP      ! soil type
    REAL                                :: FVEG         ! vegetation fraction [-]
    REAL                                :: FVGMAX       ! annual max vegetation fraction []
    REAL                                :: TBOT         ! deep soil temperature [K]
    REAL                                :: T_ML         ! temperature valid at mid-levels [K]
    REAL                                :: Q_ML         ! water vapor mixing ratio [kg/kg_dry]
    REAL                                :: U_ML         ! U wind component [m/s]
    REAL                                :: V_ML         ! V wind component [m/s]
    REAL                                :: SWDN         ! solar down at surface [W m-2]
    REAL                                :: LWDN         ! longwave down at surface [W m-2]
    REAL                                :: P_ML         ! pressure, valid at interface [Pa]
    REAL                                :: PSFC         ! surface pressure [Pa]
    REAL                                :: PRCP         ! precipitation entering land model [mm]
    REAL,   DIMENSION( 1:NSOIL)         :: HDIV         ! 

! INOUT (with generic LSM equivalent)

    REAL                                :: FSH          ! total sensible heat (w/m2) [+ to atm]
    REAL                                :: SSOIL        ! soil heat heat (w/m2) 
    REAL                                :: SALB         ! surface albedo (-)
    REAL                                :: FSNO         ! snow cover fraction (-)
    REAL,   DIMENSION( 1:NSOIL)         :: SMCEQ        ! eq vol. soil moisture (m3/m3)
    REAL,   DIMENSION( 1:NSOIL)         :: SMC          ! vol. soil moisture (m3/m3)
    REAL,   DIMENSION( 1:NSOIL)         :: SMH2O        ! vol. soil liquid water (m3/m3)
    REAL,   DIMENSION(-2:NSOIL)         :: STC          ! snow/soil tmperatures
    REAL                                :: SWE          ! snow water equivalent (mm)
    REAL                                :: SNDPTH       ! snow depth (m)
    REAL                                :: EMISSI       ! net surface emissivity
    REAL                                :: QSFC1D       ! bulk surface specific humidity

! INOUT (with no Noah LSM equivalent)

    INTEGER                             :: ISNOW        ! actual no. of snow layers
    REAL                                :: TV           ! vegetation canopy temperature
    REAL                                :: TG           ! ground surface temperature
    REAL                                :: CANICE       ! canopy-intercepted ice (mm)
    REAL                                :: CANLIQ       ! canopy-intercepted liquid water (mm)
    REAL                                :: EAH          ! canopy air vapor pressure (pa)
    REAL                                :: TAH          ! canopy air temperature (k)
    REAL                                :: CM           ! momentum drag coefficient
    REAL                                :: CH           ! sensible heat exchange coefficient
    REAL                                :: FWET         ! wetted or snowed fraction of the canopy (-)
    REAL                                :: SNEQVO       ! snow mass at last time step(mm h2o)
    REAL                                :: ALBOLD       ! snow albedo at last time step (-)
    REAL                                :: QSNOW        ! snowfall on the ground [mm/s]
    REAL                                :: WSLAKE       ! lake water storage [mm]
    REAL                                :: ZWT          ! water table depth [m]
    REAL                                :: WA           ! water in the "aquifer" [mm]
    REAL                                :: WT           ! groundwater storage [mm]
    REAL                                :: SMCWTD       ! soil moisture content in the layer to the water table when deep
    REAL                                :: DEEPRECH     ! recharge to the water table when deep
    REAL                                :: RECH         ! recharge to the water table (diagnostic)  
    REAL, DIMENSION(-2:NSOIL)           :: ZSNSO        ! snow layer depth [m]
    REAL, DIMENSION(-2:              0) :: SNICE        ! snow layer ice [mm]
    REAL, DIMENSION(-2:              0) :: SNLIQ        ! snow layer liquid water [mm]
    REAL                                :: LFMASS       ! leaf mass [g/m2]
    REAL                                :: RTMASS       ! mass of fine roots [g/m2]
    REAL                                :: STMASS       ! stem mass [g/m2]
    REAL                                :: WOOD         ! mass of wood (incl. woody roots) [g/m2]
    REAL                                :: STBLCP       ! stable carbon in deep soil [g/m2]
    REAL                                :: FASTCP       ! short-lived carbon, shallow soil [g/m2]
    REAL                                :: PLAI         ! leaf area index
    REAL                                :: PSAI         ! stem area index
    REAL                                :: TAUSS        ! non-dimensional snow age

! OUT (with no Noah LSM equivalent)

    REAL                                :: T2MV         ! 2m temperature of vegetation part
    REAL                                :: T2MB         ! 2m temperature of bare ground part
    REAL                                :: Q2MV         ! 2m mixing ratio of vegetation part
    REAL                                :: Q2MB         ! 2m mixing ratio of bare ground part
    REAL                                :: TRAD         ! surface radiative temperature (k)
    REAL                                :: NEE          ! net ecosys exchange (g/m2/s CO2)
    REAL                                :: GPP          ! gross primary assimilation [g/m2/s C]
    REAL                                :: NPP          ! net primary productivity [g/m2/s C]
    REAL                                :: FVEGMP       ! greenness vegetation fraction [-]
    REAL                                :: RUNSF        ! surface runoff [mm/s]
    REAL                                :: RUNSB        ! subsurface runoff [mm/s]
    REAL                                :: ECAN         ! evaporation of intercepted water (mm/s)
    REAL                                :: ETRAN        ! transpiration rate (mm/s)
    REAL                                :: ESOIL        ! soil surface evaporation rate (mm/s]
    REAL                                :: FSA          ! total absorbed solar radiation (w/m2)
    REAL                                :: FIRA         ! total net longwave rad (w/m2) [+ to atm]
    REAL                                :: APAR         ! photosyn active energy by canopy (w/m2)
    REAL                                :: PSN          ! total photosynthesis (umol co2/m2/s) [+]
    REAL                                :: SAV          ! solar rad absorbed by veg. (w/m2)
    REAL                                :: SAG          ! solar rad absorbed by ground (w/m2)
    REAL                                :: RSSUN        ! sunlit leaf stomatal resistance (s/m)
    REAL                                :: RSSHA        ! shaded leaf stomatal resistance (s/m)
    REAL                                :: BGAP         ! between gap fraction
    REAL                                :: WGAP         ! within gap fraction
    REAL                                :: TGV          ! under canopy ground temperature[K]
    REAL                                :: TGB          ! bare ground temperature [K]
    REAL                                :: CHV          ! sensible heat exchange coefficient vegetated
    REAL                                :: CHB          ! sensible heat exchange coefficient bare-ground
    REAL                                :: IRC          ! canopy net LW rad. [w/m2] [+ to atm]
    REAL                                :: IRG          ! veg ground net LW rad. [w/m2] [+ to atm]
    REAL                                :: SHC          ! canopy sen. heat [w/m2]   [+ to atm]
    REAL                                :: SHG          ! veg ground sen. heat [w/m2]   [+ to atm]
    REAL                                :: EVG          ! veg ground evap. heat [w/m2]  [+ to atm]
    REAL                                :: GHV          ! veg ground heat flux [w/m2]  [+ to soil]
    REAL                                :: IRB          ! bare net longwave rad. [w/m2] [+ to atm]
    REAL                                :: SHB          ! bare sensible heat [w/m2]     [+ to atm]
    REAL                                :: EVB          ! bare evaporation heat [w/m2]  [+ to atm]
    REAL                                :: GHB          ! bare ground heat flux [w/m2] [+ to soil]
    REAL                                :: TR           ! transpiration [w/m2]  [+ to atm]
    REAL                                :: EVC          ! canopy evaporation heat [w/m2]  [+ to atm]
    REAL                                :: CHLEAF       ! leaf exchange coefficient 
    REAL                                :: CHUC         ! under canopy exchange coefficient 
    REAL                                :: CHV2         ! veg 2m exchange coefficient 
    REAL                                :: CHB2         ! bare 2m exchange coefficient 

! Intermediate terms

    REAL                                :: FPICE        ! snow fraction of precip
    REAL                                :: FCEV         ! canopy evaporation heat (w/m2) [+ to atm]
    REAL                                :: FGEV         ! ground evaporation heat (w/m2) [+ to atm]
    REAL                                :: FCTR         ! transpiration heat flux (w/m2) [+ to atm]
    REAL                                :: QSNBOT       ! snowmelt out bottom of pack [mm/s]
    REAL                                :: PONDING      ! snowmelt with no pack [mm]
    REAL                                :: PONDING1     ! snowmelt with no pack [mm]
    REAL                                :: PONDING2     ! snowmelt with no pack [mm]

! Local terms

    REAL                                :: FSR          ! total reflected solar radiation (w/m2)
    REAL, DIMENSION(-2:0)               :: FICEOLD      ! snow layer ice fraction []
    REAL                                :: CO2PP        ! CO2 partial pressure [Pa]
    REAL                                :: O2PP         ! O2 partial pressure [Pa]
    REAL, DIMENSION(1:NSOIL)            :: ZSOIL        ! depth to soil interfaces [m]
    REAL                                :: FOLN         ! nitrogen saturation [%]

    REAL                                :: QC           ! cloud specific humidity for MYJ [not used]
    REAL                                :: PBLH         ! PBL height for MYJ [not used]
    REAL                                :: DZ8W1D       ! model level heights for MYJ [not used]

    INTEGER                             :: I
    INTEGER                             :: J
    INTEGER                             :: K
    INTEGER                             :: ICE
    INTEGER                             :: SLOPETYP
    LOGICAL                             :: IPRINT

    INTEGER                             :: ISC          ! soil color index
    INTEGER                             :: IST          ! surface type 1-soil; 2-lake
    INTEGER                             :: YEARLEN
    INTEGER                             :: ITIME
    REAL,intent(inout) :: si0,si1
    REAL :: QSPRING,zwt0,minzwt
    REAL, DIMENSION(1:NSOIL) :: SLDPTH
    REAL :: DZWT
    REAL :: ERRWAT
    INTEGER :: ROOT_DEPTH
    !Calculate cosz
    CALL CALC_DECLIN (NOWDATE, LAT, LON, COSZ)

    !Set parameters
    CALL REDPRM (VEGTYP,SOILTYP,SLOPETYP,ZSOIL,NSOIL,ISURBAN)

    ROOT_DEPTH = NROOT 

    if (itime .eq. 0) then 
      !Define the equilibrium soil water content
      CALL EQSMOISTURE(NSOIL ,  ZSOIL , SMCMAX , SMCWLT ,DWSAT, DKSAT  ,BEXP  , SMCEQ)
    end if

    !Update the water table depth
    call UPDATEWTD  (NSOIL,SLDPTH,SMCEQ,SMCMAX,SMCWLT,-PSISAT,BEXP,I,J,DZWT,ZWT ,SMC, SMH2O ,SMCWTD,QSPRING,&
                     DT,DEEPRECH,DKSAT)
    DEEPRECH = 0.0  

    !Calculate the soil moisture deficit
    zwt0 = zwt
    !call Calculate_Deficit(si0,smcmax,smcwtd,zwt,sldpth,smc,nsoil)
    CALL NOAHMP_SFLX (&
            I       , J       , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
            DT      , DX      , DZ8W1D  , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration 
            FVEG    , FVGMAX  , VEGTYP  , ISURBAN , ICE     , IST     , & ! IN : Vegetation/Soil characteristics
            ISC     , SMCEQ   ,                                         & ! IN : Vegetation/Soil characteristics
            IZ0TLND ,                                                   & ! IN : User options
            T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
            QC      , SWDN    , LWDN    , PRCP    , TBOT    , CO2PP   , & ! IN : Forcing
            O2PP    , FOLN    , FICEOLD , PBLH    , Z_ML    ,           & ! IN : Forcing
            HDIV, & ! IN : Horizontal
            ALBOLD  , SNEQVO  ,                                         & ! IN/OUT : 
            STC     , SMH2O   , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT : 
            CANLIQ  , CANICE  , TV      , TG      , QSFC1D  , QSNOW   , & ! IN/OUT : 
            ISNOW   , ZSNSO   , SNDPTH  , SWE     , SNICE   , SNLIQ   , & ! IN/OUT : 
            ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT : 
            STMASS  , WOOD    , STBLCP  , FASTCP  , PLAI    , PSAI    , & ! IN/OUT : 
            CM      , CH      , TAUSS   ,                               & ! IN/OUT : 
            SMCWTD  ,DEEPRECH , RECH    ,                               & ! IN/OUT :
            FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT : 
            FGEV    , FCTR    , ECAN    , ETRAN   , ESOIL   , TRAD    , & ! OUT : 
            TGB     , TGV     , T2MV    , T2MB    , Q2MV    , Q2MB    , & ! OUT : 
            RUNSF   , RUNSB   , APAR    , PSN     , SAV     , SAG     , & ! OUT : 
            FSNO    , NEE     , GPP     , NPP     , FVEGMP  , SALB    , & ! OUT : 
            QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT : 
            BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT : 
            SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
	    GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
	    CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE)

     !Add in QSPRING
     runsb = runsb+1000*qspring/dt! (Dynamic TOPMODEL ONLY)

end subroutine run_model_cell

!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM( MMINLU, MMINSL, VEGPARM_FILE, SOILPARM_FILE,GENPARM_FILE,&
           BB0,DRYSMC0,F110,MAXSMC0,REFSMC0,SATPSI0,SATDK0,SATDW0,WLTSMC0,QTZ0,NHRU)
!-----------------------------------------------------------------

  USE module_sf_noahmplsm
  IMPLICIT NONE

  integer :: nhru
  real,dimension(nhru) :: bb0,drysmc0,f110,maxsmc0,refsmc0,satpsi0,satdk0
  real,dimension(nhru) :: satdw0,wltsmc0,qtz0
  !real,dimension(nhru) :: nrotbl ! Noemi
  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  character(len=256) :: vegparm_file,soilparm_file,genparm_file
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message
  logical, external :: wrf_dm_on_monitor

!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  IF ( wrf_dm_on_monitor() ) THEN

     OPEN(19, FILE=VEGPARM_FILE,FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF


     LUMATCH=0

     FIND_LUTYPE : DO WHILE (LUMATCH == 0)
        READ (19,*,END=2002)
        READ (19,*,END=2002)LUTYPE
        READ (19,*)LUCATS,IINDEX

        IF(LUTYPE.EQ.MMINLU)THEN
           !WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
           !CALL wrf_message( mess )
           LUMATCH=1
        ELSE
           !call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
           DO LC = 1, LUCATS+12
              read(19,*)
           ENDDO
        ENDIF
     ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(SHDTBL)       < LUCATS .OR. &
          SIZE(NROTBL)       < LUCATS .OR. &
          SIZE(RSTBL)        < LUCATS .OR. &
          SIZE(RGLTBL)       < LUCATS .OR. &
          SIZE(HSTBL)        < LUCATS .OR. &
          SIZE(SNUPTBL)      < LUCATS .OR. &
          SIZE(MAXALB)       < LUCATS .OR. &
          SIZE(LAIMINTBL)    < LUCATS .OR. &
          SIZE(LAIMAXTBL)    < LUCATS .OR. &
          SIZE(Z0MINTBL)     < LUCATS .OR. &
          SIZE(Z0MAXTBL)     < LUCATS .OR. &
          SIZE(ALBEDOMINTBL) < LUCATS .OR. &
          SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
          SIZE(EMISSMINTBL ) < LUCATS .OR. &
          SIZE(EMISSMAXTBL ) < LUCATS ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
     ENDIF

     IF(LUTYPE.EQ.MMINLU)THEN
        DO LC=1,LUCATS
           READ (19,*)IINDEX,SHDTBL(LC),                        &
                NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
                SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
                LAIMAXTBL(LC),EMISSMINTBL(LC),             &
                EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
                ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
        ENDDO
!
        READ (19,*)
        READ (19,*)TOPT_DATA
        READ (19,*)
        READ (19,*)CMCMAX_DATA
        READ (19,*)
        READ (19,*)CFACTR_DATA
        READ (19,*)
        READ (19,*)RSMAX_DATA
        READ (19,*)
        READ (19,*)BARE
        READ (19,*)
        READ (19,*)NATURAL
     ENDIF
!
2002 CONTINUE

     CLOSE (19)
     IF (LUMATCH == 0) then
        CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
     ENDIF
  ENDIF

  
  CALL wrf_dm_bcast_string  ( LUTYPE  , 4 )
  CALL wrf_dm_bcast_integer ( LUCATS  , 1 )
  CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
  CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
  CALL wrf_dm_bcast_real    ( SHDTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( NROTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( RSTBL   , NLUS )
  CALL wrf_dm_bcast_real    ( RGLTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( HSTBL   , NLUS )
  CALL wrf_dm_bcast_real    ( SNUPTBL , NLUS )
  CALL wrf_dm_bcast_real    ( LAIMINTBL    , NLUS )
  CALL wrf_dm_bcast_real    ( LAIMAXTBL    , NLUS )
  CALL wrf_dm_bcast_real    ( Z0MINTBL     , NLUS )
  CALL wrf_dm_bcast_real    ( Z0MAXTBL     , NLUS )
  CALL wrf_dm_bcast_real    ( EMISSMINTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( EMISSMAXTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( ALBEDOMINTBL , NLUS )
  CALL wrf_dm_bcast_real    ( ALBEDOMAXTBL , NLUS )
  CALL wrf_dm_bcast_real    ( MAXALB  , NLUS )
  CALL wrf_dm_bcast_real    ( TOPT_DATA    , 1 )
  CALL wrf_dm_bcast_real    ( CMCMAX_DATA  , 1 )
  CALL wrf_dm_bcast_real    ( CFACTR_DATA  , 1 )
  CALL wrf_dm_bcast_real    ( RSMAX_DATA  , 1 )
  CALL wrf_dm_bcast_integer ( BARE    , 1 )
  CALL wrf_dm_bcast_integer ( NATURAL , 1 )

  
!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
     !print*,'NoahMP: Defining the soil properties'
     slcats = nhru
     bb(1:nhru) = bb0(1:nhru)
     drysmc(1:nhru) = drysmc0(1:nhru)
     f11(1:nhru) = f110(1:nhru)
     maxsmc(1:nhru) = maxsmc0(1:nhru)
     refsmc(1:nhru) = refsmc0(1:nhru)
     satpsi(1:nhru) = satpsi0(1:nhru)
     satdk(1:nhru) = satdk0(1:nhru)
     satdw(1:nhru) = satdw0(1:nhru)
     wltsmc(1:nhru) = wltsmc0(1:nhru)
     qtz(1:nhru) = qtz0(1:nhru)

  CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
  CALL wrf_dm_bcast_string  ( SLTYPE  , 4 )
  CALL wrf_dm_bcast_string  ( MMINSL  , 4 )  ! since this is reset above, see oct2 ^
  CALL wrf_dm_bcast_integer ( SLCATS  , 1 )
  CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
  CALL wrf_dm_bcast_real    ( BB      , NSLTYPE )
  CALL wrf_dm_bcast_real    ( DRYSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( F11     , NSLTYPE )
  CALL wrf_dm_bcast_real    ( MAXSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( REFSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATPSI  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDK   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDW   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( WLTSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( QTZ     , NSLTYPE )

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIl TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE=GENPARM_FILE,FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

     READ (19,*)
     READ (19,*)
     READ (19,*) NUM_SLOPE

     SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
        CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
     ENDIF

     DO LC=1,SLPCATS
        READ (19,*)SLOPE_DATA(LC)
     ENDDO

     READ (19,*)
     READ (19,*)SBETA_DATA
     READ (19,*)
     READ (19,*)FXEXP_DATA
     READ (19,*)
     READ (19,*)CSOIL_DATA
     READ (19,*)
     READ (19,*)SALP_DATA
     READ (19,*)
     READ (19,*)REFDK_DATA
     READ (19,*)
     READ (19,*)REFKDT_DATA
     READ (19,*)
     READ (19,*)FRZK_DATA
     READ (19,*)
     READ (19,*)ZBOT_DATA
     READ (19,*)
     READ (19,*)CZIL_DATA
     READ (19,*)
     READ (19,*)SMLOW_DATA
     READ (19,*)
     READ (19,*)SMHIGH_DATA
     READ (19,*)
     READ (19,*)LVCOEF_DATA
     CLOSE (19)
  ENDIF

  CALL wrf_dm_bcast_integer ( NUM_SLOPE    ,  1 )
  CALL wrf_dm_bcast_integer ( SLPCATS      ,  1 )
  CALL wrf_dm_bcast_real    ( SLOPE_DATA   ,  NSLOPE )
  CALL wrf_dm_bcast_real    ( SBETA_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( FXEXP_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( CSOIL_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SALP_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( REFDK_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( REFKDT_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( FRZK_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( ZBOT_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( CZIL_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( SMLOW_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SMHIGH_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( LVCOEF_DATA  ,  1 )


!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------
!
logical function wrf_dm_on_monitor() result(l)
  l = .TRUE.
  return
end function wrf_dm_on_monitor
!
SUBROUTINE CALC_DECLIN ( NOWDATE, LATITUDE, LONGITUDE, COSZ)

  USE MODULE_DATE_UTILITIES
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

   REAL, PARAMETER :: DEGRAD = 3.14159265/180.
   REAL, PARAMETER :: DPD    = 360./365.
   REAL, PARAMETER :: RADDEG = 180./3.14159265
! !ARGUMENTS:
   CHARACTER(LEN=19), INTENT(IN)  :: NOWDATE    ! YYYY-MM-DD_HH:MM:SS
   REAL,              INTENT(IN)  :: LATITUDE
   REAL,              INTENT(IN)  :: LONGITUDE
   REAL,              INTENT(OUT) :: COSZ
   REAL                           :: JULIAN
   REAL                           :: HRANG
   REAL                           :: DECLIN
   REAL                           :: OBECL
   REAL                           :: SINOB
   REAL                           :: SXLONG
   REAL                           :: ARG
   REAL                           :: TLOCTIM
   INTEGER                        :: IDAY
   INTEGER                        :: IHOUR
   INTEGER                        :: IMINUTE
   INTEGER                        :: ISECOND

   CALL GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
   READ(NOWDATE(12:13), *) IHOUR
   READ(NOWDATE(15:16), *) IMINUTE
   READ(NOWDATE(18:19), *) ISECOND
   JULIAN = REAL(IDAY) + REAL(IHOUR)/24.
!
! FOR SHORT WAVE RADIATION

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + RADDEG*LONGITUDE/15.0 ! LOCAL TIME IN HOURS
   TLOCTIM = AMOD(TLOCTIM+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE)*SIN(DECLIN)+COS(LATITUDE)*COS(DECLIN)*COS(HRANG)

 END SUBROUTINE CALC_DECLIN

! ----------------------------------------------------------------------
 SUBROUTINE EQSMOISTURE(NSOIL  ,  ZSOIL , SMCMAX , SMCWLT, DWSAT , DKSAT ,BEXP , & !in
                         SMCEQ                          )  !out
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
  REAL,                            INTENT(IN) :: SMCMAX , SMCWLT, BEXP , DWSAT, DKSAT
!output
  REAL,  DIMENSION(      1:NSOIL), INTENT(OUT) :: SMCEQ  !equilibrium soil water  content [m3/m3]
!local
  INTEGER                                     :: K , ITER
  REAL                                        :: DDZ , SMC, FUNC, DFUNC , AA, BB , EXPON, DX

!gmmcompute equilibrium soil moisture content for the layer when wtd=zsoil(k)


   DO K=1,NSOIL

            IF ( K == 1 )THEN
                DDZ = -ZSOIL(K+1) * 0.5
            ELSEIF ( K < NSOIL ) THEN
                DDZ = ( ZSOIL(K-1) - ZSOIL(K+1) ) * 0.5
            ELSE
                DDZ = ZSOIL(K-1) - ZSOIL(K)
            ENDIF

!use Newton-Raphson method to find eq soil moisture

            EXPON = BEXP +1.
            AA = DWSAT/DDZ
            BB = DKSAT / SMCMAX ** EXPON

            SMC = 0.5 * SMCMAX

         DO ITER = 1, 100
            FUNC = (SMC - SMCMAX) * AA +  BB * SMC ** EXPON
            DFUNC = AA + BB * EXPON * SMC ** BEXP

            DX = FUNC/DFUNC
            SMC = SMC - DX
            IF ( ABS (DX) < 1.E-6)EXIT
         ENDDO

!             SMCEQ(K) = MIN(MAX(SMC,SMCWLT),SMCMAX*0.99)
             SMCEQ(K) = MIN(MAX(SMC,1.E-4),SMCMAX*0.99)
   ENDDO

 END  SUBROUTINE EQSMOISTURE

! ==================================================================================================
! ----------------------------------------------------------------------
  SUBROUTINE UPDATEWTD  (NSOIL,  DZS,  SMCEQ                ,& !in
                         SMCMAX, SMCWLT, PSISAT, BEXP ,ILOC ,JLOC  ,& !in
                         TOTWATER, WTD ,SMC, SH2O ,SMCWTD          ,& !inout
                         QSPRING                                 ,& !out
                         DELTAT,DEEPRECH,DKSAT)  !other
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  INTEGER,                         INTENT(IN) :: ILOC, JLOC
  REAL,                         INTENT(IN)    :: SMCMAX
  REAL,                         INTENT(IN)    :: SMCWLT
  REAL,                         INTENT(IN)    :: PSISAT
  REAL,                         INTENT(IN)    :: BEXP
  REAL,  DIMENSION(       0:NSOIL) :: ZSOIL !depth of soil layer-bottom [m]
  REAL,  DIMENSION(       1:NSOIL), INTENT(IN) :: SMCEQ  !equilibrium soil water  content [m3/m3]
  REAL,  DIMENSION(       1:NSOIL), INTENT(IN) :: DZS ! soil layer thickness [m]
! input-output
  REAL                           , INTENT(INOUT) :: TOTWATER
  REAL                           , INTENT(INOUT) :: WTD
  REAL                           , INTENT(INOUT) :: SMCWTD
  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC
  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O
! output
  REAL                           , INTENT(OUT) :: QSPRING
!local
  INTEGER                                     :: K
  INTEGER                                     :: K1
  INTEGER                                     :: IWTD
  INTEGER                                     :: KWTD
  REAL                                        :: MAXWATUP, MAXWATDW ,WTDOLD
  REAL                                        :: WGPMID
  REAL                                        :: SYIELDDW
  REAL                                        :: DZUP
  REAL                                        :: SMCEQDEEP
  REAL, DIMENSION(       1:NSOIL)             :: SICE
  REAL :: DELTAT,DEEPRECH,DDZ,SMCWTDMID,PSI,WCNDDEEP,WFLUXDEEP,WPLUS,WMINUS,DKSAT,wb0,wb1,totwater0
!deep recharge variables
! -------------------------------------------------------------

  !Define depths
  ZSOIL(0) = 0.
  ZSOIL(1) = -DZS(1)
  DO K = 2, NSOIL
    ZSOIL(K)         = -DZS(K) + ZSOIL(K-1)
  END DO

  !deeprech = 0.0
  !for deep water table calculate recharge
  !IF(WTD < ZSOIL(NSOIL)-DZS(NSOIL))THEN
  !!IF(1 .eq. 2)THEN!WTD < ZSOIL(NSOIL)-DZS(NSOIL))THEN
  IF(WTD < ZSOIL(NSOIL)-DZS(NSOIL))THEN
  !IF(WTD < ZSOIL(NSOIL))THEN
  ! !assume all liquid if the wtd is deep
   DDZ = ZSOIL(NSOIL)-WTD
   !print*,WTD,ZSOIL(NSOIL)-DZS(NSOIL),DDZ
   SMCWTDMID = 0.5 * (SMCWTD + SMCMAX )
   PSI = PSISAT * ( SMCMAX / SMCWTD ) ** BEXP
   WCNDDEEP = DKSAT * ( SMCWTDMID / SMCMAX ) ** (2.0*BEXP + 3.0)
   WFLUXDEEP =  - DELTAT * WCNDDEEP * ( (PSISAT-PSI) / DDZ - 1.)
   !update deep soil moisture
   SMCWTD = SMCWTD  + (DEEPRECH -  WFLUXDEEP)  / DDZ
   WPLUS       = MAX((SMCWTD-SMCMAX), 0.0) * DDZ
   WMINUS       = MAX((1.E-4-SMCWTD), 0.0) * DDZ
   WFLUXDEEP = WFLUXDEEP + WPLUS - WMINUS
   DEEPRECH = WFLUXDEEP
   !if ((WTD .lt. -2.0) .and. (DEEPRECH .lt. 0.0))then
   ! DEEPRECH = 0.0
   !else 
   SMCWTD = MAX( MIN(SMCWTD,SMCMAX) , 1.E-4)
   !endif
  ENDIF
  
 
  !Add deeprech to the dzwt 
  totwater = totwater + deeprech
  totwater0 = totwater

 
  QSPRING=0.

  SICE = SMC - SH2O

iwtd=1
!if (totwater .ne. 0)print*,'part1',abs(sum(smc*zsoil(1:nsoil))),wtd,smcwtd*(zsoil(nsoil)-wtd)
!wb0 = sum(dzs*smc)+smcwtd*(-wtd - sum(dzs*smc))+smcmax*(100000+wtd)
wb0 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs)))
!wb0 = sum(dzs*smc)!+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(100000-max(-wtd,sum(dzs)))

!case 1: totwater > 0 (water table going up):
IF(totwater.gt.0.)then


         if(wtd.ge.zsoil(nsoil))then

            do k=nsoil-1,1,-1
              if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k
            kwtd=iwtd+1

            !max water that fits in the layer
            maxwatup=dzs(kwtd)*(smcmax-smc(kwtd))

            if(totwater.le.maxwatup)then
               smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
               smc(kwtd) = min(smc(kwtd),smcmax)
               if(smc(kwtd).gt.smceq(kwtd))wtd = min ( ( smc(kwtd)*dzs(kwtd) &
                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                     ( smcmax-smceq(kwtd) ) , zsoil(iwtd) )
               totwater=0.
            else   !water enough to saturate the layer
              smc(kwtd) = smcmax
              totwater=totwater-maxwatup
              k1=iwtd
              do k=k1,0,-1
                 wtd = zsoil(k)
                 iwtd=k-1
                 if(k.eq.0)exit
                 maxwatup=dzs(k)*(smcmax-smc(k))
                 if(totwater.le.maxwatup)then
                   smc(k) = smc(k) + totwater / dzs(k)
                   smc(k) = min(smc(k),smcmax)
                   if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                     - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                     ( smcmax-smceq(k) ) , zsoil(iwtd) )
                   totwater=0.
                   exit
                 else
                    smc(k) = smcmax
                    totwater=totwater-maxwatup
                 endif

              enddo

            endif

         elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then ! wtd below bottom of soil model

            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
                smceqdeep = max(smceqdeep,1.E-4)

            maxwatup=(smcmax-smcwtd)*dzs(nsoil)

            if(totwater.le.maxwatup)then
                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                if(smcwtd.gt.smceqdeep)wtd = min( ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep ) , zsoil(nsoil) )
                totwater=0.
            else
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))
                    if(totwater.le.maxwatup)then
                     smc(k) = min(smc(k) + totwater / dzs(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = min ( ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) ) , zsoil(iwtd) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                enddo
             endif

!deep water table
       else

            maxwatup=(smcmax-smcwtd)*(zsoil(nsoil)-dzs(nsoil)-wtd)
            if(totwater.le.maxwatup)then
               wtd = wtd + totwater/(smcmax-smcwtd)
               totwater=0.
            else
               totwater=totwater-maxwatup
               wtd=zsoil(nsoil)-dzs(nsoil)
               maxwatup=(smcmax-smcwtd)*dzs(nsoil)
              if(totwater.le.maxwatup)then

            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
                smceqdeep = max(smceqdeep,1.E-4)

                smcwtd = smcwtd + totwater / dzs(nsoil)
                smcwtd = min(smcwtd,smcmax)
                wtd = ( smcwtd*dzs(nsoil) &
                 - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                     ( smcmax-smceqdeep )
                totwater=0.
                smcwtd=smcmax
                totwater=totwater-maxwatup
                do k=nsoil,0,-1
                    wtd=zsoil(k)
                    iwtd=k-1
                    if(k.eq.0)exit
                    maxwatup=dzs(k)*(smcmax-smc(k))

                    if(totwater.le.maxwatup)then
                     smc(k) = smc(k) + totwater / dzs(k)
                     smc(k) = min(smc(k),smcmax)
                     if(smc(k).gt.smceq(k))wtd = ( smc(k)*dzs(k) &
                        - smceq(k)*zsoil(iwtd) + smcmax*zsoil(k) ) / &
                           ( smcmax-smceq(k) )
                     totwater=0.
                     exit
                    else
                     smc(k) = smcmax
                     totwater=totwater-maxwatup
                    endif
                   enddo
               endif
             endif
         endif

!water springing at the surface
        qspring=totwater
        wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs)))
        !print*,'check1',1000*(wb1 - wb0 - totwater0)

!case 2: totwater < 0 (water table going down):
ELSEIF(totwater.lt.0.)then

        wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs)))
        !print*,'check2.1',1000*(wb1 - wb0)
         if(wtd.ge.zsoil(nsoil))then !wtd in the resolved layers

            do k=nsoil-1,1,-1 !WHY nsoil -1???
               if(wtd.lt.zsoil(k))exit
            enddo
            iwtd=k

               k1=iwtd+1
               do kwtd=k1,nsoil

!max water that the layer can yield
                  maxwatdw=dzs(kwtd)*(smc(kwtd)-max(smceq(kwtd),sice(kwtd)))

                  if(-totwater.le.maxwatdw)then
                        smc(kwtd) = smc(kwtd) + totwater / dzs(kwtd)
                        if(smc(kwtd).gt.smceq(kwtd))then
                              wtd = ( smc(kwtd)*dzs(kwtd) &
                                 - smceq(kwtd)*zsoil(iwtd) + smcmax*zsoil(kwtd) ) / &
                                 ( smcmax-smceq(kwtd) )
                         else
                              wtd=zsoil(kwtd)
                              iwtd=iwtd+1
                         endif
                         totwater=0.
                         wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs))) + totwater - totwater0
                         !print*,'check2.105',1000*(wb1 - wb0)
                         exit
                   else
                         wtd = zsoil(kwtd)
                         iwtd=iwtd+1
                         if(maxwatdw.ge.0.)then
                            smc(kwtd) = smc(kwtd) - maxwatdw / dzs(kwtd)
                            totwater = totwater + maxwatdw
                         endif
                   endif

                enddo
        wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs))) + totwater - totwater0
        !wb1 = sum(dzs*smc) + totwater - totwater0
        !print*,'check2.11',1000*(wb1 - wb0)

               if(iwtd.eq.nsoil.and.totwater.lt.0.)then
            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

                  maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

                  if(-totwater.le.maxwatdw)then
                       smcwtd = smcwtd + totwater / dzs(nsoil)
                       wtd = min( ( smcwtd*dzs(nsoil) &
                           - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                            ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )
                       totwater = 0.0

                  else
                   !wb0 = sum(dzs*smc)+smcwtd*(-wtd - sum(dzs*smc))+smcmax*(10000+wtd)
                       !wtd=zsoil(nsoil)-dzs(nsoil)
                       !wtd=zsoil(nsoil)+dzs(nsoil)
                       !dzup = abs(totwater)/(smcmax-smceqdeep)
                       !smcwtd = smcwtd + totwater / dzs(nsoil)
!and now even further down
                       !dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
                       dzup = abs(totwater)/(smcmax-smceqdeep)
                       !totwater = totwater + (smceqdeep-smcwtd)*dzs(nsoil)
                       totwater = 0.0
                       wtd=wtd-dzup
                       smcwtd=smceqdeep
                       !print*,smcwtd,wtd,dzup

                  endif

                endif
        wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs))) + totwater - totwater0
        !print*,'check2.2',1000*(wb1 - wb0)


        elseif(wtd.ge.zsoil(nsoil)-dzs(nsoil))then

!if wtd was already below the bottom of the resolved soil crust
            !gmmequilibrium soil moisture content
               smceqdeep = smcmax * ( psisat / &
                           (psisat - dzs(nsoil)) ) ** (1./bexp)
!               smceqdeep = max(smceqdeep,smcwlt)
               smceqdeep = max(smceqdeep,1.E-4)

            maxwatdw=dzs(nsoil)*(smcwtd-smceqdeep)

            if(-totwater.le.maxwatdw)then

               smcwtd = smcwtd + totwater / dzs(nsoil)
               wtd = min( ( smcwtd*dzs(nsoil) & !MAX
                    - smceqdeep*zsoil(nsoil) + smcmax*(zsoil(nsoil)-dzs(nsoil)) ) / &
                    ( smcmax-smceqdeep ) , zsoil(nsoil)-dzs(nsoil) )

            else

               !wtd=zsoil(nsoil)-dzs(nsoil)
               !smcwtd = smcwtd + totwater / dzs(nsoil)
!and now even further down
               !dzup=(smceqdeep-smcwtd)*dzs(nsoil)/(smcmax-smceqdeep)
               !wtd=wtd-dzup
               !smcwtd=smceqdeep
               dzup = abs(totwater)/(smcmax-smceqdeep)
               totwater = 0.0
               wtd=wtd-dzup
               smcwtd=smceqdeep

             endif

         else
!gmmequilibrium soil moisture content
               wgpmid = smcmax * ( psisat / &
                    (psisat - (zsoil(nsoil)-wtd)) ) ** (1./bexp)
               !wgpmid=max(wgpmid,smcwlt)
               wgpmid=max(wgpmid,1.E-4)
               syielddw=smcmax-wgpmid
               wtdold=wtd
               wtd = wtdold + totwater/syielddw
!update wtdwgp
               smcwtd = (smcwtd*(zsoil(nsoil)-wtdold)+wgpmid*(wtdold-wtd) ) / (zsoil(nsoil)-wtd)

          endif

          qspring=0.

ENDIF

         SH2O = SMC - SICE
         wb1 = sum(dzs*smc)+smcwtd*(max(-wtd - sum(dzs),0.0))+smcmax*(1000.0-max(-wtd,sum(dzs)))
         !print*,'final',1000*(wb1 - wb0 - totwater0 + qspring),1000*totwater0,1000*totwater,-wtd
  
END SUBROUTINE UPDATEWTD

module model

 implicit none

 !Declare variables

 !Miscellanous parameters
 real :: badval = 0.0!-1.E36 
 real :: co2 = 395.e-06
 real :: o2 = 0.209
  
 !General parameters
 integer :: nsoil ! Number of soil layers
 integer :: ncells ! Number of grid cells
 integer :: nsnow = 3 ! Number of snow layers
 real :: dt ! Timestep in seconds
 real :: dx ! Spatial resolution (m)
 real :: julian ! Julian day
 integer :: yearlen
 logical :: iz0tlnd
 integer :: itime
 character(len=256) :: llanduse,lsoil
 character(len=256) :: vegparm_file,soilparm_file,genparm_file,mptable_file
 character(len=19) :: nowdate

 !Model options
 integer :: idveg ! dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1      
 integer :: iopt_crs ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
 integer :: iopt_btr ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
 integer :: iopt_run ! runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
 integer :: iopt_sfc ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
 integer :: iopt_frz ! supercooled liquid water (1-> NY06; 2->Koren99)
 integer :: iopt_inf ! frozen soil permeability (1-> NY06; 2->Koren99)
 integer :: iopt_rad ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
 integer :: iopt_alb ! snow surface albedo (1->BATS; 2->CLASS)
 integer :: iopt_snf ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)]
 integer :: iopt_tbot ! lower boundary of soil temperature (1->zero-flux; 2->Noah) 
 integer :: iopt_stc ! snow/soil temperature time scheme (only layer 1) 1 -> semi-implicit; 2 -> full implicit (original Noah)

 !Single layer
 integer,dimension(:),allocatable :: isnow ! actual number of snow layers
 real,dimension(:),allocatable :: z_ml ! model height [m]
 real,dimension(:),allocatable :: lwdn       ! Downward longwave radiation flux at surface [W/m2]
 real,dimension(:),allocatable :: swdn      ! Downward shortwave radiation flux at surface [W/m2]
 real,dimension(:),allocatable :: p_ml    !  Surface pressure defined at intermediate level [Pa]
 real,dimension(:),allocatable :: psfc    ! Surface pressure [Pa]
 real,dimension(:),allocatable :: prcp       ! Precipitation rate (kg m-2 s-1)
 real,dimension(:),allocatable :: t_ml     ! Air temperature (K) [Forcing]
 real,dimension(:),allocatable :: q_ml         ! Surface specific humidity (kg kg-1)
 real,dimension(:),allocatable :: u_ml       ! West-to-east component of the surface [m/s]
 real,dimension(:),allocatable :: v_ml       ! North-to-south component of the surface [m/s]
 real,dimension(:),allocatable :: fsh ! total sensible heat (w/m2)
 real,dimension(:),allocatable :: ssoil ! soil heat (w/m2)
 real,dimension(:),allocatable :: salb ! surface albedo
 real,dimension(:),allocatable :: fsno ! snow cover fraction
 real,dimension(:),allocatable :: swe ! snow water equivalent (mm)
 real,dimension(:),allocatable :: sndpth ! snow depth (m)
 real,dimension(:),allocatable :: emissi ! net surface emissivity
 real,dimension(:),allocatable :: qsfc1d ! bulk surface specific humidity (kg/kg)
 real,dimension(:),allocatable :: tv ! vegetation canopy temperature
 real,dimension(:),allocatable :: tg ! ground surface temperature
 real,dimension(:),allocatable :: canice ! canopy-intercepted ice (mm)
 real,dimension(:),allocatable :: canliq ! canopy-intercepted liquid water (mm)
 real,dimension(:),allocatable :: eah ! canopy air vapor pressure (pa)
 real,dimension(:),allocatable :: tah ! canopy air temperature (K)
 real,dimension(:),allocatable :: cm ! momentum drag coefficient
 real,dimension(:),allocatable :: ch ! sensible heat exchange coefficient
 real,dimension(:),allocatable :: fwet ! weted or snowed fraction of the canopy
 real,dimension(:),allocatable :: sneqvo ! snow mass at last time step (mm h2o)
 real,dimension(:),allocatable :: albold ! snow albedo at last time step
 real,dimension(:),allocatable :: qsnow ! snowfall on the ground (mm/s)
 real,dimension(:),allocatable :: wslake ! lake water storage (mm)
 real,dimension(:),allocatable :: zwt ! water table depth (m)
 real,dimension(:),allocatable :: dzwt ! change in water table depth (m)
 real,dimension(:),allocatable :: wa ! water in the aquifer (mm)
 real,dimension(:),allocatable :: wt ! groundwater storage (mm)
 real,dimension(:),allocatable :: smcwtd ! soil moisture content in the transmission zone
 real,dimension(:),allocatable :: deeprech ! recharge to the water table when deep
 real,dimension(:),allocatable :: rech ! recharge to the water table (diagnostic)
 real,dimension(:),allocatable :: lfmass ! leaf mass (g/m2)
 real,dimension(:),allocatable :: rtmass ! mass of fine roots (g/m2)
 real,dimension(:),allocatable :: stmass ! stem mass (g/m2)
 real,dimension(:),allocatable :: wood ! mass of wood (incl. woody roots) (g/m2)
 real,dimension(:),allocatable :: stblcp ! stable carbon in deep soil (g/m2)
 real,dimension(:),allocatable :: fastcp ! short lived carbon, shallow soil (g/m2)
 real,dimension(:),allocatable :: plai ! leaf area index
 real,dimension(:),allocatable :: psai ! stem area index
 real,dimension(:),allocatable :: tauss ! non-dimensional snow age
 real,dimension(:),allocatable :: t2mv ! 2m temperature of vegetation part
 real,dimension(:),allocatable :: t2mb ! 2m temperature of bare ground part
 real,dimension(:),allocatable :: q2mv ! 2m mixing ratio of vegetation part
 real,dimension(:),allocatable :: q2mb ! 2m mixing ratio of bare ground part
 real,dimension(:),allocatable :: trad ! surface radiative temperature (k)
 real,dimension(:),allocatable :: nee ! net ecosys exchange (g/m2/s CO2)
 real,dimension(:),allocatable :: gpp ! gross primary assimilation [g/m2/s C]
 real,dimension(:),allocatable :: npp ! net primary productivity [g/m2/s C]
 real,dimension(:),allocatable :: fvegmp ! greenness vegetation fraction [-]
 real,dimension(:),allocatable :: runsf ! surface runoff [mm/s]
 real,dimension(:),allocatable :: runsb ! subsurface runoff [mm/s]
 real,dimension(:),allocatable :: ecan ! evaporation of intercepted water (mm/s)
 real,dimension(:),allocatable :: etran ! transpiration rate (mm/s)
 real,dimension(:),allocatable :: esoil ! soil surface evaporation rate (mm/s]
 real,dimension(:),allocatable :: fsa ! total absorbed solar radiation (w/m2)
 real,dimension(:),allocatable :: fira ! total net longwave rad (w/m2) [+ to atm]
 real,dimension(:),allocatable :: apar ! photosyn active energy by canopy (w/m2)
 real,dimension(:),allocatable :: psn ! total photosynthesis (umol co2/m2/s) [+]
 real,dimension(:),allocatable :: sav ! solar rad absorbed by veg. (w/m2)
 real,dimension(:),allocatable :: sag ! solar rad absorbed by ground (w/m2)
 real,dimension(:),allocatable :: rssun ! sunlit leaf stomatal resistance (s/m)
 real,dimension(:),allocatable :: rssha ! shaded leaf stomatal resistance (s/m)
 real,dimension(:),allocatable :: bgap ! between gap fraction
 real,dimension(:),allocatable :: wgap ! within gap fraction
 real,dimension(:),allocatable :: tgv ! under canopy ground temperature [K]
 real,dimension(:),allocatable :: tgb ! bare ground temperature [K]
 real,dimension(:),allocatable :: chv ! sensible heat exchange coefficient vegetated
 real,dimension(:),allocatable :: chb ! sensible heat exchange coefficient bare-ground
 real,dimension(:),allocatable :: irc ! canopy net LW rad. [w/m2] [+ to atm]
 real,dimension(:),allocatable :: irg ! veg ground net LW rad. [w/m2] [+ to atm]
 real,dimension(:),allocatable :: shc ! canopy sen. heat [w/m2]   [+ to atm]
 real,dimension(:),allocatable :: shg ! veg ground sen. heat [w/m2]   [+ to atm]
 real,dimension(:),allocatable :: evg ! veg ground evap. heat [w/m2]  [+ to atm]
 real,dimension(:),allocatable :: ghv ! veg ground heat flux [w/m2]  [+ to soil]
 real,dimension(:),allocatable :: irb ! bare net longwave rad. [w/m2] [+ to atm]
 real,dimension(:),allocatable :: shb ! bare sensible heat [w/m2]     [+ to atm]
 real,dimension(:),allocatable :: evb ! bare evaporation heat [w/m2]  [+ to atm]
 real,dimension(:),allocatable :: ghb ! bare ground heat flux [w/m2] [+ to soil]
 real,dimension(:),allocatable :: tr ! transpiration [w/m2]  [+ to atm]
 real,dimension(:),allocatable :: evc ! canopy evaporation heat [w/m2]  [+ to atm]
 real,dimension(:),allocatable :: chleaf ! leaf exchange coefficient
 real,dimension(:),allocatable :: chuc ! under canopy exchange coefficient
 real,dimension(:),allocatable :: chv2 ! veg 2m exchange coefficient
 real,dimension(:),allocatable :: chb2 ! bare 2m exchange coefficient 
 real,dimension(:),allocatable :: cosz ! cosine zenith angle
 real,dimension(:),allocatable :: lat  ! latitude [rad]
 real,dimension(:),allocatable :: lon  ! longitude [rad]
 real,dimension(:),allocatable :: fveg ! vegetation fraction
 real,dimension(:),allocatable :: fvgmax ! annual max vegetation
 real,dimension(:),allocatable :: fpice ! snow fraction of precip
 real,dimension(:),allocatable :: fcev ! canopy evaporation heat (w/m2) [+ to atm]
 real,dimension(:),allocatable :: fgev ! ground evaporation heat (w/m2) [+ to atm]
 real,dimension(:),allocatable :: fctr ! transpiration heat flux (w/m2) [+ to atm]
 real,dimension(:),allocatable :: qsnbot ! snowmelt out bottom of pack [mm/s]
 real,dimension(:),allocatable :: ponding ! snowmelt with no pack [mm]
 real,dimension(:),allocatable :: ponding1 ! snowmelt with no pack [mm]
 real,dimension(:),allocatable :: ponding2 ! snowmelt with no pack [mm]
 real,dimension(:),allocatable :: fsr ! total reflected solar radiation (w/m2)
 real,dimension(:),allocatable :: co2pp ! CO2 partial pressure [Pa]
 real,dimension(:),allocatable :: o2pp ! O2 partial pressure [Pa]
 real,dimension(:),allocatable :: foln ! nitrogen saturation [%
 real,dimension(:),allocatable :: tbot ! deep soil temperature [K]
 integer,dimension(:),allocatable :: isurban ! if cell is urban
 integer,dimension(:),allocatable :: slopetyp
 integer,dimension(:),allocatable :: soiltyp ! soil type
 integer,dimension(:),allocatable :: vegtyp ! vegetation type
 integer,dimension(:),allocatable :: ice ! glacier cell
 integer,dimension(:),allocatable :: isc ! soil color index
 integer,dimension(:),allocatable :: ist ! surface type 1-soil; 1-lake
 real,dimension(:),allocatable :: smcmax ! saturated soil moisture
 real,dimension(:),allocatable :: smcdry ! residual soil moisture
 real,dimension(:),allocatable :: smcref ! field capacity soil moisture
 real,dimension(:),allocatable :: errwat ! water balance error
 real,dimension(:),allocatable :: si0 ! initial soil moisture deficit
 real,dimension(:),allocatable :: si1 ! final soil moisture deficit
 real,dimension(:),allocatable :: zwt0 ! initial water table depth
 real,dimension(:),allocatable :: minzwt ! minimum water table depth

 !Multi layer
 real,dimension(:,:),allocatable :: stc ! snow/soil temperatures
 real,dimension(:,:),allocatable :: sh2o ! vol. soil liquid water (m3/m3)
 real,dimension(:,:),allocatable :: smc ! vol. soil moisture (m3/m3)
 real,dimension(:,:),allocatable :: smceq ! eq vol. soil moisture (m3/m3)
 real,dimension(:,:),allocatable :: zsnso ! snow layer depth (m)
 real,dimension(:,:),allocatable :: snice ! snow layer ice (mm)
 real,dimension(:,:),allocatable :: snliq ! snow layer liquid water (mm)
 real,dimension(:,:),allocatable :: ficeold ! snow layer ice fraction 
 real,dimension(:,:),allocatable :: zsoil ! depth to soil interfaces
 real,dimension(:,:),allocatable :: sldpth ! soil layer thickness

contains

 subroutine finalize()

  !Single layer
  deallocate(isnow) ! actual number of snow layers
  deallocate(z_ml) ! model height [m]
  deallocate(lwdn)       ! Downward longwave radiation flux at surface [W/m2]
  deallocate(swdn)      ! Downward shortwave radiation flux at surface [W/m2]
  deallocate(p_ml)    !  Surface pressure defined at intermediate level [Pa]
  deallocate(psfc)    ! Surface pressure [Pa]
  deallocate(prcp)       ! Precipitation rate (kg m-2 s-1)
  deallocate(t_ml)     ! Air temperature (K) [Forcing]
  deallocate(q_ml)         ! Surface specific humidity (kg kg-1)
  deallocate(u_ml)       ! West-to-east component of the surface [m/s]
  deallocate(v_ml)       ! North-to-south component of the surface [m/s]
  deallocate(fsh) ! total sensible heat (w/m2)
  deallocate(ssoil) ! soil heat (w/m2)
  deallocate(salb) ! surface albedo
  deallocate(fsno) ! snow cover fraction
  deallocate(swe) ! snow water equivalent (mm)
  deallocate(sndpth) ! snow depth (m)
  deallocate(emissi) ! net surface emissivity
  deallocate(qsfc1d) ! bulk surface specific humidity (kg/kg)
  deallocate(tv) ! vegetation canopy temperature
  deallocate(tg) ! ground surface temperature
  deallocate(canice) ! canopy-intercepted ice (mm)
  deallocate(canliq) ! canopy-intercepted liquid water (mm)
  deallocate(eah) ! canopy air vapor pressure (pa)
  deallocate(tah) ! canopy air temperature (K)
  deallocate(cm) ! momentum drag coefficient
  deallocate(ch) ! sensible heat exchange coefficient
  deallocate(fwet) ! weted or snowed fraction of the canopy
  deallocate(sneqvo) ! snow mass at last time step (mm h2o)
  deallocate(albold) ! snow albedo at last time step
  deallocate(qsnow) ! snowfall on the ground (mm/s)
  deallocate(wslake) ! lake water storage (mm)
  deallocate(zwt) ! water table depth (m)
  deallocate(dzwt) ! change in water table depth (m)
  deallocate(wa) ! water in the aquifer (mm)
  deallocate(wt) ! groundwater storage (mm)
  deallocate(smcwtd) ! soil moisture content in the transmission zone
  deallocate(deeprech) ! recharge to the water table when deep
  deallocate(rech) ! recharge to the water table (diagnostic)
  deallocate(lfmass) ! leaf mass (g/m2)
  deallocate(rtmass) ! mass of fine roots (g/m2)
  deallocate(stmass) ! stem mass (g/m2)
  deallocate(wood) ! mass of wood (incl. woody roots) (g/m2)
  deallocate(stblcp) ! stable carbon in deep soil (g/m2)
  deallocate(fastcp) ! short lived carbon, shallow soil (g/m2)
  deallocate(plai) ! leaf area index
  deallocate(psai) ! stem area index
  deallocate(tauss) ! non-dimensional snow age
  deallocate(t2mv) ! 2m temperature of vegetation part
  deallocate(t2mb) ! 2m temperature of bare ground part
  deallocate(q2mv) ! 2m mixing ratio of vegetation part
  deallocate(q2mb) ! 2m mixing ratio of bare ground part
  deallocate(trad) ! surface radiative temperature (k)
  deallocate(nee) ! net ecosys exchange (g/m2/s CO2)
  deallocate(gpp) ! gross primary assimilation [g/m2/s C]
  deallocate(npp) ! net primary productivity [g/m2/s C]
  deallocate(fvegmp) ! greenness vegetation fraction [-]
  deallocate(runsf) ! surface runoff [mm/s]
  deallocate(runsb) ! subsurface runoff [mm/s]
  deallocate(ecan) ! evaporation of intercepted water (mm/s)
  deallocate(etran) ! transpiration rate (mm/s)
  deallocate(esoil) ! soil surface evaporation rate (mm/s]
  deallocate(fsa) ! total absorbed solar radiation (w/m2)
  deallocate(fira) ! total net longwave rad (w/m2) [+ to atm]
  deallocate(apar) ! photosyn active energy by canopy (w/m2)
  deallocate(psn) ! total photosynthesis (umol co2/m2/s) [+]
  deallocate(sav) ! solar rad absorbed by veg. (w/m2)
  deallocate(sag) ! solar rad absorbed by ground (w/m2)
  deallocate(rssun) ! sunlit leaf stomatal resistance (s/m)
  deallocate(rssha) ! shaded leaf stomatal resistance (s/m)
  deallocate(bgap) ! between gap fraction
  deallocate(wgap) ! within gap fraction
  deallocate(tgv) ! under canopy ground temperature [K]
  deallocate(tgb) ! bare ground temperature [K]
  deallocate(chv) ! sensible heat exchange coefficient vegetated
  deallocate(chb) ! sensible heat exchange coefficient bare-ground
  deallocate(irc) ! canopy net LW rad. [w/m2] [+ to atm]
  deallocate(irg) ! veg ground net LW rad. [w/m2] [+ to atm]
  deallocate(shc) ! canopy sen. heat [w/m2]   [+ to atm]
  deallocate(shg) ! veg ground sen. heat [w/m2]   [+ to atm]
  deallocate(evg) ! veg ground evap. heat [w/m2]  [+ to atm]
  deallocate(ghv) ! veg ground heat flux [w/m2]  [+ to soil]
  deallocate(irb) ! bare net longwave rad. [w/m2] [+ to atm]
  deallocate(shb) ! bare sensible heat [w/m2]     [+ to atm]
  deallocate(evb) ! bare evaporation heat [w/m2]  [+ to atm]
  deallocate(ghb) ! bare ground heat flux [w/m2] [+ to soil]
  deallocate(tr) ! transpiration [w/m2]  [+ to atm]
  deallocate(evc) ! canopy evaporation heat [w/m2]  [+ to atm]
  deallocate(chleaf) ! leaf exchange coefficient
  deallocate(chuc) ! under canopy exchange coefficient
  deallocate(chv2) ! veg 2m exchange coefficient
  deallocate(chb2) ! bare 2m exchange coefficient 
  deallocate(cosz) ! cosine zenith angle
  deallocate(lat)  ! latitude [rad]
  deallocate(lon)  ! longitude [rad]
  deallocate(fveg) ! vegetation fraction
  deallocate(fvgmax) ! annual max vegetation
  deallocate(fpice) ! snow fraction of precip
  deallocate(fcev) ! canopy evaporation heat (w/m2) [+ to atm]
  deallocate(fgev) ! ground evaporation heat (w/m2) [+ to atm]
  deallocate(fctr) ! transpiration heat flux (w/m2) [+ to atm]
  deallocate(qsnbot) ! snowmelt out bottom of pack [mm/s]
  deallocate(ponding) ! snowmelt with no pack [mm]
  deallocate(ponding1) ! snowmelt with no pack [mm]
  deallocate(ponding2) ! snowmelt with no pack [mm]
  deallocate(fsr) ! total reflected solar radiation (w/m2)
  deallocate(co2pp) ! CO2 partial pressure [Pa]
  deallocate(o2pp) ! O2 partial pressure [Pa]
  deallocate(foln) ! nitrogen saturation [%
  deallocate(tbot) ! deep soil temperature [K]
  deallocate(isurban) ! if cell is urban
  deallocate(slopetyp)
  deallocate(soiltyp) ! soil type
  deallocate(vegtyp) ! vegetation type
  deallocate(ice) ! glacier cell
  deallocate(isc) ! soil color index
  deallocate(ist) ! surface type
  deallocate(smcmax) ! saturated soil moisture
  deallocate(smcdry) ! residual soil moisture
  deallocate(smcref)
  deallocate(errwat) 
  deallocate(si0) 
  deallocate(si1) 
  deallocate(zwt0)
  deallocate(minzwt)

  !Multi layer
  deallocate(stc)
  deallocate(zsoil)
  deallocate(sh2o)
  deallocate(smc)
  deallocate(smceq) ! eq vol. soil moisture (m3/m3)
  deallocate(zsnso) ! snow layer depth (m)
  deallocate(snice)! snow layer ice (mm)
  deallocate(snliq) ! snow layer liquid water (mm)
  deallocate(ficeold) ! snow layer ice fraction 
  deallocate(sldpth) !soil layer thickness

 end subroutine finalize

 subroutine initialize()
 
  USE module_sf_noahmplsm,only : read_mp_veg_parameters,noahmp_options
  implicit none
  integer :: i
  !Allocate all memory

  !Single layer
  allocate(isnow(ncells)) ! actual number of snow layers
  allocate(z_ml(ncells)) ! model height [m]
  allocate(lwdn(ncells))       ! Downward longwave radiation flux at surface [W/m2]
  allocate(swdn(ncells))      ! Downward shortwave radiation flux at surface [W/m2]
  allocate(p_ml(ncells))    !  Surface pressure defined at intermediate level [Pa]
  allocate(psfc(ncells))    ! Surface pressure [Pa]
  allocate(prcp(ncells))       ! Precipitation rate (kg m-2 s-1)
  allocate(t_ml(ncells))     ! Air temperature (K) [Forcing]
  allocate(q_ml(ncells))         ! Surface specific humidity (kg kg-1)
  allocate(u_ml(ncells))       ! West-to-east component of the surface [m/s]
  allocate(v_ml(ncells))       ! North-to-south component of the surface [m/s]
  allocate(fsh(ncells)) ! total sensible heat (w/m2)
  allocate(ssoil(ncells)) ! soil heat (w/m2)
  allocate(salb(ncells)) ! surface albedo
  allocate(fsno(ncells)) ! snow cover fraction
  allocate(swe(ncells)) ! snow water equivalent (mm)
  allocate(sndpth(ncells)) ! snow depth (m)
  allocate(emissi(ncells)) ! net surface emissivity
  allocate(qsfc1d(ncells)) ! bulk surface specific humidity (kg/kg)
  allocate(tv(ncells)) ! vegetation canopy temperature
  allocate(tg(ncells)) ! ground surface temperature
  allocate(canice(ncells)) ! canopy-intercepted ice (mm)
  allocate(canliq(ncells)) ! canopy-intercepted liquid water (mm)
  allocate(eah(ncells)) ! canopy air vapor pressure (pa)
  allocate(tah(ncells)) ! canopy air temperature (K)
  allocate(cm(ncells)) ! momentum drag coefficient
  allocate(ch(ncells)) ! sensible heat exchange coefficient
  allocate(fwet(ncells)) ! weted or snowed fraction of the canopy
  allocate(sneqvo(ncells)) ! snow mass at last time step (mm h2o)
  allocate(albold(ncells)) ! snow albedo at last time step
  allocate(qsnow(ncells)) ! snowfall on the ground (mm/s)
  allocate(wslake(ncells)) ! lake water storage (mm)
  allocate(zwt(ncells)) ! water table depth (m)
  allocate(dzwt(ncells)) ! change in water table depth (m)
  allocate(wa(ncells)) ! water in the aquifer (mm)
  allocate(wt(ncells)) ! groundwater storage (mm)
  allocate(smcwtd(ncells)) ! soil moisture content in the transmission zone
  allocate(deeprech(ncells)) ! recharge to the water table when deep
  allocate(rech(ncells)) ! recharge to the water table (diagnostic)
  allocate(lfmass(ncells)) ! leaf mass (g/m2)
  allocate(rtmass(ncells)) ! mass of fine roots (g/m2)
  allocate(stmass(ncells)) ! stem mass (g/m2)
  allocate(wood(ncells)) ! mass of wood (incl. woody roots) (g/m2)
  allocate(stblcp(ncells)) ! stable carbon in deep soil (g/m2)
  allocate(fastcp(ncells)) ! short lived carbon, shallow soil (g/m2)
  allocate(plai(ncells)) ! leaf area index
  allocate(psai(ncells)) ! stem area index
  allocate(tauss(ncells)) ! non-dimensional snow age
  allocate(t2mv(ncells)) ! 2m temperature of vegetation part
  allocate(t2mb(ncells)) ! 2m temperature of bare ground part
  allocate(q2mv(ncells)) ! 2m mixing ratio of vegetation part
  allocate(q2mb(ncells)) ! 2m mixing ratio of bare ground part
  allocate(trad(ncells)) ! surface radiative temperature (k)
  allocate(nee(ncells)) ! net ecosys exchange (g/m2/s CO2)
  allocate(gpp(ncells)) ! gross primary assimilation [g/m2/s C]
  allocate(npp(ncells)) ! net primary productivity [g/m2/s C]
  allocate(fvegmp(ncells)) ! greenness vegetation fraction [-]
  allocate(runsf(ncells)) ! surface runoff [mm/s]
  allocate(runsb(ncells)) ! subsurface runoff [mm/s]
  allocate(ecan(ncells)) ! evaporation of intercepted water (mm/s)
  allocate(etran(ncells)) ! transpiration rate (mm/s)
  allocate(esoil(ncells)) ! soil surface evaporation rate (mm/s]
  allocate(fsa(ncells)) ! total absorbed solar radiation (w/m2)
  allocate(fira(ncells)) ! total net longwave rad (w/m2) [+ to atm]
  allocate(apar(ncells)) ! photosyn active energy by canopy (w/m2)
  allocate(psn(ncells)) ! total photosynthesis (umol co2/m2/s) [+]
  allocate(sav(ncells)) ! solar rad absorbed by veg. (w/m2)
  allocate(sag(ncells)) ! solar rad absorbed by ground (w/m2)
  allocate(rssun(ncells)) ! sunlit leaf stomatal resistance (s/m)
  allocate(rssha(ncells)) ! shaded leaf stomatal resistance (s/m)
  allocate(bgap(ncells)) ! between gap fraction
  allocate(wgap(ncells)) ! within gap fraction
  allocate(tgv(ncells)) ! under canopy ground temperature [K]
  allocate(tgb(ncells)) ! bare ground temperature [K]
  allocate(chv(ncells)) ! sensible heat exchange coefficient vegetated
  allocate(chb(ncells)) ! sensible heat exchange coefficient bare-ground
  allocate(irc(ncells)) ! canopy net LW rad. [w/m2] [+ to atm]
  allocate(irg(ncells)) ! veg ground net LW rad. [w/m2] [+ to atm]
  allocate(shc(ncells)) ! canopy sen. heat [w/m2]   [+ to atm]
  allocate(shg(ncells)) ! veg ground sen. heat [w/m2]   [+ to atm]
  allocate(evg(ncells)) ! veg ground evap. heat [w/m2]  [+ to atm]
  allocate(ghv(ncells)) ! veg ground heat flux [w/m2]  [+ to soil]
  allocate(irb(ncells)) ! bare net longwave rad. [w/m2] [+ to atm]
  allocate(shb(ncells)) ! bare sensible heat [w/m2]     [+ to atm]
  allocate(evb(ncells)) ! bare evaporation heat [w/m2]  [+ to atm]
  allocate(ghb(ncells)) ! bare ground heat flux [w/m2] [+ to soil]
  allocate(tr(ncells)) ! transpiration [w/m2]  [+ to atm]
  allocate(evc(ncells)) ! canopy evaporation heat [w/m2]  [+ to atm]
  allocate(chleaf(ncells)) ! leaf exchange coefficient
  allocate(chuc(ncells)) ! under canopy exchange coefficient
  allocate(chv2(ncells)) ! veg 2m exchange coefficient
  allocate(chb2(ncells)) ! bare 2m exchange coefficient 
  allocate(cosz(ncells)) ! cosine zenith angle
  allocate(lat(ncells))  ! latitude [rad]
  allocate(lon(ncells))  ! longitude [rad]
  allocate(fveg(ncells)) ! vegetation fraction
  allocate(fvgmax(ncells)) ! annual max vegetation
  allocate(fpice(ncells)) ! snow fraction of precip
  allocate(fcev(ncells)) ! canopy evaporation heat (w/m2) [+ to atm]
  allocate(fgev(ncells)) ! ground evaporation heat (w/m2) [+ to atm]
  allocate(fctr(ncells)) ! transpiration heat flux (w/m2) [+ to atm]
  allocate(qsnbot(ncells)) ! snowmelt out bottom of pack [mm/s]
  allocate(ponding(ncells)) ! snowmelt with no pack [mm]
  allocate(ponding1(ncells)) ! snowmelt with no pack [mm]
  allocate(ponding2(ncells)) ! snowmelt with no pack [mm]
  allocate(fsr(ncells)) ! total reflected solar radiation (w/m2)
  allocate(co2pp(ncells)) ! CO2 partial pressure [Pa]
  allocate(o2pp(ncells)) ! O2 partial pressure [Pa]
  allocate(foln(ncells)) ! nitrogen saturation [%
  allocate(tbot(ncells)) ! deep soil temperature [K]
  allocate(isurban(ncells)) ! if cell is urban
  allocate(slopetyp(ncells))
  allocate(soiltyp(ncells)) ! soil type
  allocate(vegtyp(ncells)) ! vegetation type
  allocate(ice(ncells)) ! glacier cell
  allocate(isc(ncells)) ! soil color index
  allocate(ist(ncells)) ! surface type
  allocate(smcmax(ncells)) ! saturated soil moisture
  allocate(smcdry(ncells)) ! residual soil moisture
  allocate(smcref(ncells))
  allocate(errwat(ncells))
  allocate(si0(ncells))
  allocate(si1(ncells))
  allocate(zwt0(ncells))
  allocate(minzwt(ncells))

  !Multi layer
  allocate(stc(ncells,-nsnow+1:nsoil))
  allocate(zsoil(ncells,nsoil))
  allocate(sh2o(ncells,nsoil))
  allocate(smc(ncells,nsoil))
  allocate(smceq(ncells,nsoil)) ! eq vol. soil moisture (m3/m3)
  allocate(zsnso(ncells,-nsnow+1:nsoil)) ! snow layer depth (m)
  allocate(snice(ncells,-nsnow+1:0))! snow layer ice (mm)
  allocate(snliq(ncells,-nsnow+1:0)) ! snow layer liquid water (mm)
  allocate(ficeold(ncells,-nsnow+1:0)) ! snow layer ice fraction 
  allocate(sldpth(ncells,nsoil)) !soil layer thickness
  
  !Set initial values
  
  !Single layer
  z_ml = badval
  lwdn = badval
  swdn = badval
  p_ml = badval
  psfc = badval
  prcp = badval
  t_ml = badval
  q_ml = badval
  u_ml = badval
  v_ml = badval
  fsh = badval
  ssoil = badval
  salb = badval 
  fsno = badval
  sndpth = badval
  emissi = badval
  qsfc1d = badval
  tv = badval
  tg = badval
  canice = badval
  canliq = badval
  eah = badval
  tah = badval
  cm = badval
  ch = badval
  fwet = badval
  sneqvo = badval
  albold = badval
  qsnow = badval
  wslake = badval
  zwt = badval
  wa = badval
  wt = badval
  smcwtd = badval
  deeprech = badval
  rech = badval
  lfmass = badval
  stmass = badval
  wood = badval
  stblcp = badval
  fastcp = badval
  plai = badval
  psai = badval
  tauss = badval
  t2mv = badval
  t2mb = badval
  q2mv = badval
  q2mb = badval
  trad = badval
  nee = badval
  gpp = badval
  npp = badval
  fvegmp = badval
  etran = badval
  esoil = badval
  fsa = badval
  fira = badval
  apar = badval
  psn = badval
  sav = badval
  sag = badval
  rssun = badval
  rssha = badval
  bgap = badval
  wgap = badval
  tgv = badval
  tgb = badval
  chv = badval
  chb = badval
  irc = badval
  irg = badval
  shc = badval
  shg = badval
  evg = badval
  ghv = badval
  irb = badval
  shb = badval
  evb = badval
  ghb = badval
  tr = badval
  evc = badval
  chleaf = badval
  chuc = badval
  chv2 = badval
  chb2 = badval
  cosz = badval
  lat = badval
  fveg = badval
  fvgmax = badval
  fpice = badval
  fcev = badval
  fgev = badval
  fctr = badval
  qsnbot = badval
  ponding = badval
  ponding1 = badval
  ponding2 = badval
  fsr = badval
  o2pp = badval
  foln = badval
  tbot= badval
  isurban = -9999
  slopetyp = -9999
  soiltyp = -9999
  vegtyp = -9999
  isnow = -9999
  ice = -9999
  isc = -9999
  ist = -9999
  si0 = badval
  si1 = badval
  zwt0 = badval
  minzwt = badval

  !Multi layer
  zsoil = badval
  sh2o = badval
  smc = badval
  smceq = badval
  zsnso = badval
  snice = badval
  snliq = badval
  ficeold = badval

  ! Read our lookup tables and parameter tables:  VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
  call soil_veg_gen_parm(llanduse,lsoil,vegparm_file,soilparm_file,genparm_file)

  ! Read the Noah-MP table
  call read_mp_veg_parameters(llanduse,mptable_file)

  ! Define the NOAH-MP options
  CALL NOAHMP_OPTIONS(IDVEG  ,IOPT_CRS  ,IOPT_BTR  ,IOPT_RUN,IOPT_SFC,IOPT_FRZ , &
                     IOPT_INF  ,IOPT_RAD  ,IOPT_ALB  ,IOPT_SNF,IOPT_TBOT,IOPT_STC )

 end subroutine

 subroutine run_model(ncores)

  implicit none
  integer :: i,ncores,omp_get_thread_num
  real*8 :: ostart,oend,ostart1,oend1,omp_get_wtime,tmp
  !real :: stc0(-nsnow+1:nsoil),sh2o0(1:nsoil),smc0(1:nsoil),sldpth0(1:nsoil),zsoil0(1:nsoil)
  !real :: ficeold0(-nsnow+1:0),zsnso0(-nsnow+1:nsoil),snice0(-nsnow+1:0),snliq0(-nsnow+1:0),smceq0(1:nsoil)
  !real :: fveg0,fvgmax0,t_ml0,p_ml0,psfc0,u_ml0,v_ml0
  !real :: q_ml0,swdn0,lwdn0,prcp0,tbot0,co2pp0
  !real :: ZWT0     , WA0      , WT0      , WSLAKE0  , LFMASS0  , RTMASS0
  !real :: STMASS0  , WOOD0    , STBLCP0  , FASTCP0  , PLAI0    , PSAI0
  !real :: CM0      , CH0      , TAUSS0               
  !real :: SMCWTD0  ,DEEPRECH0 , RECH0                               
  !real :: FSA0     , FSR0     , FIRA0    , FSH0     , SSOIL0   , FCEV0  
  !real :: FGEV0    , FCTR0    , ECAN0    , ETRAN0   , ESOIL0   , TRAD0  
  !real :: TGB0     , TGV0     , T2MV0    , T2MB0    , Q2MV0    , Q2MB0 
  !real :: RUNSF0   , RUNSB0   , APAR0    , PSN0    , SAV0     , SAG0   
  !real :: FSNO0    , NEE0     , GPP0     , NPP0     , FVEGMP0  , SALB0  
  !real :: QSNBOT0  , PONDING0 , PONDING10, PONDING20, RSSUN0   , RSSHA0  
  !real :: BGAP0    , WGAP0    , CHV0     , CHB0     , EMISSI0 
  !real :: SHG0     , SHC0     , SHB0     , EVG0     , EVB0     , GHV0  
  !real :: GHB0     , IRG0     , IRC0     , IRB0     , TR0      , EVC0    
  !real :: CHLEAF0  , CHUC0    , CHV20    , CHB20    , FPICE0
  !real :: LAT0     , COSZ0
  !integer :: VEGTYP0  , ISURBAN0 , ICE0 , IST0, ISC0, SOILTYP0, SLOPETYP0
  !real :: O2PP0    , FOLN0    , Z_ML0, ALBOLD0  , SNEQVO0
  !real :: TAH0     , EAH0     , FWET0
  !real :: CANLIQ0  , CANICE0  , TV0      , TG0      , QSFC1D0  , QSNOW0 
  !real :: ISNOW0   , SNDPTH0  , SWE0
  !real :: LON0,DZWT0

  !real :: lat,yearlen,julian,cosz
  !integer :: nsoil
  call OMP_SET_NUM_THREADS(ncores)
  ostart = omp_get_wtime()
  !$OMP PARALLEL
  !$OMP DO
  !!$OMP DO PRIVATE(BEXP,SMCDRY,F1,SMCMAX,SMCREF,PSISAT,DKSAT,&
  !!$OMP DWSAT,SMCWLT,QUARTZ,SLOPE,CSOIL,ZBOT,CZIL,KDT,FRZX,NROOT,RGL,RSMIN,HS,&
  !!$OMP RSMAX,TOPT)
  !!$OMP ISWATER,ISBARREN,ISSNOW,EBLFOREST,CH2OP,DLEAF,Z0MVT,HVT,HVB,DEN,RC,SAIM,LAIM,SLA,DILEFC,DILEFW,FRAGR,LTOVRC,C3PSN,&
  !!$OMP KC25,AKC,KO25,AKO,VCMX25,AVCMX,BP,MP,QE25,AQE,RMF25,RMS25,RMR25,ARM,FOLNMX,TMIN,XL,RHOL,RHOS,TAUL,TAUS,RMP,CWPVT,WRRAT,&
  !!$OMP WDPPOL,TDLEF,IK,IM,TMP10,TMP11,TMP12,TMP13,TMP14,TMP15,TMP16,slarea,eps)

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
            zwt0(i),minzwt(i)) ! OTHER
   !ostart1 = omp_get_wtime()
   !Copy cell variables
   !stc0 = STC(i,:); sh2o0 = SH2O(i,:); smc0 = SMC(i,:); zsnso0 = zsnso(i,:); snice0 = snice(i,:); snliq0 = snliq(i,:)
   !sldpth0 = sldpth(i,:); zsoil0 = zsoil(i,:); ficeold0 = ficeold(i,:); smceq0 = smceq(i,:)
   !fveg0 = fveg(i); fvgmax0 = fvgmax(i); t_ml0 = t_ml(i); p_ml0 = p_ml(i); psfc0 = psfc(i); u_ml0 = u_ml(i); v_ml0 = v_ml(i)
   !q_ml0 = q_ml(i); swdn0 = swdn(i); lwdn0 = lwdn(i); prcp0 = prcp(i); tbot0 = tbot(i); co2pp0 = co2pp(i)
   !ZWT0 = zwt(i); WA0 = wa(i); WT0 = wt(i); WSLAKE0 = wslake(i); LFMASS0 = lfmass(i); RTMASS0 = rtmass(i)
   !STMASS0 = stmass(i); WOOD0 = wood(i); STBLCP0 = stblcp(i); FASTCP0 =fastcp(i); PLAI0 = plai(i); PSAI0 = psai(i)
   !CM0 = cm(i); CH0 =  ch(i); TAUSS0 = tauss(i)
   !SMCWTD0 = smcwtd(i); DEEPRECH0 = deeprech(i); RECH0 = rech(i)
   !FSA0 = fsa(i); FSR0 = fsr(i); FIRA0 = fira(i); FSH0 = fsh(i); SSOIL0 = ssoil(i); FCEV0 = fcev(i)
   !FGEV0 = fgev(i); FCTR0 = fctr(i); ECAN0 = ecan(i); ETRAN0 = etran(i); ESOIL0 = esoil(i); TRAD0 = trad(i)
   !TGB0 = tgb(i); TGV0 = tgv(i); T2MV0 = t2mv(i); T2MB0 = t2mb(i); Q2MV0 = q2mv(i); Q2MB0 = q2mb(i)
   !RUNSF0 = runsf(i); RUNSB0 = runsb(i); APAR0 = apar(i); PSN0 = psn(i); SAV0 = sav(i); SAG0 = sag(i)
   !FSNO0 = fsno(i); NEE0 = nee(i); GPP0 = gpp(i); NPP0 = npp(i); FVEGMP0 = fvegmp(i); SALB0 = salb(i)
   !QSNBOT0 = qsnbot(i); PONDING0 = ponding(i); PONDING10 = ponding1(i); PONDING20 = ponding2(i); RSSUN0 = rssun(i); RSSHA0 = rssha(i)
   !BGAP0 = bgap(i); WGAP0 = wgap(i); CHV0 = chv(i); CHB0 = chb(i); EMISSI0 = emissi(i)
   !SHG0 = shg(i); SHC0 = shc(I); SHB0 = shb(i); EVG0 = evg(i); EVB0 = evb(i); GHV0 = ghv(i)
   !GHB0 = ghb(i); IRG0 = irg(i); IRC0 = irc(i); IRB0 = irb(i); TR0 =  tr(i); EVC0 = evc(i)
   !CHLEAF0 = chleaf(i); CHUC0 = chuc(i); CHV20 = chv2(i); CHB20 = chb2(i); FPICE0 = fpice(i)
   !LAT0 = lat(i); COSZ0 = cosz(i)
   !VEGTYP0 = vegtyp(i); ISURBAN0 = isurban(i); ICE0 = ice(i); IST0 = ist(i); ISC0 = isc(i)
   !O2PP0 = o2pp(i); FOLN0 = foln(i); Z_ML0 = z_ml(i); ALBOLD0 = albold(i); SNEQVO0 = sneqvo(i)
   !TAH0 = tah(i); EAH0 = eah(i); FWET0 = fwet(i)
   !CANLIQ0 = canliq(i); CANICE0 = canice(i); TV0 = tv(i); TG0 = tg(i); QSFC1D0 = qsfc1d(i); QSNOW0 = qsnow(i)
   !ISNOW0 = isnow(i); SNDPTH0 =  sndpth(i); SWE0 = swe(i)
   !SOILTYP0 = soiltyp(i); SLOPETYP0 = slopetyp(i);LON0 = lon(i);DZWT0 = dzwt(i)

   !call run_model_cell(&
   !         LAT0     , YEARLEN , JULIAN  , COSZ0    , & ! IN : Time/Space-related
   !         DT      , DX      , NSOIL   , ZSOIL0   , NSNOW,  & ! IN : Model configuration 
   !         fveg0    , fvgmax0  , VEGTYP0  , ISURBAN0 , ICE0     , IST0     , & ! IN : Vegetation/Soil characteristics
   !         ISC0     , SMCEQ0   ,                                         & ! IN : Vegetation/Soil characteristics
   !         IZ0TLND ,                                                   & ! IN : User options
   !         t_ml0    , p_ml0    , psfc0    , u_ml0    , v_ml0    , q_ml0    , & ! IN : Forcing
   !         swdn0    , lwdn0    , prcp0    , tbot0    , co2pp0   , & ! IN : Forcing
   !         O2PP0    , FOLN0    , FICEOLD0 , Z_ML0    ,           & ! IN : Forcing
   !         ALBOLD0  , SNEQVO0  ,                                         & ! IN/OUT : 
   !         stc0    , sh2o0   , smc0     , TAH0     , EAH0     , FWET0    , & ! IN/OUT : 
   !         CANLIQ0  , CANICE0  , TV0      , TG0      , QSFC1D0  , QSNOW0   , & ! IN/OUT : 
   !         ISNOW0   , ZSNSO0   , SNDPTH0  , SWE0     , SNICE0   , SNLIQ0   , & ! IN/OUT : 
   !         ZWT0     , WA0      , WT0      , WSLAKE0  , LFMASS0  , RTMASS0  , & ! IN/OUT : 
   !         STMASS0  , WOOD0    , STBLCP0  , FASTCP0  , PLAI0    , PSAI0    , & ! IN/OUT : 
   !         CM0      , CH0      , TAUSS0   ,                               & ! IN/OUT : 
   !         SMCWTD0  ,DEEPRECH0 , RECH0    ,                               & ! IN/OUT :
   !         FSA0     , FSR0     , FIRA0    , FSH0     , SSOIL0   , FCEV0    , & ! OUT : 
   !         FGEV0    , FCTR0    , ECAN0    , ETRAN0   , ESOIL0   , TRAD0    , & ! OUT : 
   !         TGB0     , TGV0     , T2MV0    , T2MB0    , Q2MV0    , Q2MB0    , & ! OUT : 
   !         RUNSF0   , RUNSB0   , APAR0    , PSN0    , SAV0     , SAG0     , & ! OUT : 
   !         FSNO0    , NEE0     , GPP0     , NPP0     , FVEGMP0  , SALB0    , & ! OUT : 
   !         QSNBOT0  , PONDING0 , PONDING10, PONDING20, RSSUN0   , RSSHA0   , & ! OUT : 
   !         BGAP0    , WGAP0    , CHV0     , CHB0     , EMISSI0  ,           & ! OUT : 
   !         SHG0     , SHC0     , SHB0     , EVG0     , EVB0     , GHV0     , & ! OUT :
   ! 	    GHB0     , IRG0     , IRC0     , IRB0     , TR0      , EVC0     , & ! OUT :
   ! 	    CHLEAF0  , CHUC0    , CHV20    , CHB20    , FPICE0, &
   !         SOILTYP0, SLOPETYP0,LON0,NOWDATE,ITIME,SLDPTH0,DZWT0) ! OTHER
   !stc(i,:) = STC0; sh2o(i,:) = SH2O0; smc(i,:) = SMC0; zsnso(i,:) = zsnso0; snice(i,:) = snice0; snliq(i,:) = snliq0
   !sldpth(i,:) = sldpth0; zsoil(i,:) = zsoil0; ficeold(i,:) = ficeold0; smceq(i,:) = smceq0
   !STC(i,:) = stc0; SH2O(i,:) = sh2o0; SMC(i,:) = smc0
   !fveg(i) = fveg0; fvgmax(i) = fvgmax0; t_ml(i) = t_ml0; p_ml(i) = p_ml0; psfc(i) = psfc0; u_ml(i) = u_ml0; v_ml(i) = v_ml0
   !q_ml(i) = q_ml0; swdn(i) = swdn0; lwdn(i) = lwdn0; prcp(i) = prcp0; tbot(i) = tbot0; co2pp(i) = co2pp0
   !oend1 = omp_get_wtime()
   !!print*,'Run hsu',oend1-ostart1
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  oend = omp_get_wtime()
  !print*,'Run all hsus',oend - ostart,ncores

 end subroutine

end module model

subroutine run_model_cell(&
            LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
            DT      , DX      , NSOIL   , ZSOIL   , NSNOW,  & ! IN : Model configuration 
            FVEG    , FVGMAX  , VEGTYP  , ISURBAN , ICE     , IST     , & ! IN : Vegetation/Soil characteristics
            ISC     , SMCEQ   ,                                         & ! IN : Vegetation/Soil characteristics
            IZ0TLND ,                                                   & ! IN : User options
            T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
            SWDN    , LWDN    , PRCP    , TBOT    , CO2PP   , & ! IN : Forcing
            O2PP    , FOLN    , FICEOLD , Z_ML    ,           & ! IN : Forcing
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
            SOILTYP, SLOPETYP,LON,NOWDATE,ITIME,SLDPTH,DZWT,ERRWAT,si0,si1,zwt0,minzwt) ! OTHER

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
    REAL :: QSPRING,si0,si1,zwt0,minzwt
    REAL, DIMENSION(1:NSOIL) :: SLDPTH
    REAL :: DZWT
    REAL :: ERRWAT

    !Calculate cosz
    !ostart = omp_get_wtime()
    CALL CALC_DECLIN (NOWDATE, LAT, LON, COSZ)
    !oend = omp_get_wtime()
    !print*,'DECLIN',oend - ostart

    !Set parameters
    !ostart = omp_get_wtime()
    CALL REDPRM (VEGTYP,SOILTYP,SLOPETYP,ZSOIL,NSOIL,ISURBAN)
    !oend = omp_get_wtime()
    !print*,'REDPRM',oend - ostart

    if (itime .eq. 0) then 
      !Define the equilibrium soil water content
      CALL EQSMOISTURE(NSOIL ,  ZSOIL , SMCMAX , SMCWLT ,DWSAT, DKSAT  ,BEXP  , SMCEQ)
      !Make sure that below the water table the layers are saturated and initialize the deep soil moisture
      !CALL INITIALIZE_GROUNDWATER(NSOIL,DT,ZSOIL,SLDPTH,ZWT,SMC,SMH2O,SMCEQ,SMCWTD,&
      !     DEEPRECH,RECH,QSPRING,DZWT)
    end if

    !Update water table depth
    !print*,'BEFORE',ZWT,DZWT,QSPRING,DEEPRECH
    
    !DZWT = 0.0
    !if (IOPT_RUN .eq. 5)then
    !PROBLEM!
    !Calculate the soil moisture deficit
    !call Calculate_Deficit(si0,smcmax,smcwtd,zwt,sldpth,smc,nsoil)

    !Update the water table depth
    call UPDATEWTD  (NSOIL,SLDPTH,SMCEQ,SMCMAX,SMCWLT,-PSISAT,BEXP,I,J,DZWT,ZWT ,SMC, SMH2O ,SMCWTD,QSPRING,&
                     DT,DEEPRECH,DKSAT)
    DEEPRECH = 0.0
    !Calculate the soil moisture deficit
    zwt0 = zwt
    call Calculate_Deficit(si0,smcmax,smcwtd,zwt,sldpth,smc,nsoil)
    !if (itime .eq. 9)return
    !endif
    !print*,'AFTER',ZWT,DZWT,QSPRING,DEEPRECH

    !Run the model
     !write(*,'("Before call to NOAHMP_SFLX")')
     !write(*,'(10x, "ICE        = ",  I10   )') ICE
     !write(*,'(10x, "IST        = ",  I10   )') IST
     !write(*,'(10x, "VEGTYPE    = ",  I10   )') VEGTYP
     !write(*,'(10x, "ISC        = ",  I10   )') ISC
     !write(*,'(10x, "NSOIL      = ",  I10   )') NSOIL
     !write(*,'(10x, "ZSOIL      = ", 7F20.10)') ZSOIL
     !write(*,'(10x, "DT         = ",  F20.10)') DT
     !write(*,'(10x, "Q2         = ",  F20.10)') Q_ML
     !write(*,'(10x, "SFCTMP     = ",  F20.10)') T_ML
     !write(*,'(10x, "UU         = ",  F20.10)') U_ML
     !write(*,'(10x, "VV         = ",  F20.10)') V_ML
     !write(*,'(10x, "SOLDN      = ",  F20.10)') SWDN
     !write(*,'(10x, "LWDN       = ",  F20.10)') LWDN
     !write(*,'(10x, "PRCP       = ",  F20.10)') PRCP
     !write(*,'(10x, "ZLVL       = ",  F20.10)') Z_ML
     !write(*,'(10x, "CO2AIR     = ",  F20.10)') CO2PP
     !write(*,'(10x, "O2AIR      = ",  F20.10)') O2PP
     !write(*,'(10x, "COSZ       = ",  F20.10)') COSZ
     !write(*,'(10x, "TBOT       = ",  F20.10)') TBOT
     !write(*,'(10x, "FOLN       = ",  F20.10)') FOLN
     !write(*,'(10x, "SFCPRS     = ",  F20.10)') SFCPRS
     !write(*,'(10x, "SHDFAC     = ",  F20.10)') FVEG
     !write(*,'(10x, "LAT        = ",  F20.10)') LAT
     !write(*,'(10x, "DZ8W       = ",  F20.10)') DZ8W1D
     !write(*,'(10x, "EAH        = ",  F20.10)') EAH
     !write(*,'(10x, "TAH        = ",  F20.10)') TAH
     !write(*,'(10x, "FWET       = ",  F20.10)') FWET
     !write(*,'(10x, "FICEOLD    = ", 7F20.10)') FICEOLD
     !write(*,'(10x, "QSNOW      = ",  F20.10)') QSNOW
     !write(*,'(10x, "SNEQVO     = ",  F20.10)') SNEQVO
     !write(*,'(10x, "ISNOW      = ",  I10   )') ISNOW
     !write(*,'(10x, "ZSNSO      = ", 7F20.10)') ZSNSO
     !write(*,'(10x, "CANLIQ     = ",  F20.10)') CANLIQ
     !write(*,'(10x, "CANICE     = ",  F20.10)') CANICE
     !write(*,'(10x, "SNOWH      = ",  F20.10)') SNOWH
     !write(*,'(10x, "SNEQV      = ",  F20.10)') SNEQV
     !write(*,'(10x, "SNICE      = ", 7F20.10)') SNICE
     !write(*,'(10x, "SNLIQ      = ", 7F20.10)') SNLIQ
     !write(*,'(10x, "TV         = ",  F20.10)') TV
     !write(*,'(10x, "TG         = ",  F20.10)') TG
     !write(*,'(10x, "STC        = ", 7F20.10)') STC
     !write(*,'(10x, "SH2O       = ", 7F20.10)') SMH2O
     !write(*,'(10x, "SMC        = ", 7F20.10)') SMC
     !write(*,'(10x, "ZWT        = ",  F20.10)') ZWT
     !write(*,'(10x, "WA         = ",  F20.10)') WA
     !write(*,'(10x, "WT         = ",  F20.10)') WT
     !write(*,'(10x, "WSLAKE     = ",  F20.10)') WSLAKE
     !write(*,'(10x, "LFMASS     = ",  F20.10)') LFMASS
     !write(*,'(10x, "RTMASS     = ",  F20.10)') RTMASS
     !write(*,'(10x, "STMASS     = ",  F20.10)') STMASS
     !write(*,'(10x, "WOOD       = ",  F20.10)') WOOD
     !write(*,'(10x, "STBLCP     = ",  F20.10)') STBLCP
     !write(*,'(10x, "FASTCP     = ",  F20.10)') FASTCP
     !write(*,'(10x, "PLAI        = ",  F20.10)') PLAI
     !write(*,'(10x, "SAI        = ",  F20.10)') PSAI
     !write(*,'(10x, "ALBOLD     = ",  F20.10)') ALBOLD
     !write(*,'(10x, "CM         = ",  F20.10)') CM
     !write(*,'(10x, "CH         = ",  F20.10)') CH
     !write(*,'(10x, "DX         = ",  F20.10)') DX
     !write(*,'(10x, "ISURBAN    = ",  I10   )') ISURBAN
     !write(*,'(10x, "IZ0TLND    = ",  I10   )') IZ0TLND
     !write(*,'(10x, "QC         = ",  F20.10)') QC
     !write(*,'(10x, "PBLH       = ",  F20.10)') PBLH
     !write(*,'(10x, "QSFC       = ",  F20.10)') QSFC1D
     !write(*,'(10x, "PSFC       = ",  F20.10)') PSFC

    !ostart = omp_get_wtime()
    CALL NOAHMP_SFLX (&
            I       , J       , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
            DT      , DX      , DZ8W1D  , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration 
            FVEG    , FVGMAX  , VEGTYP  , ISURBAN , ICE     , IST     , & ! IN : Vegetation/Soil characteristics
            ISC     , SMCEQ   ,                                         & ! IN : Vegetation/Soil characteristics
            IZ0TLND ,                                                   & ! IN : User options
            T_ML    , P_ML    , PSFC    , U_ML    , V_ML    , Q_ML    , & ! IN : Forcing
            QC      , SWDN    , LWDN    , PRCP    , TBOT    , CO2PP   , & ! IN : Forcing
            O2PP    , FOLN    , FICEOLD , PBLH    , Z_ML    ,           & ! IN : Forcing
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
	    CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE, ERRWAT,minzwt)

    !oend = omp_get_wtime()
    !print*,'FLX',oend - ostart

     !Add in QSPRING
     !runsb = runsb + 1000*qspring/dt
     runsb = 1000*qspring/dt
     !runsf = runsf + 1000*qspring/dt
     !print*,runsb,qspring

    !Calculate the soil moisture deficit
    call Calculate_Deficit(si1,smcmax,smcwtd,zwt,sldpth,smc,nsoil)

end subroutine run_model_cell

!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM( MMINLU, MMINSL, VEGPARM_FILE, SOILPARM_FILE,GENPARM_FILE)
!-----------------------------------------------------------------

  USE module_sf_noahmplsm
  IMPLICIT NONE

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
           WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
           !CALL wrf_message( mess )
           LUMATCH=1
        ELSE
           call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
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
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE=SOILPARM_FILE,FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

     WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICAION = ', TRIM ( MMINSL )
     !CALL wrf_message( mess )

     LUMATCH=0

     READ (19,*)
     READ (19,2000,END=2003)SLTYPE
2000 FORMAT (A4)
     READ (19,*)SLCATS,IINDEX
     IF(SLTYPE.EQ.MMINSL)THEN
        WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
             SLCATS,' CATEGORIES'
        !CALL wrf_message ( mess )
        LUMATCH=1
     ENDIF
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(BB    ) < SLCATS .OR. &
          SIZE(DRYSMC) < SLCATS .OR. &
          SIZE(F11   ) < SLCATS .OR. &
          SIZE(MAXSMC) < SLCATS .OR. &
          SIZE(REFSMC) < SLCATS .OR. &
          SIZE(SATPSI) < SLCATS .OR. &
          SIZE(SATDK ) < SLCATS .OR. &
          SIZE(SATDW ) < SLCATS .OR. &
          SIZE(WLTSMC) < SLCATS .OR. &
          SIZE(QTZ   ) < SLCATS  ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
     ENDIF
     IF(SLTYPE.EQ.MMINSL)THEN
        DO LC=1,SLCATS
           READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
                REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                WLTSMC(LC), QTZ(LC)
        ENDDO
     ENDIF

2003 CONTINUE

     CLOSE (19)
  ENDIF

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

logical function wrf_dm_on_monitor() result(l)
  l = .TRUE.
  return
end function wrf_dm_on_monitor

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
  !print*,DEEPRECH
  !if (DEEPRECH .gt. 10.0)stop
  
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
              else
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

SUBROUTINE INITIALIZE_GROUNDWATER(NSOIL,DELTAT,ZSOIL,DZS,WTD,SMC,SH2O,SMCEQ,SMCWTD,&
           DEEPRECH,RECH,QSPRING,DZWT)

! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------

    INTEGER, INTENT(IN)                              :: NSOIL
    REAL,   INTENT(IN)                               :: DELTAT,DZWT
    REAL,    INTENT(IN), DIMENSION(1:NSOIL)          :: ZSOIL,DZS
    REAL, INTENT(INOUT) :: WTD
    REAL,INTENT(INOUT),dimension(1:NSOIL) :: SMC,SH2O,SMCEQ
    REAL,    INTENT(INOUT) :: SMCWTD,DEEPRECH,RECH,QSPRING
    INTEGER  :: K,ITER
    REAL :: BX,SMCMAX,PSISAT,SMCWLT,DWSAT,DKSAT
    REAL :: FRLIQ,SMCEQDEEP
    REAL :: AA,BBB,CC,DD,DX,FUNC,DFUNC,DDZ,EXPON,SMV,FLUX

            !initialize equilibrium soil moisture for water table diagnostic
            CALL EQSMOISTURE(NSOIL ,  ZSOIL , SMCMAX , SMCWLT ,DWSAT, DKSAT  ,BX  , & !in
                                     SMCEQ                          )  !out

            !make sure that below the water table the layers are saturated and initialize the deep soil moisture
            IF(WTD < ZSOIL(NSOIL)-DZS(NSOIL)) THEN

!initialize deep soil moisture so that the flux compensates qlat+qrf
!use Newton-Raphson method to find soil moisture

                         EXPON = 2. * BX + 3.
                         DDZ = ZSOIL(NSOIL) - WTD
                         CC = PSISAT/DDZ
                         FLUX = (DZWT)/DELTAT

                         SMV = 0.5 * SMCMAX

                         DO ITER = 1, 100
                           DD = (SMV+SMCMAX)/(2.*SMCMAX)
                           AA = -DKSAT * DD  ** EXPON
                           BBB = CC * ( (SMCMAX/SMV)**BX - 1. ) + 1.
                           FUNC =  AA * BBB - FLUX
                           DFUNC = -DKSAT * (EXPON/(2.*SMCMAX)) * DD ** (EXPON - 1.) * BBB &
                                   + AA * CC * (-BX) * SMCMAX ** BX * SMV ** (-BX-1.)

                           DX = FUNC/DFUNC
                           SMV = SMV - DX
                           IF ( ABS (DX) < 1.E-6)EXIT
                         ENDDO

                  SMCWTD = MAX(SMV,1.E-4)

             ELSEIF(WTD < ZSOIL(NSOIL))THEN
                  SMCEQDEEP = SMCMAX * ( PSISAT / ( PSISAT - DZS(NSOIL) ) ) ** (1./BX)
                  SMCEQDEEP = MAX(SMCEQDEEP,1.E-4)
                  SMCWTD = SMCMAX * ( WTD -  (ZSOIL(NSOIL)-DZS(NSOIL))) + &
                                  SMCEQDEEP * (ZSOIL(NSOIL) - WTD)

             ELSE !water table within the resolved layers
                  SMCWTD = SMCMAX
                  DO K=NSOIL,2,-1
                     IF(WTD .GE. ZSOIL(K-1))THEN
                          FRLIQ = SH2O(K) / SMC(K)
                          SMC(K) = SMCMAX
                          SH2O(K) = SMCMAX * FRLIQ
                     ELSE
                          IF(SMC(K).LT.SMCEQ(K))THEN
                              WTD = ZSOIL(K)
                          ELSE
                              WTD = ( SMC(K)*DZS(K) - SMCEQ(K)*ZSOIL(K-1) + SMCMAX*ZSOIL(K) ) / &
                                         (SMCMAX - SMCEQ(K))
                          ENDIF
                          EXIT
                     ENDIF
                  ENDDO
             ENDIF

!zero out some variables

             DEEPRECH = 0.
             RECH = 0.
             QSPRING = 0.

END SUBROUTINE INITIALIZE_GROUNDWATER

subroutine Calculate_Deficit(si,smcmax,smcwtd,zwt,sldpth,smc,nsoil)

 implicit none
 integer :: nsoil
 real,intent(out) :: si
 real,intent(in) :: smcmax,smcwtd,zwt
 real,intent(in),dimension(nsoil) :: sldpth,smc
 integer :: isoil

 si = 0.0
 !Iterate though all soil layers
 do isoil = 1,nsoil

  !Calculate the empty space
  si = si + smcmax*sldpth(isoil) - smc(isoil)*sldpth(isoil)

 enddo

 !Add in the transmission zone deficit
 if (abs(zwt) .gt. sum(sldpth)) si = si + (abs(zwt) - sum(sldpth))*(smcmax - smcwtd)

end subroutine Calculate_Deficit

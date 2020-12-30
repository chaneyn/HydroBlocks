subroutine initialize_parameters(MMINLU,SNOW,SNOWH,CANWAT,ISLTYP,IVGTYP,XLAT, &
       TSLB,SMOIS,SH2O,DZS,FNDSOILW,FNDSNOWH,&
       TSK,isnowxy,tvxy,tgxy,canicexy,TMN,XICE,&
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy  ,wslakexy ,zwtxy    ,waxy     , &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy,croptype ,cropcat  ,t2mvxy   ,t2mbxy   ,chstarxy,            &
       restart,                 &
       allowed_to_read , iopt_run,  iopt_crop,                        &
       sf_urban_physics,                         &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &
       rechclim,                                                             &
       gecros_state,nsoil,nhru)

 USE module_sf_noahmpdrv, only: noahmp_init, noahmplsm, noahmp_urban,groundwater_init

 implicit none

    INTEGER, INTENT(INOUT)    ::     nhru
    INTEGER, INTENT(INOUT)       ::     NSOIL, iopt_run, iopt_crop
    LOGICAL      ::     restart,allowed_to_read
    INTEGER, INTENT(INOUT)       ::     sf_urban_physics
    REAL,    DIMENSION( NSOIL), INTENT(INOUT)    ::     DZS  ! Thickness of the soil layers [m]
    REAL,    INTENT(INOUT) ::     DX, DY
    REAL,    DIMENSION(nhru) ,  INTENT(INOUT) ::MSFTX,MSFTY
    REAL,    DIMENSION(nhru, NSOIL) ,INTENT(INOUT)::SMOIS,SH2O,TSLB
    REAL,    DIMENSION(nhru),INTENT(INOUT) :: SNOW,SNOWH,CANWAT
  INTEGER, DIMENSION(nhru),                      &
         &   INTENT(INOUT)       ::     ISLTYP,  &
                                     IVGTYP

    LOGICAL, INTENT(INOUT)       ::     FNDSOILW,FNDSNOWH

    character(len=*), intent(in) :: MMINLU
    REAL, DIMENSION(nhru), INTENT(INOUT) :: XLAT         !latitude
    REAL, DIMENSION(nhru), INTENT(INOUT) :: TSK         !skin temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: TMN         !deep soil temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: XICE         !sea ice fraction
    INTEGER, DIMENSION(nhru), INTENT(INOUT) :: isnowxy     !actual no. of snow layers
    REAL, DIMENSION(nhru,-2:NSOIL), INTENT(INOUT) :: zsnsoxy  !snow layer depth [m]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: tsnoxy !snow temperature [K]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: snicexy  !snow layer ice [mm]
    REAL, DIMENSION(nhru,-2:0), INTENT(INOUT) :: snliqxy  !snow layer liquid water [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tvxy        !vegetation canopy temperature
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tgxy        !ground surface temperature
    REAL, DIMENSION(nhru), INTENT(INOUT) :: canicexy !canopy-intercepted ice (mm)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: canliqxy !canopy-intercepted liquid water (mm)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: eahxy       !canopy air vapor pressure (pa)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: tahxy       !canopy air temperature (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: cmxy        !momentum drag coefficient
    REAL, DIMENSION(nhru), INTENT(INOUT) :: chxy        !sensible heat exchange coefficient
    REAL, DIMENSION(nhru), INTENT(INOUT) :: fwetxy      !wetted or snowed fraction of the canopy (-)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: sneqvoxy    !snow mass at last time step(mm h2o)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: alboldxy    !snow albedo at last time step (-)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: qsnowxy     !snowfall on the ground [mm/s]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: wslakexy    !lake water storage [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: zwtxy       !water table depth [m]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: waxy        !water in the "aquifer" [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: wtxy        !groundwater storage [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: lfmassxy    !leaf mass [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: rtmassxy    !mass of fine roots [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: stmassxy    !stem mass [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: woodxy      !mass of wood (incl. woody roots) [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: grainxy     !mass of grain [g/m2] !XING
    REAL, DIMENSION(nhru), INTENT(INOUT) :: gddxy       !growing degree days !XING
    REAL, DIMENSION(nhru), INTENT(INOUT) :: stblcpxy    !stable carbon in deep soil [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: fastcpxy    !short-lived carbon, shallow soil [g/m2]
    REAL, DIMENSION(nhru), INTENT(INOUT) :: xsaixy      !stem area index
    REAL, DIMENSION(nhru), INTENT(INOUT) :: lai         !leaf area index
    INTEGER, DIMENSION(nhru), INTENT(INOUT) :: cropcat
    REAL   , DIMENSION(nhru,5), INTENT(IN ) :: croptype
! IOPT_RUN = 5 option

    REAL, DIMENSION(nhru,1:nsoil), INTENT(INOUT)  :: smoiseq !equilibrium soil moisture content [m3m-3]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: smcwtdxy !deep soil moisture content [m3m-3]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: deeprechxy !deep recharge [m]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: rechxy !accumulated recharge [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qrfsxy !accumulated flux from groundwater to rivers [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qspringsxy !accumulated seeping water [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: qslatxy !accumulated lateral flow [mm]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: areaxy !grid cell area [m2]
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: FDEPTHXY !efolding depth for transmissivity (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: HT !terrain height (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: RIVERBEDXY !riverbed depth (m)
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: EQZWT !equilibrium water table depth (m)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: RIVERCONDXY !river conductance
    REAL, DIMENSION(nhru), INTENT(INOUT) :: PEXPXY !factor for river conductance
    REAL, DIMENSION(nhru), INTENT(INOUT)  :: rechclim
    REAL, DIMENSION(nhru,60), INTENT(INOUT) :: gecros_state ! Optional gecros crop
    INTEGER,  INTENT(INOUT)  :: STEPWTD
    REAL, INTENT(INOUT)  :: DT, WTDDT
!jref:start
    REAL, DIMENSION(nhru), INTENT(INOUT) :: t2mvxy        !2m temperature vegetation part (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: t2mbxy        !2m temperature bare ground part (k)
    REAL, DIMENSION(nhru), INTENT(INOUT) :: chstarxy        !dummy
!jref:end

integer :: ids,ide, jds,jde, kds,kde, ims,ime, jms,jme, kms,kme, its,ite
integer :: jts,jte, kts,kte

ids = 1
ide = nhru+1
jds = 1
jde = 2
kds = 1
kde = 2
ims = 1
ime = nhru
jms = 1
jme = 1
kms = 1
kme = 2
its = 1
ite = nhru
jts = 1
jte = 1
kts = 1
kte = 2

!Copy the 


 ! Call NOAHM_INIT
print*,'IN'
print*,'LLANDUSE ',MMINLU
print*,'SNOW ',SNOW
print*,'SNOWH ',SNOWH
print*,'CANWAT ',CANWAT
print*,'ISLTYP ',ISLTYP
print*,'IVGTYP ',IVGTYP
print*,'XLAT ',XLAT
print*,'TSLB ',TSLB
print*,'SMOIS ',SMOIS
print*,'SH2O ',SH2O
print*,'DZS ',DZS
print*,'FNDSOILW ',FNDSOILW
print*,'FNDSNOWH ',FNDSNOWH
print*,'TSK ',TSK
print*,'ISNOWXY ',ISNOWXY
print*,'TVXY ',TVXY
print*,'TGXY ',TGXY
print*,'CANICEXY ',CANICEXY
print*,'TMN ',TMN
print*,'XICE ',XICE
print*,'CANLIQXY',CANLIQXY
print*,'EAHXY',EAHXY
print*,'TAHXY',TAHXY
print*,'CMXY',CMXY
print*,'CHXY',CHXY
print*,'FWETXY',FWETXY
print*,'SNEQVOXY',SNEQVOXY
print*,'ALBOLDXY',ALBOLDXY
print*,'QSNOWXY',QSNOWXY
print*,'WSLAKEXY',WSLAKEXY
print*,'ZWTXY',ZWTXY
print*,'WAXY',WAXY
print*,'WOODXY',WOODXY
print*,'STBLCPXY',STBLCPXY
print*,'FASTCPXY',FASTCPXY
print*,'XSAIXY',XSAIXY
print*,'LAI',LAI
print*,'GRAINXY',GRAINXY
print*,'GDDXY',GDDXY
print*,'CROPTYPE',CROPTYPE
print*,'CROPCAT',CROPCAT
print*,'T2MVXY',T2MVXY
print*,'T2MBXY',T2MBXY
print*,'CHSTARXY',CHSTARXY
print*,'NSOIL',NSOIL
print*,'runoff_option',iopt_run
print*,'crop_option',iopt_crop
print*,'sf_urban_physics',sf_urban_physics
print*,'ids',ids
print*,'ide',ide
print*,'jds',jds
print*,'jde',jde
print*,'kds',kds
print*,'kde',kde
print*,'ims',ims
print*,'ime',ime
print*,'jms',jms
print*,'jme',jme
print*,'kms',kms
print*,'kme',kme
print*,'its',its
print*,'ite',ite
print*,'jts',jts
print*,'jte',jte
print*,'kts',kts
print*,'kte',kte
print*,'smoiseq',smoiseq
print*,'smcwtdxy',smcwtdxy
print*,'rechxy',rechxy
print*,'deeprechxy',deeprechxy
print*,'areaxy',areaxy
print*,'dx',dx
print*,'dy',dy
print*,'msftx',msftx
print*,'msfty',msfty
print*,'wtddt',wtddt
print*,'stepwtd',stepwtd
print*,'dt',dt
print*,'qrfsxy',qrfsxy
print*,'qspringsxy',qspringsxy
print*,'qslatxy',qslatxy
print*,'fdepthxy',fdepthxy
print*,'terrain',ht
print*,'riverbedxy',riverbedxy
print*,'eqzwt',eqzwt
print*,'rivercondxy',rivercondxy
print*,'pexpxy',pexpxy
print*,'rechclim',rechclim
print*,'gecros_state',gecros_state

 call NOAHMP_INIT ( MMINLU, SNOW , SNOWH , CANWAT , ISLTYP ,   IVGTYP,XLAT, &
       TSLB , SMOIS , SH2O , DZS , FNDSOILW , FNDSNOWH ,             &
       TSK, isnowxy , tvxy     ,tgxy     ,canicexy ,         TMN,     XICE,   &
       canliqxy ,eahxy    ,tahxy    ,cmxy     ,chxy     ,                     &
       fwetxy   ,sneqvoxy ,alboldxy ,qsnowxy  ,wslakexy ,zwtxy    ,waxy     , &
       wtxy     ,tsnoxy   ,zsnsoxy  ,snicexy  ,snliqxy  ,lfmassxy ,rtmassxy , &
       stmassxy ,woodxy   ,stblcpxy ,fastcpxy ,xsaixy   ,lai      ,           &
       grainxy  ,gddxy    ,                                                   &
       croptype ,cropcat  ,                      &
       t2mvxy   ,t2mbxy   ,chstarxy,            &
       NSOIL, restart,                 &
       allowed_to_read , iopt_run,  iopt_crop,                        &
       sf_urban_physics,                         &
       ids,ide, jds,jde, kds,kde,                &
       ims,ime, jms,jme, kms,kme,                &
       its,ite, jts,jte, kts,kte,                &
       smoiseq  ,smcwtdxy ,rechxy   ,deeprechxy, areaxy, dx, dy, msftx, msfty,&
       wtddt    ,stepwtd  ,dt       ,qrfsxy     ,qspringsxy  , qslatxy    ,  &
       fdepthxy ,ht     ,riverbedxy ,eqzwt     ,rivercondxy ,pexpxy       ,  &
       rechclim,                                                             &
       gecros_state)

print*,'OUT'
print*,'LLANDUSE ',MMINLU
print*,'SNOW ',SNOW
print*,'SNOWH ',SNOWH
print*,'CANWAT ',CANWAT
print*,'ISLTYP ',ISLTYP
print*,'IVGTYP ',IVGTYP
print*,'XLAT ',XLAT
print*,'TSLB ',TSLB
print*,'SMOIS ',SMOIS
print*,'SH2O ',SH2O
print*,'DZS ',DZS
print*,'FNDSOILW ',FNDSOILW
print*,'FNDSNOWH ',FNDSNOWH
print*,'TSK ',TSK
print*,'ISNOWXY ',ISNOWXY
print*,'TVXY ',TVXY
print*,'TGXY ',TGXY
print*,'CANICEXY ',CANICEXY
print*,'TMN ',TMN
print*,'XICE ',XICE
print*,'CANLIQXY',CANLIQXY
print*,'EAHXY',EAHXY
print*,'TAHXY',TAHXY
print*,'CMXY',CMXY
print*,'CHXY',CHXY
print*,'FWETXY',FWETXY
print*,'SNEQVOXY',SNEQVOXY
print*,'ALBOLDXY',ALBOLDXY
print*,'QSNOWXY',QSNOWXY
print*,'WSLAKEXY',WSLAKEXY
print*,'ZWTXY',ZWTXY
print*,'WAXY',WAXY
print*,'WOODXY',WOODXY
print*,'STBLCPXY',STBLCPXY
print*,'FASTCPXY',FASTCPXY
print*,'XSAIXY',XSAIXY
print*,'LAI',LAI
print*,'GRAINXY',GRAINXY
print*,'GDDXY',GDDXY
print*,'CROPTYPE',CROPTYPE
print*,'CROPCAT',CROPCAT
print*,'T2MVXY',T2MVXY
print*,'T2MBXY',T2MBXY
print*,'CHSTARXY',CHSTARXY
print*,'NSOIL',NSOIL
print*,'runoff_option',iopt_run
print*,'crop_option',iopt_crop
print*,'sf_urban_physics',sf_urban_physics
print*,'ids',ids
print*,'ide',ide
print*,'jds',jds
print*,'jde',jde
print*,'kds',kds
print*,'kde',kde
print*,'ims',ims
print*,'ime',ime
print*,'jms',jms
print*,'jme',jme
print*,'kms',kms
print*,'kme',kme
print*,'its',its
print*,'ite',ite
print*,'jts',jts
print*,'jte',jte
print*,'kts',kts
print*,'kte',kte
print*,'smoiseq',smoiseq
print*,'smcwtdxy',smcwtdxy
print*,'rechxy',rechxy
print*,'deeprechxy',deeprechxy
print*,'areaxy',areaxy
print*,'dx',dx
print*,'dy',dy
print*,'msftx',msftx
print*,'msfty',msfty
print*,'wtddt',wtddt
print*,'stepwtd',stepwtd
print*,'dt',dt
print*,'qrfsxy',qrfsxy
print*,'qspringsxy',qspringsxy
print*,'qslatxy',qslatxy
print*,'fdepthxy',fdepthxy
print*,'terrain',ht
print*,'riverbedxy',riverbedxy
print*,'eqzwt',eqzwt
print*,'rivercondxy',rivercondxy
print*,'pexpxy',pexpxy
print*,'rechclim',rechclim
print*,'gecros_state',gecros_state

end subroutine


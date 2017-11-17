module SolarRadiationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar radiation transfer through canopy
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use decompMod, only : bounds_type
  use pftconMod, only : pftcon
  use PatchType, only : patch
  use SurfaceAlbedoType, only : surfalb_type
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SolarRadiation           ! Main driver for radiative transfer
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: NormanRadiation         ! Norman radiative transfer
  private :: GoudriaanRadiation      ! Goudriaan radiative transfer
  private :: TwoStreamRadiation      ! Two-stream approximation radiative transfer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SolarRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
  surfalb_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Solar radiation transfer through canopy
    !
    ! !USES:
    use clm_varpar, only : numrad, nlevcan, isun, isha, ivis
    use clm_varcon, only : pi => rpi
    use clm_varctl, only : light
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp               ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)         ! CLM patch filter for non-snow-covered vegetation
    type(surfalb_type), intent(in) :: surfalb_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                        ! Filter index
    integer  :: p                                        ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                        ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                                       ! Aboveground layer index
    integer  :: ib                                       ! Waveband index
    integer  :: lev                                      ! Canopy level index
    integer  :: j                                        ! Sky angle index
    real(r8) :: angle                                    ! Sky angles (5, 15, 25, 35, 45, 55, 65, 75, 85 degrees)
    real(r8) :: gdirj                                    ! Relative projected area of leaf elements in the direction of sky angle
    real(r8) :: cumpai                                   ! Cumulative plant area index (m2/m2)
    real(r8) :: chil(bounds%begp:bounds%endp)            ! Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
    real(r8) :: phi1(bounds%begp:bounds%endp)            ! Term in Ross-Goudriaan function for gdir
    real(r8) :: phi2(bounds%begp:bounds%endp)            ! Term in Ross-Goudriaan function for gdir
    real(r8) :: gdir(bounds%begp:bounds%endp)            ! Relative projected area of leaf elements in the direction of solar beam
    real(r8) :: kb(bounds%begp:bounds%endp)              ! Direct beam extinction coefficient
    real(r8) :: wl(bounds%begp:bounds%endp)              ! Leaf fraction of canopy
    real(r8) :: ws(bounds%begp:bounds%endp)              ! Stem fraction of canopy
    real(r8) :: rho(bounds%begp:bounds%endp,1:numrad)    ! Leaf/stem reflectance
    real(r8) :: tau(bounds%begp:bounds%endp,1:numrad)    ! Leaf/stem transmittance
    real(r8) :: omega(bounds%begp:bounds%endp,1:numrad)  ! Leaf/stem scattering coefficient

    ! For Norman radiation
    real(r8) :: tbj(bounds%begp:bounds%endp,0:nlevcan)   ! Exponential transmittance of direct beam onto canopy layer
    real(r8) :: tb(bounds%begp:bounds%endp,1:nlevcan)    ! Exponential transmittance of direct beam through a single leaf layer

    ! For Goudriaan radiation
    real(r8) :: albvegh                                    ! Vegetation albedo, horizontal leaves
    real(r8) :: albvegb                                    ! Direct beam vegetation albedo, non-horizontal leaves
    real(r8) :: albvegd                                    ! Diffuse vegetation albedo, non-horizontal leaves
    real(r8) :: kbj                                        ! kb for sky angle j
    real(r8) :: albvegbj                                   ! albvegb for sky angle j
    real(r8) :: kd(bounds%begp:bounds%endp)                ! Diffuse radiation extinction coefficient
    real(r8) :: kbm(bounds%begp:bounds%endp,1:numrad)      ! kb adjusted for scattering
    real(r8) :: kdm(bounds%begp:bounds%endp,1:numrad)      ! kd adjusted for scattering
    real(r8) :: albcanb(bounds%begp:bounds%endp,1:numrad)  ! Direct beam albedo above canopy
    real(r8) :: albcand(bounds%begp:bounds%endp,1:numrad)  ! Diffuse albedo above canopy

    ! For two-stream radiation
    real(r8) :: asu                                        ! Single scattering albedo
    real(r8) :: tmp0,tmp1,tmp2                             ! Intermediate variables
    real(r8) :: avmu(bounds%begp:bounds%endp)              ! Average inverse diffuse optical depth per unit leaf area
    real(r8) :: betad(bounds%begp:bounds%endp,1:numrad)    ! Upscatter parameter for diffuse radiation
    real(r8) :: betab(bounds%begp:bounds%endp,1:numrad)    ! Upscatter parameter for direct beam radiation
    !---------------------------------------------------------------------

    associate ( &
                                                ! *** Input ***
    xl         => pftcon%xl                , &  ! Departure of leaf angle from spherical orientation (-)
    clump_fac  => pftcon%clump_fac         , &  ! Foliage clumping index (-)
    rhol       => pftcon%rhol              , &  ! Leaf reflectance (-)
    taul       => pftcon%taul              , &  ! Leaf transmittance (-)
    rhos       => pftcon%rhos              , &  ! Stem reflectance (-)
    taus       => pftcon%taus              , &  ! Stem transmittance (-)
    albsoib    => surfalb_inst%albgrd_col  , &  ! Direct beam albedo of ground (soil)
    albsoid    => surfalb_inst%albgri_col  , &  ! Diffuse albedo of ground (soil)
    lai        => mlcanopy_inst%lai        , &  ! Leaf area index of canopy (m2/m2)
    sai        => mlcanopy_inst%sai        , &  ! Stem area index of canopy (m2/m2)
    ncan       => mlcanopy_inst%ncan       , &  ! Number of aboveground layers
    nbot       => mlcanopy_inst%nbot       , &  ! Index for bottom leaf layer
    ntop       => mlcanopy_inst%ntop       , &  ! Index for top leaf layer
    dpai       => mlcanopy_inst%dpai       , &  ! Layer plant area index (m2/m2)
    sumpai     => mlcanopy_inst%sumpai     , &  ! Cumulative plant area index (m2/m2)
    solar_zen  => mlcanopy_inst%solar_zen  , &  ! Solar zenith angle (radians)
                                                ! *** Output ***
    fracsun    => mlcanopy_inst%fracsun    , &  ! Sunlit fraction of canopy layer
    fracsha    => mlcanopy_inst%fracsha    , &  ! Shaded fraction of canopy layer
    swleaf     => mlcanopy_inst%swleaf     , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    apar       => mlcanopy_inst%apar       , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    td         => mlcanopy_inst%td           &  ! Exponential transmittance of diffuse radiation through a single leaf layer
    )

    !---------------------------------------------------------------------
    ! Weight reflectance and transmittance by lai and sai and calculate
    ! leaf scattering coefficient
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       if ((lai(p)+sai(p)) > 0._r8) then
          wl(p) = lai(p) / (lai(p)+sai(p))
          ws(p) = sai(p) / (lai(p)+sai(p))
       else
          wl(p) = 0._r8
          ws(p) = 0._r8
       end if
    end do

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          rho(p,ib) = max(rhol(patch%itype(p),ib)*wl(p) + rhos(patch%itype(p),ib)*ws(p), 1.e-06_r8)
          tau(p,ib) = max(taul(patch%itype(p),ib)*wl(p) + taus(patch%itype(p),ib)*ws(p), 1.e-06_r8)
          omega(p,ib) = rho(p,ib) + tau(p,ib)
       end do
    end do

    !---------------------------------------------------------------------
    ! Direct beam extinction coefficient
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       chil(p) = min(max(xl(patch%itype(p)), -0.4_r8), 0.6_r8)
       if (abs(chil(p)) <= 0.01_r8) chil(p) = 0.01_r8

       phi1(p) = 0.5_r8 - 0.633_r8*chil(p) - 0.330_r8*chil(p)*chil(p)
       phi2(p) = 0.877_r8 * (1._r8 - 2._r8*phi1(p))

!      phi1(p) = 0._r8       ! Horizontal leaves
!      phi2(p) = 1._r8       ! Horizontal leaves

       gdir(p) = phi1(p) + phi2(p) * cos(solar_zen(p))
       kb(p) = gdir(p) / cos(solar_zen(p))
       kb(p) = min(kb(p), 40._r8)

    end do

    !---------------------------------------------------------------------
    ! Sunlit and shaded fraction of leaf layer
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Zero out for all layers

       do ic = 1, ncan(p)
          fracsun(p,ic) = 0._r8
          fracsha(p,ic) = 0._r8
       end do

       ! Leaf layers

       do ic = nbot(p), ntop(p)
          fracsun(p,ic) = clump_fac(patch%itype(p)) * exp(-kb(p) * sumpai(p,ic) * clump_fac(patch%itype(p)))
          fracsha(p,ic) = 1._r8 - fracsun(p,ic)
       end do

    end do

    !---------------------------------------------------------------------
    ! Diffuse transmittance for a single layer estimated for nine sky
    ! angles in increments of 10 degrees (also needed for longwave radiation)
    !---------------------------------------------------------------------

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Zero out for all layers

       do ic = 1, ncan(p)
          td(p,ic) = 0._r8
       end do

       ! Leaf layers

       do ic = nbot(p), ntop(p)
          do j = 1, 9
             angle = (5._r8 + (j - 1) * 10._r8) * pi / 180._r8
             gdirj = phi1(p) + phi2(p) * cos(angle)
             td(p,ic) = td(p,ic) &
                      + exp(-gdirj / cos(angle) * dpai(p,ic) * clump_fac(patch%itype(p))) * sin(angle) * cos(angle)
          end do
          td(p,ic) = td(p,ic) * 2._r8 * (10._r8 * pi / 180._r8)
       end do

    end do

    !---------------------------------------------------------------------
    ! Parameters for Norman radiative transfer
    !---------------------------------------------------------------------

    ! Direct beam transmittance (tb) for a single leaf layer

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Zero out for all layers

       do ic = 1, ncan(p)
          tb(p,ic) = 0._r8
       end do

       ! Leaf layers

       do ic = nbot(p), ntop(p)
          tb(p,ic) = exp(-kb(p) * dpai(p,ic) * clump_fac(patch%itype(p)))
       end do

    end do

    ! Direct beam transmittance (tbj) uses cumulative plant area index
    ! above layer j to give unscattered direct beam onto layer j

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Zero out for all layers

       do ic = 0, ncan(p)
          tbj(p,ic) = 0._r8
       end do

       ! Leaf layers and ground

       cumpai = 0._r8
       tbj(p,ntop(p)) = 1._r8
       do ic = ntop(p), nbot(p), -1
          cumpai = cumpai + dpai(p,ic)
          if (ic > nbot(p)) then
             lev = ic - 1   ! transmittance onto leaf layer below
          else if (ic == nbot(p)) then
             lev = 0        ! transmittance onto ground
          end if
          tbj(p,lev) = exp(-kb(p) * cumpai * clump_fac(patch%itype(p)))
       end do
    end do

    !---------------------------------------------------------------------
    ! Parameters for Goudriaan radiative transfer
    !---------------------------------------------------------------------

    ! Diffuse extinction coefficient for canopy

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       kd(p) = 0._r8
       do j = 1, 9
          angle = (5._r8 + (j - 1) * 10._r8) * pi / 180._r8
          gdirj = phi1(p) + phi2(p) * cos(angle)
          kd(p) = kd(p) + exp(-gdirj / cos(angle) * (lai(p)+sai(p))) * sin(angle) * cos(angle)
       end do
       kd(p) = kd(p) * 2._r8 * (10._r8 * pi / 180._r8)
       if ((lai(p)+sai(p)) > 0._r8) then
          kd(p) = -log(kd(p)) / (lai(p)+sai(p))
       else
          kd(p) = 0._r8
       end if
    end do

    ! Adjust extinction coefficents for scattering and calculate canopy albedos

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          c = patch%column(p)

          ! Adjust for scattering

          kbm(p,ib) = kb(p) * sqrt(1._r8 - omega(p,ib))
          kdm(p,ib) = kd(p) * sqrt(1._r8 - omega(p,ib))

          ! Vegetation albedo, horizontal leaves

          albvegh = (1._r8 - sqrt(1._r8 - omega(p,ib))) / (1._r8 + sqrt(1._r8 - omega(p,ib)))

          ! Direct beam vegetation albedo, non-horizontal leaves
      
          albvegb = 2._r8 * kb(p) / (kb(p) + kd(p)) * albvegh

          ! Diffuse vegetation albedo, non-horizontal leaves

          albvegd = 0._r8
          do j = 1, 9
             angle = (5._r8 + (j - 1) * 10._r8) * pi / 180._r8
             gdirj = phi1(p) + phi2(p) * cos(angle)
             kbj = gdirj / cos(angle)
             albvegbj = 2._r8 * kbj / (kbj + kd(p)) * albvegh
             albvegd = albvegd + albvegbj * sin(angle) * cos(angle)
          end do
          albvegd = albvegd * 2._r8 * (10._r8 * pi / 180._r8)
   
          ! Effective canopy albedo, including soil
      
          albcanb(p,ib) = albvegb + (albsoib(c,ib) - albvegb) &
                        * exp(-2._r8 * kbm(p,ib) * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))
          albcand(p,ib) = albvegd + (albsoid(c,ib) - albvegd) &
                        * exp(-2._r8 * kdm(p,ib) * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))

       end do
    end do

    !---------------------------------------------------------------------
    ! Parameters for two-stream radiative transfer
    !---------------------------------------------------------------------

    ! avmu - average inverse diffuse optical depth per unit leaf area

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       avmu(p) = ( 1._r8 - phi1(p)/phi2(p) * log((phi1(p)+phi2(p))/phi1(p)) ) / phi2(p)
    end do

    ! betad - upscatter parameter for diffuse radiation
    ! betab - upscatter parameter for direct beam radiation

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)

          ! upscatter parameter for diffuse radiation

          betad(p,ib) = 0.5_r8 / omega(p,ib) &
                      * ( rho(p,ib) + tau(p,ib) + (rho(p,ib)-tau(p,ib)) * ((1._r8+chil(p))/2._r8)**2 )

          ! upscatter parameter for direct beam radiation

          tmp0 = gdir(p) + phi2(p) * cos(solar_zen(p))
          tmp1 = phi1(p) * cos(solar_zen(p))
          tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
          asu = 0.5_r8 * omega(p,ib) * gdir(p) / tmp0 * tmp2
          betab(p,ib) = (1._r8 + avmu(p)*kb(p)) / (omega(p,ib)*avmu(p)*kb(p)) * asu

       end do
    end do

    !---------------------------------------------------------------------
    ! Calculate radiative transfer through canopy
    !---------------------------------------------------------------------

    if (light == 1) then
       call NormanRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
       rho, tau, omega, tb, td, tbj, surfalb_inst, mlcanopy_inst)
    else if (light == 2) then
       call GoudriaanRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
       omega, kb, kbm, kdm, albcanb, albcand, mlcanopy_inst)
    else if (light == 3) then
       call TwoStreamRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
       omega, kb, avmu, betad, betab, surfalb_inst, mlcanopy_inst)
    end if

    ! APAR per unit sunlit and shaded leaf area

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       do ic = 1, ncan(p)
          apar(p,ic,isun) = swleaf(p,ic,isun,ivis) * 4.6_r8
          apar(p,ic,isha) = swleaf(p,ic,isha,ivis) * 4.6_r8
       end do
    end do

    end associate
  end subroutine SolarRadiation

  !-----------------------------------------------------------------------
  subroutine NormanRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
  rho, tau, omega, tb, td, tbj, surfalb_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using Norman (1979)
    !
    ! !USES:
    use clm_varpar, only : numrad, nlevcan, isun, isha
    use clm_varctl, only : iulog
    use MathToolsMod, only : tridiag
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in) :: num_exposedvegp                          ! Number of non-snow-covered veg points in CLM patch filter
    integer,  intent(in) :: filter_exposedvegp(:)                    ! CLM patch filter for non-snow-covered vegetation
    real(r8), intent(in) :: rho(bounds%begp:bounds%endp,1:numrad)    ! Leaf/stem reflectance
    real(r8), intent(in) :: tau(bounds%begp:bounds%endp,1:numrad)    ! Leaf/stem transmittance
    real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:numrad)  ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: tb(bounds%begp:bounds%endp,1:nlevcan)    ! Exponential transmittance of direct beam through a single leaf layer
    real(r8), intent(in) :: td(bounds%begp:bounds%endp,1:nlevcan)    ! Exponential transmittance of diffuse through a single leaf layer
    real(r8), intent(in) :: tbj(bounds%begp:bounds%endp,0:nlevcan)   ! Exponential transmittance of direct beam onto canopy layer
    type(surfalb_type), intent(in) :: surfalb_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                                ! Filter index
    integer  :: p                                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                                               ! Aboveground layer index
    integer  :: icm1                                             ! Layer below ic (ic-1)
    integer  :: ib                                               ! Waveband index
    real(r8) :: suminc                                           ! Incident radiation for energy conservation check
    real(r8) :: sumref                                           ! Reflected radiation for energy conservation check
    real(r8) :: sumabs                                           ! Absorbed radiation for energy conservation check
    real(r8) :: error                                            ! Error check
    real(r8) :: trand                                            ! Term for diffuse radiation transmitted by layer
    real(r8) :: refld                                            ! Term for diffuse radiation reflected by layer
    integer  :: m                                                ! Index to the tridiagonal matrix
    real(r8) :: aic, bic                                         ! Intermediate terms for tridiagonal matrix
    real(r8) :: eic, fic                                         ! Intermediate terms for tridiagonal matrix
    integer, parameter :: neq = (nlevcan+1)*2                    ! Number of tridiagonal equations to solve
    real(r8) :: atri(neq), btri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: ctri(neq), dtri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: utri(neq)                                        ! Tridiagonal solution
    real(r8) :: swbeam                                           ! Direct beam solar flux onto canopy layer (W/m2 ground)
    real(r8) :: swabsb                                           ! Absorbed direct beam solar flux for canopy layer (W/m2 ground)
    real(r8) :: swabsd                                           ! Absorbed diffuse solar flux for canopy layer (W/m2 ground)
    real(r8) :: swsun                                            ! Absorbed solar radiation, sunlit fraction of layer (W/m2 ground)
    real(r8) :: swsha                                            ! Absorbed solar radiation, shaded fraction of layer (W/m2 ground)
    real(r8) :: swup(bounds%begp:bounds%endp,0:nlevcan,numrad)   ! Upward diffuse solar flux above canopy layer (W/m2 ground)
    real(r8) :: swdn(bounds%begp:bounds%endp,0:nlevcan,numrad)   ! Downward diffuse solar flux onto canopy layer (W/m2 ground)
    !---------------------------------------------------------------------

    associate ( &
                                               ! *** Input ***
    ncan      => mlcanopy_inst%ncan       , &  ! Number of aboveground layers
    nbot      => mlcanopy_inst%nbot       , &  ! Index for bottom leaf layer
    ntop      => mlcanopy_inst%ntop       , &  ! Index for top leaf layer
    dpai      => mlcanopy_inst%dpai       , &  ! Layer plant area index (m2/m2)
    swskyb    => mlcanopy_inst%swskyb     , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd    => mlcanopy_inst%swskyd     , &  ! Atmospheric diffuse solar radiation (W/m2)
    fracsun   => mlcanopy_inst%fracsun    , &  ! Sunlit fraction of canopy layer
    fracsha   => mlcanopy_inst%fracsha    , &  ! Shaded fraction of canopy layer
    albsoib   => surfalb_inst%albgrd_col  , &  ! Direct beam albedo of ground (soil)
    albsoid   => surfalb_inst%albgri_col  , &  ! Diffuse albedo of ground (soil)
                                               ! *** Output ***
    swleaf    => mlcanopy_inst%swleaf     , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swveg     => mlcanopy_inst%swveg      , &  ! Absorbed solar radiation, vegetation (W/m2)
    swvegsun  => mlcanopy_inst%swvegsun   , &  ! Absorbed solar radiation, sunlit canopy (W/m2)
    swvegsha  => mlcanopy_inst%swvegsha   , &  ! Absorbed solar radiation, shaded canopy (W/m2)
    albcan    => mlcanopy_inst%albcan     , &  ! Albedo above canopy
    swsoi     => mlcanopy_inst%swsoi        &  ! Absorbed solar radiation, ground (W/m2)
    )

    !---------------------------------------------------------------------
    ! Set up and solve tridiagonal system of equations for radiative fluxes
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          c = patch%column(p)

          ! Zero out radiative fluxes for all layers

          swup(p,0,ib) = 0._r8
          swdn(p,0,ib) = 0._r8

          do ic = 1, ncan(p)
             swup(p,ic,ib) = 0._r8
             swdn(p,ic,ib) = 0._r8
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! There are two equations for each leaf layer and the soil. The first
          ! equation is the upward flux and the second equation is the downward flux.

          m = 0

          ! Soil: upward flux

          m = m + 1 
          atri(m) = 0._r8 
          btri(m) = 1._r8 
          ctri(m) = -albsoid(c,ib) 
          dtri(m) = swskyb(p,ib) * tbj(p,0) * albsoib(c,ib) 

          ! Soil: downward flux

          refld = (1._r8 - td(p,nbot(p))) * rho(p,ib) 
          trand = (1._r8 - td(p,nbot(p))) * tau(p,ib) + td(p,nbot(p)) 
          aic = refld - trand * trand / refld 
          bic = trand / refld 

          m = m + 1 
          atri(m) = -aic 
          btri(m) = 1._r8 
          ctri(m) = -bic 
          dtri(m) = swskyb(p,ib) * tbj(p,nbot(p)) * (1._r8 - tb(p,nbot(p))) * (tau(p,ib) - rho(p,ib) * bic)

          ! Leaf layers, excluding top layer

          do ic = nbot(p), ntop(p)-1

             ! Upward flux

             refld = (1._r8 - td(p,ic)) * rho(p,ib)
             trand = (1._r8 - td(p,ic)) * tau(p,ib) + td(p,ic)
             fic = refld - trand * trand / refld
             eic = trand / refld

             m = m + 1
             atri(m) = -eic
             btri(m) = 1._r8
             ctri(m) = -fic
             dtri(m) = swskyb(p,ib) * tbj(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ib) - tau(p,ib) * eic)

             ! Downward flux

             refld = (1._r8 - td(p,ic+1)) * rho(p,ib)
             trand = (1._r8 - td(p,ic+1)) * tau(p,ib) + td(p,ic+1)
             aic = refld - trand * trand / refld
             bic = trand / refld

             m = m + 1
             atri(m) = -aic
             btri(m) = 1._r8
             ctri(m) = -bic
             dtri(m) = swskyb(p,ib) * tbj(p,ic+1) * (1._r8 - tb(p,ic+1)) * (tau(p,ib) - rho(p,ib) * bic)

          end do

          ! Top canopy layer: upward flux

          ic = ntop(p)
          refld = (1._r8 - td(p,ic)) * rho(p,ib)
          trand = (1._r8 - td(p,ic)) * tau(p,ib) + td(p,ic)
          fic = refld - trand * trand / refld
          eic = trand / refld

          m = m + 1
          atri(m) = -eic
          btri(m) = 1._r8
          ctri(m) = -fic
          dtri(m) = swskyb(p,ib) * tbj(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ib) - tau(p,ib) * eic)

          ! Top canopy layer: downward flux

          m = m + 1
          atri(m) = 0._r8
          btri(m) = 1._r8
          ctri(m) = 0._r8
          dtri(m) = swskyd(p,ib)

          ! Solve tridiagonal system of equations for upward and downward fluxes

          call tridiag (atri, btri, ctri, dtri, utri, m)

          ! Now copy the solution (utri) to the upward (swup) and downward (swdn)
          ! fluxes for each layer
          ! swup =  Upward diffuse flux above layer
          ! swdn =  Downward diffuse flux onto layer

          m = 0

          ! Soil fluxes

          m = m + 1
          swup(p,0,ib) = utri(m)
          m = m + 1
          swdn(p,0,ib) = utri(m)

          ! Leaf layer fluxes

          do ic = nbot(p), ntop(p)
             m = m + 1
             swup(p,ic,ib) = utri(m)
             m = m + 1
             swdn(p,ic,ib) = utri(m)
          end do

       end do             ! end patch loop
    end do                ! end waveband loop

    !---------------------------------------------------------------------
    ! Compute fluxes
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          c = patch%column(p)

          ! Solar radiation absorbed by ground (soil)

          swbeam = tbj(p,0) * swskyb(p,ib)
          swabsb = swbeam * (1._r8 - albsoib(c,ib))
          swabsd = swdn(p,0,ib) * (1._r8 - albsoid(c,ib))
          swsoi(p,ib) = swabsb + swabsd

          ! Leaf layer fluxes

          swveg(p,ib) = 0._r8
          swvegsun(p,ib) = 0._r8
          swvegsha(p,ib) = 0._r8

          do ic = nbot(p), ntop(p)

             ! Downward direct beam incident on layer and absorbed direct
             ! beam and diffuse for layer. Note special case for first
             ! leaf layer, where the upward flux from below is from the ground.
             ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
             ! canopy layers without leaves.

             swbeam = tbj(p,ic) * swskyb(p,ib)
             swabsb = swbeam * (1._r8 - tb(p,ic)) * (1._r8 - omega(p,ib))
             if (ic == nbot(p)) then
                icm1 = 0
             else
                icm1 = ic - 1
             end if
             swabsd = (swdn(p,ic,ib) + swup(p,icm1,ib)) * (1._r8 - td(p,ic)) * (1._r8 - omega(p,ib))

             ! Absorbed radiation for shaded and sunlit portions of layer

             swsha = swabsd * fracsha(p,ic)
             swsun = swabsd * fracsun(p,ic) + swabsb 

             ! Per unit sunlit and shaded leaf area

             swleaf(p,ic,isun,ib) = swsun / (fracsun(p,ic) * dpai(p,ic))
             swleaf(p,ic,isha,ib) = swsha / (fracsha(p,ic) * dpai(p,ic))

             ! Sum solar radiation absorbed by vegetation and sunlit/shaded leaves

             swveg(p,ib) = swveg(p,ib) + (swabsb + swabsd)
             swvegsun(p,ib) = swvegsun(p,ib) + swsun
             swvegsha(p,ib) = swvegsha(p,ib) + swsha

          end do  

          ! Albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          if (suminc > 0._r8) then
             albcan(p,ib) = swup(p,ntop(p),ib) / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Conservation check for total radiation balance: absorbed = incoming - outgoing

          sumref = albcan(p,ib) * (swskyb(p,ib) + swskyd(p,ib))
          sumabs = suminc - sumref
          error = sumabs - (swveg(p,ib) + swsoi(p,ib))
          if (abs(error) > 1.e-03) then
             call endrun (msg='ERROR: NormanRadiation: total solar conservation error')
          end if

          ! Sunlit and shaded absorption

          error = (swvegsun(p,ib) + swvegsha(p,ib)) - swveg(p,ib) 
          if (abs(error) > 1.e-03) then
             call endrun (msg='ERROR: NormanRadiation: sunlit/shade solar conservation error')
          end if

       end do             ! end patch loop
    end do                ! end waveband loop

    end associate
  end subroutine NormanRadiation

  !-----------------------------------------------------------------------
  subroutine GoudriaanRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
  omega, kb, kbm, kdm, albcanb, albcand, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using Goudriaan (1977)
    ! as described by Goudriaan and van Laar (1994)
    !
    ! !USES:
    use clm_varpar, only : numrad, isun, isha
    use clm_varctl, only : iulog
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in) :: num_exposedvegp                            ! Number of non-snow-covered veg points in CLM patch filter
    integer,  intent(in) :: filter_exposedvegp(:)                      ! CLM patch filter for non-snow-covered vegetation
    real(r8), intent(in) :: omega(bounds%begp:bounds%endp,1:numrad)    ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: kb(bounds%begp:bounds%endp)                ! Direct beam extinction coefficient
    real(r8), intent(in) :: kbm(bounds%begp:bounds%endp,1:numrad)      ! kb adjusted for scattering
    real(r8), intent(in) :: kdm(bounds%begp:bounds%endp,1:numrad)      ! kd adjusted for scattering
    real(r8), intent(in) :: albcanb(bounds%begp:bounds%endp,1:numrad)  ! Direct beam albedo above canopy
    real(r8), intent(in) :: albcand(bounds%begp:bounds%endp,1:numrad)  ! Diffuse albedo above canopy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f          ! Filter index
    integer  :: p          ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic         ! Aboveground layer index
    integer  :: ib         ! Waveband index
    real(r8) :: suminc     ! Incident radiation for energy conservation check
    real(r8) :: sumref     ! Reflected radiation for energy conservation check
    real(r8) :: sumabs     ! Absorbed radiation for energy conservation check

    ! Fluxes per unit leaf area (W/m2 leaf)
    real(r8) :: ild        ! Absorbed diffuse flux per unit leaf area (W/m2)
    real(r8) :: ilb        ! Absorbed dir beam (total, with scattering) per unit leaf area (W/m2)
    real(r8) :: ilbb       ! Absorbed direct beam (unscattered) per unit leaf area (W/m2)
    real(r8) :: ilbs       ! Absorbed direct beam (scattered) per unit leaf area (W/m2)
    real(r8) :: ilsun      ! Absorbed solar rad (sunlit leaves) per unit sunlit leaf area (W/m2)
    real(r8) :: ilsha      ! Absorbed solar rad (shaded leaves) per unit shaded leaf area (W/m2)

    ! Fluxes per unit ground area (W/m2 ground area)
    real(r8) :: icsun      ! Absorbed solar radiation, sunlit canopy (W/m2)
    real(r8) :: icsha      ! Absorbed solar radiation, shaded canopy (W/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                 ! *** Input ***
    clump_fac  => pftcon%clump_fac          , &  ! Foliage clumping index (-)
    lai        => mlcanopy_inst%lai         , &  ! Leaf area index of canopy (m2/m2)
    sai        => mlcanopy_inst%sai         , &  ! Stem area index of canopy (m2/m2)
    ncan       => mlcanopy_inst%ncan        , &  ! Number of aboveground layers
    nbot       => mlcanopy_inst%nbot        , &  ! Index for bottom leaf layer
    ntop       => mlcanopy_inst%ntop        , &  ! Index for top leaf layer
    dpai       => mlcanopy_inst%dpai        , &  ! Layer plant area index (m2/m2)
    sumpai     => mlcanopy_inst%sumpai      , &  ! Cumulative plant area index (m2/m2)
    swskyb     => mlcanopy_inst%swskyb      , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd     => mlcanopy_inst%swskyd      , &  ! Atmospheric diffuse solar radiation (W/m2)
    fracsun    => mlcanopy_inst%fracsun     , &  ! Sunlit fraction of canopy layer
    fracsha    => mlcanopy_inst%fracsha     , &  ! Shaded fraction of canopy layer
                                                 ! *** Output ***
    swleaf     => mlcanopy_inst%swleaf      , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swveg      => mlcanopy_inst%swveg       , &  ! Absorbed solar radiation, vegetation (W/m2)
    swvegsun   => mlcanopy_inst%swvegsun    , &  ! Absorbed solar radiation, sunlit canopy (W/m2)
    swvegsha   => mlcanopy_inst%swvegsha    , &  ! Absorbed solar radiation, shaded canopy (W/m2)
    albcan     => mlcanopy_inst%albcan      , &  ! Albedo above canopy
    swsoi      => mlcanopy_inst%swsoi         &  ! Absorbed solar radiation, ground (W/m2)
    )

    !---------------------------------------------------------------------
    ! Multi-layer radiative transfer
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)

          ! Zero out fluxes for all layers

          do ic = 1, ncan(p)
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! Calculate fluxes for leaf layers

          do ic = nbot(p), ntop(p)
         
             ! ild - absorbed diffuse flux per unit leaf area at cumulative LAI, 
             ! average for all leaves (J / m2 leaf / s)

             ild = (1._r8 - albcand(p,ib)) * swskyd(p,ib) * kdm(p,ib) * clump_fac(patch%itype(p)) &
                 * exp(-kdm(p,ib) * sumpai(p,ic) * clump_fac(patch%itype(p)))
            
             ! ilb - absorbed direct beam flux (total with scattering) per unit leaf area 
             ! at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilb = (1._r8 - albcanb(p,ib)) * swskyb(p,ib) * kbm(p,ib) * clump_fac(patch%itype(p)) &
                 * exp(-kbm(p,ib) * sumpai(p,ic) * clump_fac(patch%itype(p)))
            
             ! ilbb - absorbed direct beam flux (unscattered direct component) per unit leaf area 
             ! at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilbb = (1._r8 - omega(p,ib)) * swskyb(p,ib) * kb(p) * clump_fac(patch%itype(p)) &
                  * exp(-kb(p) * sumpai(p,ic) * clump_fac(patch%itype(p)))
            
             ! ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area 
             ! at cumulative LAI, average for all leaves (J / m2 leaf / s)

             ilbs = ilb - ilbb
            
             ! ilsha - total absorbed flux (shaded leaves) per unit shaded leaf area 
             ! at cumulative LAI (J / m2 leaf / s)

             ilsha = ild + ilbs
            
             ! ilsun - total absorbed flux (sunlit leaves) per unit sunlit leaf area 
             ! at cumulative LAI (J / m2 leaf / s)

             ilsun = ilsha + kb(p) * (1._r8 - omega(p,ib)) * swskyb(p,ib)

             ! Solar radiation absorbed by sunlit and shaded leaf

             swleaf(p,ic,isun,ib) = ilsun
             swleaf(p,ic,isha,ib) = ilsha

          end do
       end do
    end do

    !---------------------------------------------------------------------
    ! Canopy summation and soil absorption
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)

          ! Canopy integration, sunlit and shaded leaves

          icsun = 0._r8
          icsha = 0._r8

          do ic = nbot(p), ntop(p)
             icsun = icsun + swleaf(p,ic,isun,ib) * fracsun(p,ic) * dpai(p,ic)
             icsha = icsha + swleaf(p,ic,isha,ib) * fracsha(p,ic) * dpai(p,ic)
          end do

          ! Solar radiation absorbed by vegetation

          swveg(p,ib) = icsun + icsha
          swvegsun(p,ib) = icsun
          swvegsha(p,ib) = icsha

          ! Solar radiation absorbed by ground (soil)

          swsoi(p,ib) = swskyb(p,ib) * (1._r8 - albcanb(p,ib)) * exp(-kbm(p,ib)*(lai(p)+sai(p))*clump_fac(patch%itype(p))) &
                      + swskyd(p,ib) * (1._r8 - albcand(p,ib)) * exp(-kdm(p,ib)*(lai(p)+sai(p))*clump_fac(patch%itype(p)))

          ! Solar radiation reflected by canopy

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumabs = swveg(p,ib) + swsoi(p,ib)
          sumref = suminc - sumabs
!         sumref = swskyb(p,ib) * albcanb(p,ib) + swskyd(p,ib) * albcand(p,ib)

          if (suminc > 0._r8) then
             albcan(p,ib) = sumref / suminc
          else
             albcan(p,ib) = 0._r8
          end if

       end do
    end do

    end associate
  end subroutine GoudriaanRadiation

  !-----------------------------------------------------------------------
  subroutine TwoStreamRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
  omega, kb, avmu, betad, betab, surfalb_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using the two-stream approximation
    !
    ! !USES:
    use clm_varpar, only : numrad, isun, isha
    use clm_varctl, only : iulog
    !
    ! Arguments
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer,  intent(in) :: num_exposedvegp                       ! Number of non-snow-covered veg points in CLM patch filter
    integer,  intent(in) :: filter_exposedvegp(:)                 ! CLM patch filter for non-snow-covered vegetation
    real(r8), intent(in) :: omega(bounds%begp:bounds%endp,numrad) ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: kb(bounds%begp:bounds%endp)           ! Direct beam extinction coefficient
    real(r8), intent(in) :: avmu(bounds%begp:bounds%endp)         ! Average inverse diffuse optical depth per unit leaf area
    real(r8), intent(in) :: betad(bounds%begp:bounds%endp,numrad) ! Upscatter parameter for diffuse radiation
    real(r8), intent(in) :: betab(bounds%begp:bounds%endp,numrad) ! Upscatter parameter for direct beam radiation
    type(surfalb_type), intent(in) :: surfalb_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! Local variable declarations
    integer  :: f                       ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: cp                      ! Column index for CLM g/l/c/p hierarchy
    integer  :: ib                      ! Waveband index
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: b,c,d,h,u,v             ! Intermediate two-stream parameters
    real(r8) :: g1,g2                   ! Intermediate two-stream parameters
    real(r8) :: s1,s2                   ! Intermediate two-stream parameters
    real(r8) :: num1,num2               ! Intermediate two-stream parameters
    real(r8) :: den1,den2               ! Intermediate two-stream parameters
    real(r8) :: n1b,n2b                 ! Two-stream parameters
    real(r8) :: n1d,n2d                 ! Two-stream parameters
    real(r8) :: a1b,a2b                 ! Parameter for sunlit/shaded leaf radiation absorption
    real(r8) :: a1d,a2d                 ! Parameter for sunlit/shaded leaf radiation absorption

    real(r8) :: iupwb0                  ! Direct beam flux scattered upward (reflected) above canopy (W/m2)
    real(r8) :: iupwb                   ! Direct beam flux scattered upward at the canopy depth (W/m2)
    real(r8) :: idwnb                   ! Direct beam flux scattered downward below canopy (W/m2)
    real(r8) :: iabsb                   ! Direct beam flux absorbed by canopy (W/m2)
    real(r8) :: iabsb_sun               ! Direct beam flux absorbed by sunlit canopy (W/m2)
    real(r8) :: iabsb_sha               ! Direct beam flux absorbed by shaded canopy (W/m2)

    real(r8) :: iupwd0                  ! Diffuse flux scattered upward (reflected) above canopy (W/m2)
    real(r8) :: iupwd                   ! Diffuse flux scattered upward at the canopy depth (W/m2)
    real(r8) :: idwnd                   ! Diffuse flux scattered downward below canopy (W/m2)
    real(r8) :: iabsd                   ! Diffuse flux absorbed by canopy (W/m2)
    real(r8) :: iabsd_sun               ! Diffuse flux absorbed by sunlit canopy (W/m2)
    real(r8) :: iabsd_sha               ! Diffuse flux absorbed by shaded canopy (W/m2)

    real(r8) :: diupwb, didwnb          ! Flux derivatives
    real(r8) :: diupwd, didwnd          ! Flux derivatives
    real(r8) :: ilbs, ild               ! Leaf fluxes (per unit leaf area)

    real(r8) :: dir,dif                 ! Direct beam and diffuse fluxes
    real(r8) :: sun,sha                 ! Sunlit and shaded fluxes
    real(r8) :: suminc                  ! Incident radiation for energy conservation check
    real(r8) :: sumref                  ! Reflected radiation for energy conservation check
    real(r8) :: sumabs                  ! Absorbed radiation for energy conservation check
    !---------------------------------------------------------------------

    associate ( &
                                                  ! *** Input ***
    clump_fac => pftcon%clump_fac            , &  ! Foliage clumping index (-)
    ncan      => mlcanopy_inst%ncan          , &  ! Number of aboveground layers
    nbot      => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop      => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    lai       => mlcanopy_inst%lai           , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai           , &  ! Stem area index of canopy (m2/m2)
    sumpai    => mlcanopy_inst%sumpai        , &  ! Cumulative plant area index (m2/m2)
    dpai      => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    swskyb    => mlcanopy_inst%swskyb        , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd    => mlcanopy_inst%swskyd        , &  ! Atmospheric diffuse solar radiation (W/m2)
    fracsun   => mlcanopy_inst%fracsun       , &  ! Sunlit fraction of canopy layer
    fracsha   => mlcanopy_inst%fracsha       , &  ! Shaded fraction of canopy layer
    albsoib   => surfalb_inst%albgrd_col     , &  ! Direct beam albedo of ground (soil)
    albsoid   => surfalb_inst%albgri_col     , &  ! Diffuse albedo of ground (soil)
                                                  ! *** Output ***
    swleaf    => mlcanopy_inst%swleaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swveg     => mlcanopy_inst%swveg         , &  ! Absorbed solar radiation, vegetation (W/m2)
    swvegsun  => mlcanopy_inst%swvegsun      , &  ! Absorbed solar radiation, sunlit canopy (W/m2)
    swvegsha  => mlcanopy_inst%swvegsha      , &  ! Absorbed solar radiation, shaded canopy (W/m2)
    swsoi     => mlcanopy_inst%swsoi         , &  ! Absorbed solar radiation, ground (W/m2)
    albcan    => mlcanopy_inst%albcan          &  ! Albedo above canopy
    )

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          cp = patch%column(p)

          !----------------------------------------------------------------
          ! Canopy fluxes using cumulative lai+sai
          !----------------------------------------------------------------

          ! Common terms

          b = (1._r8 - (1._r8 - betad(p,ib)) * omega(p,ib)) / avmu(p)
          c = betad(p,ib) * omega(p,ib) / avmu(p)
          h = sqrt(b*b - c*c)
          u = (h - b - c) / (2._r8 * h)
          v = (h + b + c) / (2._r8 * h)
          d = omega(p,ib) * kb(p) * swskyb(p,ib) / (h*h - kb(p)*kb(p))
          g1 = (betab(p,ib) * kb(p) - b * betab(p,ib) - c * (1._r8 - betab(p,ib))) * d
          g2 = ((1._r8 - betab(p,ib)) * kb(p) + c * betab(p,ib) + b * (1._r8 - betab(p,ib))) * d
          s1 = exp(-h * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))
          s2 = exp(-kb(p) * (lai(p)+sai(p)) * clump_fac(patch%itype(p)))

          ! Terms for direct beam radiation

          num1 = v * (g1 + g2 * albsoid(cp,ib) + albsoib(cp,ib) * swskyb(p,ib)) * s2
          num2 = g2 * (u + v * albsoid(cp,ib)) * s1
          den1 = v * (v + u * albsoid(cp,ib)) / s1
          den2 = u * (u + v * albsoid(cp,ib)) * s1
          n2b = (num1 - num2) / (den1 - den2)
          n1b = (g2 - n2b * u) / v

          a1b = -g1 *      (1._r8 - s2*s2) / (2._r8 * kb(p)) &
              +  n1b * u * (1._r8 - s2*s1) / (kb(p) + h) + n2b * v * (1._r8 - s2/s1) / (kb(p) - h)
          a2b =  g2 *      (1._r8 - s2*s2) / (2._r8 * kb(p)) &
              -  n1b * v * (1._r8 - s2*s1) / (kb(p) + h) - n2b * u * (1._r8 - s2/s1) / (kb(p) - h)

          ! Direct beam radiative fluxes

          iupwb0 = -g1 + n1b * u + n2b * v
          iupwb = -g1 * s2 + n1b * u * s1 + n2b * v / s1
          idwnb =  g2 * s2 - n1b * v * s1 - n2b * u / s1
          iabsb = swskyb(p,ib) - iupwb0 - (1._r8 - albsoid(cp,ib)) * idwnb - (1._r8 - albsoib(cp,ib)) * swskyb(p,ib) * s2
          iabsb_sun = (1._r8 - omega(p,ib)) * ((1._r8 - s2) * swskyb(p,ib) + 1._r8 / avmu(p) * (a1b + a2b) * clump_fac(patch%itype(p)))
          iabsb_sha = iabsb - iabsb_sun

          ! Terms for diffuse radiation
 
          num1 = swskyd(p,ib) * (u + v * albsoid(cp,ib)) * s1
          den1 = v * (v + u * albsoid(cp,ib)) / s1
          den2 = u * (u + v * albsoid(cp,ib)) * s1
          n2d = num1 / (den1 - den2)
          n1d = -(swskyd(p,ib) + n2d * u) / v

          a1d =  n1d * u * (1._r8 - s2*s1) / (kb(p) + h) + n2d * v * (1._r8 - s2/s1) / (kb(p) - h)
          a2d = -n1d * v * (1._r8 - s2*s1) / (kb(p) + h) - n2d * u * (1._r8 - s2/s1) / (kb(p) - h)

          ! Diffuse radiative fluxes

          iupwd0 = n1d * u + n2d * v
          iupwd =  n1d * u * s1 + n2d * v / s1
          idwnd = -n1d * v * s1 - n2d * u / s1
          iabsd = swskyd(p,ib) - iupwd0 - (1._r8 - albsoid(cp,ib)) * idwnd
          iabsd_sun = (1._r8 - omega(p,ib)) / avmu(p) * (a1d + a2d) * clump_fac(patch%itype(p))
          iabsd_sha = iabsd - iabsd_sun

          !----------------------------------------------------------------
          ! Save necessary radiative fluxes
          !----------------------------------------------------------------

          ! Albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = iupwb0 + iupwd0
          if (suminc > 0._r8) then
             albcan(p,ib) = sumref / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Solar radiation absorbed by canopy

          swveg(p,ib) = iabsb +  iabsd
          swvegsun(p,ib) = iabsb_sun + iabsd_sun
          swvegsha(p,ib) = iabsb_sha + iabsd_sha

          ! Solar radiation absorbed by ground (soil)

          dir = swskyb(p,ib) * s2 * (1._r8 - albsoib(cp,ib))
          dif = (idwnb + idwnd) * (1._r8 - albsoid(cp,ib))
          swsoi(p,ib) = dir + dif

          ! Conservation check: total incident = total reflected + total absorbed

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = iupwb0 + iupwd0
          sumabs = swveg(p,ib) + swsoi(p,ib)

          if (abs(suminc - (sumabs+sumref)) >= 1.e-06_r8) then
             call endrun (msg='ERROR: TwoStreamRadiation: total solar radiation conservation error')
          end if

          !----------------------------------------------------------------
          ! Repeat two-stream calculations for each leaf layer to
          ! calculate leaf fluxes
          !----------------------------------------------------------------

          ! Zero out fluxes for all layers

          do ic = 1, ncan(p)
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! Calculate fluxes for leaf layers

          do ic = nbot(p), ntop(p)

             ! s1 and s2 depend on cumulative plant area index

             s1 = exp(-h * sumpai(p,ic) * clump_fac(patch%itype(p)))
             s2 = exp(-kb(p) * sumpai(p,ic) * clump_fac(patch%itype(p)))

             ! ilbs - absorbed direct beam flux (scattered direct component) per unit leaf area
             ! at cumulative LAI, average for all leaves (J / m2 leaf / s)

             diupwb =  kb(p) * g1 * s2 - h * n1b * u * s1 + h * n2b * v / s1
             didwnb = -kb(p) * g2 * s2 + h * n1b * v * s1 - h * n2b * u / s1
             ilbs = (omega(p,ib) * kb(p) * swskyb(p,ib) * s2 + (diupwb - didwnb)) * clump_fac(patch%itype(p))

             ! ild - absorbed diffuse flux per unit leaf area at cumulative LAI,
             ! average for all leaves (J / m2 leaf / s)

             diupwd = -h * n1d * u * s1 + h * n2d * v / s1
             didwnd =  h * n1d * v * s1 - h * n2d * u / s1
             ild = (diupwd - didwnd) * clump_fac(patch%itype(p))

             ! Save leaf fluxes per unit sunlit and shaded leaf area (W/m2 leaf)

             swleaf(p,ic,isun,ib) = (1._r8 - omega(p,ib)) * kb(p) * swskyb(p,ib) + (ilbs + ild)
             swleaf(p,ic,isha,ib) =  ilbs + ild

          end do

       end do      ! end grid point loop
    end do         ! end waveband loop

    !---------------------------------------------------------------------
    ! Adjust leaf fluxes as needed. The sum of the fluxes for sunlit and
    ! shaded leaves should equal the total absorbed by the canopy, but may
    ! not because of inaccuracies in the flux derivatives (this is a small
    ! error if the dpai increment is small). Normalize these fluxes to sum
    ! to the canopy absorption.
    !---------------------------------------------------------------------

    do ib = 1, numrad
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)

          ! Sum canopy absorption (W/m2 ground) using leaf fluxes per unit sunlit
          ! and shaded leaf area (W/m2 leaf)

          sumabs = 0._r8
          do ic = nbot(p), ntop(p)
             sun = swleaf(p,ic,isun,ib) * fracsun(p,ic) * dpai(p,ic)
             sha = swleaf(p,ic,isha,ib) * fracsha(p,ic) * dpai(p,ic)
             sumabs = sumabs + sun + sha
          end do

          ! Normalize profile

          if (sumabs > 0.0) then
             do ic = nbot(p), ntop(p)
                swleaf(p,ic,isun,ib) = swleaf(p,ic,isun,ib) * swveg(p,ib) / sumabs
                swleaf(p,ic,isha,ib) = swleaf(p,ic,isha,ib) * swveg(p,ib) / sumabs
             end do
          end if

       end do
    end do

    end associate
  end subroutine TwoStreamRadiation

end module SolarRadiationMod

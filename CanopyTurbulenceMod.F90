module CanopyTurbulenceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy turbulence, aeorodynamic conductances, and scalar profiles using above-
  ! and within-canopy coupling with a roughness sublayer (RSL) parameterization 
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use decompMod , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  real(r8), parameter :: zetam = -1.574_r8 ! CLM: transition point of flux-gradient relation (wind profile)
  real(r8), parameter :: zetah = -0.465_r8 ! CLM: transition point of flux-gradient relation (temperature profile)
  real(r8), parameter :: zetamaxstable = 2.0_r8 ! CLM: Maximum value for zeta under stable conditions
! real(r8), parameter :: zetamaxstable = 0.5_r8 ! CLM: Maximum value for zeta under stable conditions

  real(r8), parameter :: cd = 0.25_r8               ! RSL leaf drag coefficient (dimensionless)
  real(r8), parameter :: beta_neutral_max = 0.35_r8 ! RSL maximum value for beta in neutral conditions
  real(r8), parameter :: cr = 0.3_r8                ! RSL parameter to calculate beta_neutral
  real(r8), parameter :: c2 = 0.5_r8                ! RSL depth scale multiplier
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyTurbulenceDummy  ! Canopy scalar profiles equal reference height values
  public :: CanopyTurbulence       ! Canopy turbulence and scalar profiles used RSL theory
  public :: LookupPsihatINI        ! Initialize the RSL psihat look-up tables
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: phim_monin_obukhov    ! Monin-Obukhov phi stability function for momentum
  private :: phic_monin_obukhov    ! Monin-Obukhov phi stability function for scalars
  private :: psim_monin_obukhov    ! Monin-Obukhov psi stability function for momentum
  private :: psic_monin_obukhov    ! Monin-Obukhov psi stability function for scalars
  private :: ObuFunc               ! Function to solve for the Obukhov length
  private :: GetBeta               ! Calculate beta = u* / u(h)
  private :: GetPrSc               ! Turbulent Prandlt number (Pr) and Schmidt number (Sc) at canopy top
  private :: GetPsiRSL             ! RSL-modified stability functions
  private :: LookupPsihat          ! Determines the RSL function psihat as provided through a look-up table
  private :: ScalarProfile         ! Calculate scalar profile using implicit solution
  private :: ScalarProfileIterative! Calculate scalar profile using iterative solution
  private :: ObuFuncCLM            ! Calculate the Obukhov length as in CLM
  private :: MoninObukIniCLM       ! CLM initialization of the Obukhov length
  private :: FrictionVelocityCLM   ! CLM friction velocity and stability dependent log-z functions
  private :: GetPsiCLM             ! CLM surface-layer stability functions
  private :: StabilityFunc1        ! CLM Monin-Obukhov stability function for momentum
  private :: StabilityFunc2        ! CLM Monin-Obukhov stability function for scalars

contains

  !-----------------------------------------------------------------------
  function phim_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for momentum
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       phi = 1._r8 / sqrt(sqrt(1._r8 - 16._r8 * zeta))
    else                             ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phim_monin_obukhov

  !-----------------------------------------------------------------------
  function phic_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       phi = 1._r8 / sqrt(1._r8 - 16._r8 * zeta)
    else                             ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phic_monin_obukhov

  !-----------------------------------------------------------------------
  function psim_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for momentum
    !
    ! !USES:
    use clm_varcon, only : pi => rpi
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x)/2._r8) + log((1._r8+x*x)/2._r8) - 2._r8*atan(x) + pi * 0.5_r8
    else                             ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psim_monin_obukhov

  !-----------------------------------------------------------------------
  function psic_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then           ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x*x)/2._r8)
    else                             ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psic_monin_obukhov

  !-----------------------------------------------------------------------
  subroutine CanopyTurbulenceDummy (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Dummy canopy turbulence - canopy scalar profiles equal reference height values
    !
    ! !USES:
    use clm_varpar, only : isun, isha
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p       ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                 ! Aboveground layer index
    !---------------------------------------------------------------------

    associate ( &
                                              ! *** Input ***
    ncan     => mlcanopy_inst%ncan       , &  ! Number of aboveground layers
    dpai     => mlcanopy_inst%dpai       , &  ! Layer plant area index (m2/m2)
    uref     => mlcanopy_inst%uref       , &  ! Wind speed at reference height (m/s)
    tref     => mlcanopy_inst%tref       , &  ! Air temperature at reference height (K)
    eref     => mlcanopy_inst%eref       , &  ! Vapor pressure at reference height (Pa)
    co2ref   => mlcanopy_inst%co2ref     , &  ! Atmospheric CO2 at reference height (umol/mol)
    tleaf    => mlcanopy_inst%tleaf      , &  ! Leaf temperature (K)
                                              ! *** Output ***
    z0mg     => mlcanopy_inst%z0mg       , &  ! Roughness length of ground (m)
    Lc       => mlcanopy_inst%Lc         , &  ! Canopy density length scale (m)
    obu      => mlcanopy_inst%obu        , &  ! Obukhov length (m)
    obuold   => mlcanopy_inst%obuold     , &  ! Obukhov length from previous iteration
    nmozsgn  => mlcanopy_inst%nmozsgn    , &  ! Number of times stability changes sign during iteration
    ustar    => mlcanopy_inst%ustar      , &  ! Friction velocity (m/s)
    tstar    => mlcanopy_inst%tstar      , &  ! Temperature scale (K)
    qstar    => mlcanopy_inst%qstar      , &  ! Water vapor scale (kg/kg)
    PrSc     => mlcanopy_inst%PrSc       , &  ! Prandtl (Schmidt) number at canopy top
    zdisp    => mlcanopy_inst%zdisp      , &  ! Displacement height (m)
    uaf      => mlcanopy_inst%uaf        , &  ! Wind speed at canopy top (m/s)
    taf      => mlcanopy_inst%taf        , &  ! Air temperature at canopy top (K)
    eaf      => mlcanopy_inst%eaf        , &  ! Vapor pressure at canopy top (Pa)
    qaf      => mlcanopy_inst%qaf        , &  ! Specific humidity at canopy top (kg/kg)
    gah      => mlcanopy_inst%gah        , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    shair    => mlcanopy_inst%shair      , &  ! Canopy air sensible heat flux (W/m2)
    etair    => mlcanopy_inst%etair      , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    stair    => mlcanopy_inst%stair      , &  ! Canopy air storage heat flux (W/m2)
    wind     => mlcanopy_inst%wind       , &  ! Wind speed profile (m/s)
    tair     => mlcanopy_inst%tair       , &  ! Air temperature profile (K)
    eair     => mlcanopy_inst%eair       , &  ! Vapor pressure profile (Pa)
    cair     => mlcanopy_inst%cair       , &  ! Atmospheric CO2 profile (umol/mol)
    tveg     => mlcanopy_inst%tveg       , &  ! Vegetation temperature profile (K)
    ga_prof  => mlcanopy_inst%ga_prof    , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    wind_most=> mlcanopy_inst%wind_most  , &  ! Wind speed profile from MOST (m/s)
    tair_most=> mlcanopy_inst%tair_most    &  ! Air temperature profile from MOST (K)
    )

    z0mg(p) = 0._r8
    Lc(p) = 0._r8
    obu(p) = 0._r8
    obuold(p) = 0._r8
    nmozsgn(p) = 0
    ustar(p) = uref(p)     ! cannot be zero for analysis package to work
    tstar(p) = 0._r8
    qstar(p) = 0._r8
    PrSc(p) = 0._r8
    zdisp(p) = 0._r8
    uaf(p) = 0._r8
    taf(p) = 0._r8
    eaf(p) = 0._r8
    qaf(p) = 0._r8
    gah(p) = 0._r8

    do ic = 0, ncan(p)
       wind(p,ic) = uref(p)
       tair(p,ic) = tref(p)
       eair(p,ic) = eref(p)
       cair(p,ic) = co2ref(p)
       ga_prof(p,ic) = (1._r8 / 10._r8) * 42.3_r8  ! Small non-zero resistance
       wind_most(p,ic) = uref(p)
       tair_most(p,ic) = tref(p)
    end do

    ! For soil surface, use a large resistance to approximate bare ground and
    ! so that soil fluxes are small

    ga_prof(p,0) = (1._r8 / 100._r8) * 42.3_r8

    tveg(p,0,isun) = tref(p)
    tveg(p,0,isha) = tref(p)
    do ic = 1, ncan(p)
       shair(p,ic) = 0._r8
       etair(p,ic) = 0._r8
       stair(p,ic) = 0._r8
       if (dpai(p,ic) > 0._r8) then
          tveg(p,ic,isun) = tleaf(p,ic,isun)
          tveg(p,ic,isha) = tleaf(p,ic,isha)
       else
          tveg(p,ic,isun) = tref(p)
          tveg(p,ic,isha) = tref(p)
       end if
    end do


    end associate
  end subroutine CanopyTurbulenceDummy

  !-----------------------------------------------------------------------
  subroutine CanopyTurbulence (p, niter, soilstate_inst, temperature_inst, &
  energyflux_inst, waterflux_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy turbulence, aeorodynamic conductances, and wind/temperature/water vapor
    ! profiles using above- and within-canopy coupling with a roughness sublayer
    ! (RSL) parameterization 
    !
    ! !USES:
    use clm_varctl, only : turb_type, use_scalars
    use clm_varcon, only : mmh2o, mmdry, vkc
    use MathToolsMod, only : hybrid
    use PatchType, only : patch
    use SoilStateType, only : soilstate_type
    use TemperatureType, only : temperature_type
    use EnergyFluxType, only : energyflux_type
    use WaterFluxType, only : waterflux_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    use pftconMod, only : pftcon
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p       ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: niter   ! Iteration index
    type(soilstate_type), intent(inout) :: soilstate_inst
    type(temperature_type), intent(inout) :: temperature_inst
    type(energyflux_type), intent(inout) :: energyflux_inst
    type(waterflux_type), intent(inout) :: waterflux_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c                  ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                 ! Canopy layer index
    integer  :: il                 ! Sunlit (1) or shaded (2) leaf index
    integer  :: iterate            ! CLM: 1 = iterate for Obukhov length. 0 = use specified Obukhov length
    real(r8) :: dummy              ! Dummy argument
    real(r8) :: lm                 ! Length scale for momentum (m)
    real(r8) :: beta               ! Value of u*/u(h) at canopy top
    real(r8) :: res                ! Resistance (s/m)
    real(r8) :: zu, zl             ! Upper and lower heights for within canopy resistances (m)
    real(r8) :: obu0, obu1         ! Initial estimates for Obukhov length (m)
    real(r8) :: z0hg               ! Roughness length of ground for sensible heat (m)
    real(r8) :: zlog_m, zlog_c     ! Log-profiles at ground
    real(r8) :: ustar_g            ! Friction velocity at ground (m/s)
    real(r8) :: eta                ! Used to limit value for beta/lm
    real(r8) :: beta_over_lm       ! beta/lm (1/m)
    real(r8) :: tol                ! Accuracy tolerance for Obukhov length (m)
    real(r8) :: zeta               ! (z - d)/L
    real(r8) :: phic               ! Flux-gradient relationship at (hc - d)/L
    real(r8) :: kc_at_hc           ! Scalar eddy diffusivity at canopy top (m2/s)

!!! NEW
    real(r8) :: dt                 ! Height below canopy top (dt = ztop - zdisp)
    real(r8) :: psim               ! psi function for momentum
    real(r8) :: psic               ! psi function for scalars
    real(r8) :: psim1, psim2
    real(r8) :: psic1, psic2
    real(r8) :: ga_add
    real(r8) :: sumres

    real(r8) :: zdisp_most, z0m_most, zeta_most, psim_most_z0m, psim_most, zlog_most, ustar_most
    real(r8) :: z0c_most, sh_most, tstar_most, psic_most_z0c, psic_most, ts_most
    !---------------------------------------------------------------------

    associate ( &
                                              ! *** Input ***
    ztop     => mlcanopy_inst%ztop       , &  ! Canopy height (m)
    ncan     => mlcanopy_inst%ncan       , &  ! Number of aboveground layers
    ntop     => mlcanopy_inst%ntop       , &  ! Index for top leaf layer
    lai      => mlcanopy_inst%lai        , &  ! Leaf area index of canopy (m2/m2)
    sai      => mlcanopy_inst%sai        , &  ! Stem area index of canopy (m2/m2)
    zs       => mlcanopy_inst%zs         , &  ! Canopy height for scalar concentration and source (m)
    co2ref   => mlcanopy_inst%co2ref     , &  ! Atmospheric CO2 at reference height (umol/mol)
    pref     => mlcanopy_inst%pref       , &  ! Air pressure at reference height (Pa)
    rhomol   => mlcanopy_inst%rhomol     , &  ! Molar density at reference height (mol/m3)
    eg       => mlcanopy_inst%eg         , &  ! Soil surface vapor pressure (Pa)
    tg       => mlcanopy_inst%tg         , &  ! Soil surface temperature (K)
    zref     => mlcanopy_inst%zref       , &  ! Reference height (m)
    uref     => mlcanopy_inst%uref       , &  ! Wind speed at reference height (m/s)
    thref    => mlcanopy_inst%thref      , &  ! Atmospheric potential temperature (K)
    cpair    => mlcanopy_inst%cpair      , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    z0mr     => pftcon%z0mr              , &  ! Ratio of momentum roughness length to canopy top height
    displar  => pftcon%displar           , &  ! Ratio of displacement height to canopy top height
                                              ! *** Output ***
    z0mg     => mlcanopy_inst%z0mg       , &  ! Roughness length of ground (m)
    Lc       => mlcanopy_inst%Lc         , &  ! Canopy density length scale (m)
    obu      => mlcanopy_inst%obu        , &  ! Obukhov length (m)
    obu_gah  => mlcanopy_inst%obu_gah    , &  ! Obukhov length used for gah (m)
    obuold   => mlcanopy_inst%obuold     , &  ! Obukhov length from previous iteration
    nmozsgn  => mlcanopy_inst%nmozsgn    , &  ! Number of times stability changes sign during iteration
    ustar    => mlcanopy_inst%ustar      , &  ! Friction velocity (m/s)
    tstar    => mlcanopy_inst%tstar      , &  ! Temperature scale (K)
    qstar    => mlcanopy_inst%qstar      , &  ! Water vapor scale (kg/kg)
    PrSc     => mlcanopy_inst%PrSc       , &  ! Prandtl (Schmidt) number at canopy top
    zdisp    => mlcanopy_inst%zdisp      , &  ! Displacement height (m)
    uaf      => mlcanopy_inst%uaf        , &  ! Wind speed at canopy top (m/s)
    taf      => mlcanopy_inst%taf        , &  ! Air temperature at canopy top (K)
    eaf      => mlcanopy_inst%eaf        , &  ! Vapor pressure at canopy top (Pa)
    qaf      => mlcanopy_inst%qaf        , &  ! Specific humidity at canopy top (kg/kg)
    wind     => mlcanopy_inst%wind       , &  ! Wind speed profile (m/s)
    wind_most=> mlcanopy_inst%wind_most  , &  ! Wind speed profile from MOST (m/s)
    tair     => mlcanopy_inst%tair       , &  ! Air temperature profile (K)
    tair_most=> mlcanopy_inst%tair_most  , &  ! Air temperature profile from MOST (K)
    eair     => mlcanopy_inst%eair       , &  ! Vapor pressure profile (Pa)
    cair     => mlcanopy_inst%cair       , &  ! Atmospheric CO2 profile (umol/mol)
    gah      => mlcanopy_inst%gah        , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    ga_prof  => mlcanopy_inst%ga_prof      &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    )

    c = patch%column(p)

    ! Initialization for first iteration

    if (niter == 1) then
       obuold(p) = 0._r8
       nmozsgn(p) = 0
    end if

    ! Set roughness length of ground

    z0mg(p) = 0.01_r8

    ! Canopy density length scale

    Lc(p) = ztop(p) / (cd * (lai(p) + sai(p)))

    ! These are not used, but are needed to pass into the hybrid root solver

    ic = 0
    il = 0

    ! Calculate Obukhov length (obu)

    select case (turb_type)

       case (1:3)

       ! Calculate Obukhov length (obu) using the function ObuFunc to iterate until
       ! the change in obu is < tol. Do not use the returned value (dummy), and
       ! instead use the value of obu calculated in the final call to ObuFunc.

       obu0 = 100._r8          ! Initial estimate for Obukhov length (m)
       obu1 = -100._r8         ! Initial estimate for Obukhov length (m)
       tol = 0.1_r8            ! Accuracy tolerance for Obukhov length (m)

       dummy = hybrid ('CanopyTurbulence', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)
       obu(p) = obu_gah(p)

       case (4)

       ! Use CLM calculation of Obukhov length (obu)

       iterate = 1
       call ObuFuncCLM (p, iterate, mlcanopy_inst)
       obu(p) = obu_gah(p)

    end select

    ! Check to see if Obukhov length is changing signs between iterations.
    ! If too many changes in sign, set it to a near-neutral value.

    if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p) + 1
    obuold(p) = obu(p)
    if (nmozsgn(p) >= 4) then
       obu(p) = -1000._r8
       select case (turb_type)
          case (1:3)
             dummy = ObuFunc (p, ic, il, mlcanopy_inst, obu(p))
             obu(p) = obu_gah(p)
          case (4)
             iterate = 0
             call ObuFuncCLM (p, iterate, mlcanopy_inst)
             obu(p) = obu_gah(p)
       end select
    end if

    ! Above-canopy wind profile: wind speed is defined at zs

    select case (turb_type)
       case (1:3)

       dt = ztop(p) - zdisp(p)
       beta = ustar(p) / uaf(p)

       do ic = ntop(p)+1, ncan(p)
          zeta = (zs(p,ic) - zdisp(p)) / obu(p)
          call GetPsiRSL (zs(p,ic), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim, psic)
          zlog_m = log((zs(p,ic)-zdisp(p)) / (ztop(p)-zdisp(p)))
          wind(p,ic) = ustar(p) / vkc * (zlog_m + psim)
       end do

    end select

    ! Above-canopy aerodynamic conductances: these are defined between
    ! zs(i) and zs(i+1). Note the special case for the top layer to the
    ! reference height.

    select case (turb_type)
       case (1:3)

       do ic = ntop(p)+1, ncan(p)-1
          zeta = (zs(p,ic) - zdisp(p)) / obu(p)
          call GetPsiRSL (zs(p,ic), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim1, psic1)
          zeta = (zs(p,ic+1) - zdisp(p)) / obu(p)
          call GetPsiRSL (zs(p,ic+1), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim2, psic2)
          ! equivalent to: -psi_c_z2 + psi_c_z1 + psi_c_rsl_z2 - psi_c_rsl_z1
          psic = psic2 - psic1
          zlog_c = log((zs(p,ic+1)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
          ga_prof(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)
       end do

       ic = ncan(p)
       zeta = (zs(p,ic) - zdisp(p)) / obu(p)
       call GetPsiRSL (zs(p,ic), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim1, psic1)
       zeta = (zref(p) - zdisp(p)) / obu(p)
       call GetPsiRSL (zref(p), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim2, psic2)
       psic = psic2 - psic1
       zlog_c = log((zref(p)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
       ga_prof(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)

    end select

    ! Within-canopy wind profile: wind speed is defined at zs

    beta = ustar(p) / uaf(p)
    lm = 2._r8 * beta**3 * Lc(p)

    select case (turb_type)
       case (1, 2)
!         eta = min (beta/lm*ztop(p), 10._r8)
          eta = min (beta/lm*ztop(p), 20._r8)
       case (3, 4)
          eta = 3._r8
    end select

    beta_over_lm = eta / ztop(p)
!   beta_over_lm = beta / lm

    do ic = 1, ntop(p)
       wind(p,ic) = uaf(p) * exp((zs(p,ic) - ztop(p)) * beta_over_lm)
       wind(p,ic) = max(wind(p,ic), 0.1_r8)
    end do

    ! Scalar eddy diffusivity at canopy top

    select case (turb_type)

       case (1:3)

       kc_at_hc = 2._r8 * beta**3 * Lc(p) * ustar(p) / PrSc(p)

       case (4)

       zeta = (ztop(p) - zdisp(p)) / obu(p)
       if (zeta < zetah) then             ! very unstable (zeta < zetah)
          phic = 0.9 * vkc**1.333_r8 * (-zeta)**(-0.333_r8)
       else if (zeta < 0._r8) then        ! unstable (zetah <= zeta < 0)
          phic = 1._r8 / sqrt(1._r8 - 16._r8 * zeta)
       else if (zeta <=  1._r8) then      ! stable (0 <= zeta <= 1)
          phic = 1._r8 + 5._r8 * zeta
       else                               ! very stable (zeta > 1)
          phic = 5._r8 + zeta
       end if
       kc_at_hc = vkc * ustar(p) * (ztop(p) - zdisp(p)) / phic

    end select

    ! Within-canopy aerodynamic conductances: these are defined between
    ! zs(i) and zs(i+1). Note the special case for the top layer to the
    ! canopy height.

    do ic = 1, ntop(p)-1
       zl = zs(p,ic) - ztop(p)
       zu = zs(p,ic+1) - ztop(p)
!      res = PrSc(p) / (beta * ustar(p)) * (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
       res = 1._r8 / (kc_at_hc * beta_over_lm) * (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
       ga_prof(p,ic) = rhomol(p) / res
    end do

    ! Special case for top layer: conductance to top of canopy

    ic = ntop(p)
    zl = zs(p,ic) - ztop(p)
    zu = ztop(p) - ztop(p)
!   res = PrSc(p) / (beta * ustar(p)) * (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
    res = 1._r8 / (kc_at_hc * beta_over_lm) * (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
    ga_prof(p,ic) = rhomol(p) / res

    ! Now add additional resistance to first atmospheric layer

    if (ncan(p) == ntop(p)) then
       ga_add = gah(p)
    else if (ncan(p) > ntop(p)) then
       select case (turb_type)
          case (1:3)
          zeta = (ztop(p) - zdisp(p)) / obu(p)
          call GetPsiRSL (ztop(p), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim1, psic1)
          zeta = (zs(p,ic+1) - zdisp(p)) / obu(p)
          call GetPsiRSL (zs(p,ic+1), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim2, psic2)
          psic = psic2 - psic1
          zlog_c = log((zs(p,ic+1)-zdisp(p)) / (ztop(p)-zdisp(p)))
          ga_add = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)
       end select
    end if
    ga_prof(p,ic) = 1._r8 / (1._r8 / ga_prof(p,ic) + 1._r8 / ga_add)

    ! Make sure above-canopy aerodynamic resistances sum to 1/gah

    sumres = 0._r8
    do ic = ntop(p)+1, ncan(p)
       sumres = sumres + 1._r8 / ga_prof(p,ic)
    end do
    sumres = sumres + 1._r8 / ga_add

    if (abs(1._r8/sumres - gah(p)) > 1.e-06_r8) then
       call endrun (msg=' ERROR: CanopyTurbulenceMod: above-canopy aerodynamic conductance error')
    end if

    ! Resistance at ground
    
!   zl = 0._r8 - ztop(p)
!   zu = zs(p,1) - ztop(p)
!   res = PrSc(p) / (beta * ustar(p)) * (exp(-zl*beta/lm) - exp(-zu*beta/lm))
!   res = ga_prof(p,1)

    z0hg = 0.1_r8 * z0mg(p)
    zlog_m = log(zs(p,1)/z0mg(p))
    zlog_c = log(zs(p,1)/z0hg)
    ustar_g = wind(p,1) * vkc / zlog_m
    ustar_g = max(ustar_g, 0.01_r8)
    res = zlog_c / (vkc * ustar_g)
    ga_prof(p,0) = rhomol(p) / res

    if (zlog_m < 0._r8 .or. zlog_c < 0._r8) then
       call endrun (msg=' ERROR: CanopyTurbulenceMod: soil roughness error')
    end if

    ! Convert resistance (s/m) to conductance (mol/m2/s)
    ! Limit resistances to < 500 s/m

    do ic = 0, ncan(p)
       res = min (rhomol(p)/ga_prof(p,ic), 500._r8)
       ga_prof(p,ic) = rhomol(p) / res
    end do

    ! Values at ground

    wind(p,0) = 0._r8
    tair(p,0) = tg(p)
    eair(p,0) = eg(p)

    ! Calculate within-canopy scalar profiles for temperature and vapor pressure
    ! using an implicit solution

     call ScalarProfile (p, niter, soilstate_inst, temperature_inst, &
     energyflux_inst, waterflux_inst, mlcanopy_inst)

    ! Calculate within-canopy scalar profiles for temperature, vapor pressure, and CO2
    ! using an iterative coupling
    !
    ! Temperature: scalar is T*cp (J/mol) and flux is J/m2/s
    ! Water vapor: scalar is e/P (mol/mol) and flux is mol/m2/s
    ! CO2: scalar is c (umol/mol) and flux is umol/m2/s

    ! call ScalarProfileIterative (ncan(p), tref(p), ga_prof(p,:), sh_prof(p,:), &
    ! 1._r8/cpair(p), tair_old(p,:), tair(p,:))

    ! call ScalarProfileIterative (ncan(p), eref(p), ga_prof(p,:), et_prof(p,:), &
    ! pref(p), eair_old(p,:), eair(p,:))

    ! call ScalarProfileIterative (ncan(p), co2ref(p), ga_prof(p,:), fc_prof(p,:), &
    ! 1._r8, cair_old(p,:), cair(p,:))

    do ic = 0, ncan(p)
       cair(p,ic) = co2ref(p)
    end do

    if (.not. use_scalars) then
       do ic = 1, ntop(p)
          wind(p,ic) = wind(p,ntop(p))
          tair(p,ic) = tair(p,ntop(p))
          eair(p,ic) = eair(p,ntop(p))
       end do
    end if

    taf(p) = tair(p,ntop(p))
    eaf(p) = eair(p,ntop(p))
    qaf(p) = mmh2o / mmdry * eaf(p) / (pref(p) - (1._r8 - mmh2o/mmdry) * eaf(p))

    ! Wind and temperature profiles using MOST

    wind_most(p,0) = -999._r8
    tair_most(p,0) = -999._r8

    zdisp_most = displar(patch%itype(p)) * ztop(p)
    z0m_most = z0mr(patch%itype(p)) * ztop(p)
    z0c_most = z0m_most

    zeta_most = (zref(p) - zdisp_most) / obu(p)
    psim_most_z0m = psim_monin_obukhov(z0m_most/obu(p))
    psim_most = -psim_monin_obukhov(zeta_most) + psim_most_z0m
    zlog_most = log((zref(p)-zdisp_most) / z0m_most)
    ustar_most = uref(p) * vkc / (zlog_most + psim_most)

    sh_most = -cpair(p) * (tair(p,ntop(p)+1) - tair(p,ntop(p))) * ga_prof(p,ntop(p))
    tstar_most = -sh_most / (rhomol(p) * cpair(p) * ustar_most)
    psic_most_z0c = psic_monin_obukhov(z0c_most/obu(p))
    psic_most = -psic_monin_obukhov(zeta_most) + psic_most_z0c
    zlog_most = log((zref(p)-zdisp_most) / z0c_most)
    ts_most = thref(p) - tstar_most / vkc * (zlog_most + psic_most)

    do ic = ncan(p), 1, -1
       if (zs(p,ic) > zdisp_most) then
          zeta_most = (zs(p,ic) - zdisp_most) / obu(p)
          psim_most = -psim_monin_obukhov(zeta_most) + psim_most_z0m
          psic_most = -psic_monin_obukhov(zeta_most) + psic_most_z0c
          zlog_most = log((zs(p,ic)-zdisp_most) / z0m_most)
          wind_most(p,ic) = ustar_most / vkc * (zlog_most + psim_most)
          zlog_most = log((zs(p,ic)-zdisp_most) / z0c_most)
          tair_most(p,ic) = ts_most + tstar_most / vkc * (zlog_most + psic_most)
       else
          wind_most(p,ic) = -999._r8
          tair_most(p,ic) = -999._r8
       end if
    end do

    end associate
  end subroutine CanopyTurbulence

  !-----------------------------------------------------------------------
  function ObuFunc (p, ic, il, mlcanopy_inst, obu_val) result(obu_dif)
    !
    ! !DESCRIPTION:
    ! This is the function to solve for the Obukhov length (obu). For the current
    ! estimate of the Obukhov length (obu_val), calculate u* , T*, and q* and then
    ! the new length (obu). The function value is the change in Obukhov length and
    ! equals zero when the the Obukhov length does not change value between iterations.
    !
    ! !USES:
    use clm_varcon, only : grav, vkc
    use clm_varctl, only : turb_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p               ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic              ! Aboveground layer index
    integer, intent(in)  :: il              ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: obu_val         ! Input value for Obukhov length (m)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: c1                          ! Parameter to calculate beta_neutral
    real(r8) :: beta_neutral                ! Neutral value for beta = u*/u(h), ~0.3-0.5
    real(r8) :: beta                        ! Value of u*/u(h) at canopy top
    real(r8) :: dt                          ! Height below canopy top (dt = ztop - zdisp)
    real(r8) :: zeta                        ! Monin-Obukhov stability parameter (z-d)/L
    real(r8) :: psim                        ! psi function for momentum
    real(r8) :: psic                        ! psi function for scalars
    real(r8) :: zlog                        ! log((zref-zdisp)/(ztop-zdisp))
    real(r8) :: tvstar                      ! Virtual potential temperature scale (K)
    real(r8) :: obu_dif                     ! Difference in Obuhkov length (m)
    !---------------------------------------------------------------------

    associate ( &
                                            ! *** Input ***
    ztop    => mlcanopy_inst%ztop      , &  ! Canopy height (m)
    z0mg    => mlcanopy_inst%z0mg      , &  ! Roughness length of ground (m)
    lai     => mlcanopy_inst%lai       , &  ! Leaf area index of canopy (m2/m2)
    sai     => mlcanopy_inst%sai       , &  ! Stem area index of canopy (m2/m2)
    Lc      => mlcanopy_inst%Lc        , &  ! Canopy density length scale (m)
    zref    => mlcanopy_inst%zref      , &  ! Reference height (m)
    uref    => mlcanopy_inst%uref      , &  ! Wind speed at reference height (m/s)
    rhomol  => mlcanopy_inst%rhomol    , &  ! Molar density at reference height (mol/m3)
    qref    => mlcanopy_inst%qref      , &  ! Specific humidity at reference height (kg/kg)
    qaf      => mlcanopy_inst%qaf      , &  ! Specific humidity at canopy top (kg/kg)
    thref   => mlcanopy_inst%thref     , &  ! Atmospheric potential temperature (K)
    thvref  => mlcanopy_inst%thvref    , &  ! Atmospheric virtual potential temperature (K)
    taf     => mlcanopy_inst%taf       , &  ! Air temperature at canopy top (K)
                                            ! *** Output ***
    ustar   => mlcanopy_inst%ustar     , &  ! Friction velocity (m/s)
    tstar   => mlcanopy_inst%tstar     , &  ! Temperature scale (K)
    qstar   => mlcanopy_inst%qstar     , &  ! Water vapor scale (kg/kg)
    gah     => mlcanopy_inst%gah       , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu     => mlcanopy_inst%obu       , &  ! Obukhov length (m)
    obu_gah => mlcanopy_inst%obu_gah   , &  ! Obukhov length used for gah (m)
    PrSc    => mlcanopy_inst%PrSc      , &  ! Prandtl (Schmidt) number at canopy top
    zdisp   => mlcanopy_inst%zdisp     , &  ! Displacement height (m)
    uaf     => mlcanopy_inst%uaf         &  ! Wind speed at canopy top (m/s)
    )

    ! Use this current value of Obukhov length

    obu(p) = obu_val

    ! Prevent near-zero value of Obukhov length

    if (abs(obu(p)) <= 0.1_r8) obu(p) = 0.1_r8

    ! Determine neutral value of beta

    select case (turb_type)
       case (1, 3)
          c1 = (vkc / log((ztop(p) + z0mg(p))/z0mg(p)))**2
          beta_neutral = min( sqrt(c1 + cr*(lai(p)+sai(p))), beta_neutral_max)
       case (2)
          beta_neutral = vkc / 2._r8
    end select

    ! Calculate beta = u* / u(h) for current Obukhov length

    call GetBeta (beta_neutral, Lc(p)/obu(p), beta)

    ! Displacement height: dt is height below canopy top (dt = ztop - zdisp)

    dt = beta**2 * Lc(p)
    dt = dt * (1._r8 - exp(-0.25_r8*(lai(p)+sai(p))/beta**2))
    dt = min(ztop(p), dt)
    zdisp(p) = ztop(p) - dt

    if ((zref(p) - zdisp(p)) < 0._r8) then
       call endrun (msg=' ERROR: CanopyTurbulenceMod: ObuFunc error, zdisp height > zref')
    end if

    ! Calculate turbulent Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length. Only need Pr because Sc = Pr.

    select case (turb_type)
       case (1, 3)
          call GetPrSc (beta_neutral, beta_neutral_max, Lc(p)/obu(p), PrSc(p))
       case (2)
          PrSc(p) = 1._r8
    end select

    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psi_hat functions.
    ! These are evaluated at the reference height and at the canopy height.
    ! Here, limit Obukhov length based on values of zeta so that extreme
    ! cases are excluded.

    zeta = (zref(p) - zdisp(p)) / obu(p)
    if (zeta >= 0._r8) then
       zeta = min(1._r8, max(zeta,0.01_r8))
    else
       zeta = max(-2._r8, min(zeta,-0.01_r8))
    end if
    obu(p) = (zref(p) - zdisp(p)) / zeta

    call GetPsiRSL (zref(p), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim, psic)

    ! Friction velocity

    zlog = log((zref(p)-zdisp(p)) / (ztop(p)-zdisp(p)))
    ustar(p) = uref(p) * vkc / (zlog + psim)

    ! Wind speed at canopy height

    uaf(p) = ustar(p) / beta

    ! Temperature scale

    tstar(p) = (thref(p) - taf(p)) * vkc / (zlog + psic)

    ! Water vapor scale - use units of specific humidity (kg/kg)

    qstar(p) = (qref(p) - qaf(p)) * vkc / (zlog + psic)

    ! Aerodynamic conductance to canopy height

    gah(p) = rhomol(p) * vkc * ustar(p) / (zlog + psic)
    obu_gah(p) = obu(p)

    ! Calculate new Obukhov length (m)

    tvstar = tstar(p) + 0.61_r8 * thref(p) * qstar(p)
    obu(p) = ustar(p)**2 * thvref(p) / (vkc * grav * tvstar)

    ! Change in Obukhov length (m)

    obu_dif = obu(p) - obu_val

    end associate
  end function ObuFunc

  !-----------------------------------------------------------------------
  subroutine GetBeta (beta_neutral, LcL, beta)
    !
    ! !DESCRIPTION:
    ! Calculate beta = u* / u(h) for current Obukhov length
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral    ! Neutral value for beta = u*/u(h), ~0.3-0.5
    real(r8), intent(in)  :: LcL             ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: beta            ! Value of u*/u(h) at canopy top
    !
    ! !LOCAL VARIABLES:
    real(r8) :: aa, bb, cc, dd, qq, rr       ! Terms for quadratic or cubic solution
    real(r8) :: LcL_val                      ! LcL with limits applied
    !---------------------------------------------------------------------

    LcL_val = LcL
    
    if (LcL_val <= 0._r8) then

       ! Unstable case: quadratic equation for beta^2 at LcL_val

       bb = 16._r8 * LcL_val * beta_neutral**4     
       beta = sqrt( 0.5_r8*(-bb + sqrt(bb**2 + 4._r8 * beta_neutral**4)) )

    else

       ! Stable case: cubic equation for beta at LcL_val

       aa = 5._r8 * LcL_val
       bb = 0._r8
       cc = 1._r8
       dd = -beta_neutral
       qq = (2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)**2 - 4._r8*(bb**2 - 3._r8*aa*cc)**3
       qq = sqrt(qq)
       rr = 0.5_r8 * (qq + 2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)
       rr = rr**(1._r8/3._r8)
       beta = -(bb+rr)/(3._r8*aa) - (bb**2 - 3._r8*aa*cc)/(3._r8*aa*rr)    

    end if

    beta = min(0.50_r8, max(beta,0.20_r8))

  end subroutine GetBeta

  !-----------------------------------------------------------------------
  subroutine GetPrSc (beta_neutral, beta_neutral_max, LcL, Pr)
    !
    ! !DESCRIPTION:
    ! Calculate turbulent Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral     ! Neutral value for beta = u*/u(h), ~0.3-0.5
    real(r8), intent(in)  :: beta_neutral_max ! Maximum value for beta in neutral conditions
    real(r8), intent(in)  :: LcL              ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: Pr               ! Prandtl number at canopy top
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: Prn = 0.5_r8       ! Neutral value for Pr
    real(r8), parameter :: Prvr = 0.3_r8      ! Magnitude of variation of Pr with stability
    real(r8), parameter :: Prsc = 2.0_r8      ! Scale of variation of Pr with stability

    real(r8) :: Sc                            ! Schmidt number at canopy top
    real(r8), parameter :: Scn = Prn          ! Neutral value for Sc
    real(r8), parameter :: Scvr = Prvr        ! Magnitude of variation of Sc with stability
    real(r8), parameter :: Scsc = Prsc        ! Scale of variation of Sc with stability
    !---------------------------------------------------------------------
    
    Pr = Prn + Prvr * tanh(Prsc*LcL)
    Pr = (1._r8 - beta_neutral/beta_neutral_max) * 1._r8 + beta_neutral/beta_neutral_max*Pr

    Sc = Scn + Scvr * tanh(Scsc*LcL)
    Sc = (1._r8 - beta_neutral/beta_neutral_max) * 1._r8 + beta_neutral/beta_neutral_max*Sc

  end subroutine GetPrSc

  !-----------------------------------------------------------------------
  subroutine GetPsiRSL (za, hc, dt, beta, zeta, obu, PrSc, psim, psic)
    !
    ! !DESCRIPTION:
    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psihat functions.
    ! These are evaluated between the height za and at the canopy height hc.
    !
    ! !USES:
    use clm_varctl, only : turb_type
    use clm_varcon, only : vkc
    use clm_varcon, only : dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: za        ! Atmospheric height (m)
    real(r8), intent(in)  :: hc        ! Canopy height (m)
    real(r8), intent(in)  :: dt        ! Height below canopy top (dt = ztop - zdisp)
    real(r8), intent(in)  :: beta      ! Value of u*/u(h) at canopy top
    real(r8), intent(in)  :: zeta      ! Monin-Obukhov stability parameter (za-d)/L
    real(r8), intent(in)  :: obu       ! Obukhov length (m)
    real(r8), intent(in)  :: PrSc      ! Turbulent Prandtl (Schmidt) number at canopy top
    real(r8), intent(out) :: psim      ! psi function for momentum including RSL influence
    real(r8), intent(out) :: psic      ! psi function for scalars  including RSL influence
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phim                   ! Monin-Obukhov phi function for momentum at canopy top
    real(r8) :: phic                   ! Monin-Obukhov phi function for scalars at canopy top
    real(r8) :: c1                     ! RSL magnitude multiplier
    real(r8) :: psihat1                ! Multiply by c1 to get the RSL psihat function evaluated at za
    real(r8) :: psihat2                ! Multiply by c1 to get the RSL psihat function evaluated at hc
    !---------------------------------------------------------------------

    ! In the RSL theory, c1 and c2 represent the scaled magnitude and
    ! height over which the RSL theory modifies MOST via the psihat functions:
    !
    !
    !     z
    !     ^
    !     |                                    . ___
    !     |                                    .  ^
    !     |                                   .   |
    !     |                                  .    |
    !     |                                 .     |
    !     |                               .       |
    !     |                             .        c2
    !     |                          .            |
    !     |                       .               |
    !     |                  .                    |
    !     |           .                           |
    !     |   .                                 _\/_
    !     -------------------------------------------> u
    !         |<-------------- c1 --------------->|

    ! Evaluate the roughness sublayer psihat function for momentum at
    ! the height za and at the canopy height hc. Values for psihat are obtained
    ! from a look-up table. Here, heights are above the canopy for compatibility
    ! with the supplied look-up table. These heights are also scaled by dt = hc-d
    ! so that the look-up table uses (za-hc)/dt and (hc-hc)/dt. Also the term
    ! dt in the integration of psihat is scaled by the Obukhov length L (dt/obu).
    ! This means that the returned psihat value needs to be scaled (multiplied) by
    ! c1 before it fully represents psihat as it appears in the RSL equations.

    select case (turb_type)
       case (1, 3)
          call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat1)
          call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psihat2)
       case (2)
          psihat1 = 0._r8
          psihat2 = 0._r8
    end select

    phim = phim_monin_obukhov(dt/obu)
    c1 = (1._r8 - vkc / (2._r8 * beta * phim)) * exp(0.5_r8*c2)
    psim = -psim_monin_obukhov(zeta) + psim_monin_obukhov(dt/obu) + c1*(psihat1 - psihat2) + vkc / beta

    ! Now do the same for heat

    select case (turb_type)
       case (1, 3)
          call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat1)
          call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psihat2)
       case (2)
          psihat1 = 0._r8
          psihat2 = 0._r8
    end select

    phic = phic_monin_obukhov(dt/obu)
    c1 = (1._r8 - PrSc*vkc / (2._r8 * beta * phic)) * exp(0.5_r8*c2)
    psic = -psic_monin_obukhov(zeta) + psic_monin_obukhov(dt/obu) + c1*(psihat1 - psihat2)

  end subroutine GetPsiRSL

  !-----------------------------------------------------------------------
  subroutine LookupPsihat (zdt, dtL, zdtgrid, dtLgrid, psigrid, psihat)
    !
    ! !DESCRIPTION:
    ! Determines the RSL function psihat as provided through a look-up table
    ! for input values of zdt and dtL. Linearly interpolates between values
    ! supplied on the look-up table grid defined by zdtgrid, dtLgrid, psigrid.
    !
    ! NOTE: The psihat presented in Harman and Finnigan (2007,2008) and Harman (2012)
    ! has been re-written in non-dimensional form such that it now appears as:
    !
    ! psihat(z) = c1 * A(z/(beta^2*Lc),(beta^2*Lc)/L)
    !
    ! This routine gets the value of A from a look-up table. Noting that dt=beta^2*Lc,
    ! this routine therefore requires values of z/dt and dt/L. In addition, this means
    ! that the returned psihat value needs to be scaled (multiplied) by c1 before it fully
    ! represents psihat as it appears in the RSL equations.
    !
    ! !USES:
    use clm_varcon, only : nZ, nL
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zdt            ! Height (above canopy) normalized by dt
    real(r8), intent(in) :: dtL            ! dt/L (displacement height/Obukhov length)
    real(r8), intent(in) :: zdtgrid(nZ,1)  ! Grid of zdt on which psihat is given
    real(r8), intent(in) :: dtLgrid(1,nL)  ! Grid of dtL on which psihat is given
    real(r8), intent(in) :: psigrid(nZ,nL) ! Grid of psihat values
    real(r8), intent(out):: psihat         ! Value of psihat
    !
    !LOCAL VARIABLES
    integer  :: ii, jj                     ! Looping indices
    integer  :: L1, L2, Z1, Z2             ! Grid indices for psihat sought
    real(r8) :: wL1, wL2, wZ1, wZ2         ! Weights for averaging
    !---------------------------------------------------------------------

    ! Find indices and weights for dtL values which bracket the specified dtL

    L1 = 0; L2 = 0
    if ( dtL <= dtLgrid(1,1) ) then
       L1 = 1
       L2 = 1
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else if ( dtL > dtLgrid(1,nL) ) then
       L1 = nL
       L2 = nL
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else
       do jj = 1, nL-1
          if ( (dtL <= dtLgrid(1,jj+1)) .and. (dtL > dtLgrid(1,jj)) ) then
             L1 = jj
             L2 = jj + 1
             wL1 = (dtLgrid(1,L2) - dtL) / (dtLgrid(1,L2) - dtLgrid(1,L1))
             wL2 = 1._r8 - wL1
          end if
       end do
    end if

    if (L1 == 0 .or. L2 == 0) then
       call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihat error, indices L1 and L2 not found')
    end if

    ! Find indices and weights for zdt values which bracket the specified zdt

    Z1 = 0; Z2 = 0
    if ( zdt > zdtgrid(1,1)) then
       Z1 = 1
       Z2 = 1
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else if (zdt < zdtgrid(nZ,1) ) then
       Z1 = nZ
       Z2 = nZ
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else
       do ii = 1, nZ-1
          if ( (zdt >= zdtgrid(ii+1,1)) .and. (zdt < zdtgrid(ii,1)) ) then
             Z1 = ii
             Z2 = ii + 1
             wZ1 = (zdt - zdtgrid(ii+1,1)) / (zdtgrid(ii,1) - zdtgrid(ii+1,1))
             wZ2 = 1._r8 - wZ1
          end if
       end do
    end if

    if (Z1 == 0 .or. Z2 == 0) then
       call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihat error, indices Z1 and Z2 not found')
    end if

    ! Calculate psihat as a weighted average of the values of psihat on the grid

    psihat = wZ1*wL1*psigrid(Z1,L1) + wZ2*wL1*psigrid(Z2,L1) &
           + wZ1*wL2*psigrid(Z1,L2) + wZ2*wL2*psigrid(Z2,L2)

  end subroutine LookupPsihat

  !-----------------------------------------------------------------------
  subroutine ScalarProfileIterative (ncan, cref, ga, s, fac, cval0, cval)
    !
    ! !DESCRIPTION:
    ! Form a vertical aerodynamic conductance network (ga) that when connected
    ! with the network of leaf-level sources (s) forms a tridiagonal system of
    ! equations that is solved to give the within-canopy vertical scalar profile (cval).
    ! Boundary conditions are the above-canopy scalar value at the reference height (cref)
    ! and the ground flux (s(0)).
    !
    ! !USES:
    use clm_varpar, only : nlevcan
    use MathToolsMod, only : tridiag
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: ncan             ! Number of layers
    real(r8), intent(in)  :: cref             ! Scalar value at reference height
    real(r8), intent(in)  :: ga(0:nlevcan)    ! Aerodynamic conductance across each canopy layer
    real(r8), intent(in)  :: s(0:nlevcan)     ! Scalar source for each canopy layer
    real(r8), intent(in)  :: fac              ! Scaling factor for fluxes
    real(r8), intent(in)  :: cval0(0:nlevcan) ! Scalar value at each canopy layer for previous timestep
    real(r8), intent(out) :: cval(0:nlevcan)  ! Scalar value at each canopy layer for current timestep
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                            ! Canopy layer index
    real(r8) :: a(ncan),b(ncan)               ! Entries in tridiagonal matrix
    real(r8) :: c(ncan),d(ncan)               ! Entries in tridiagonal matrix
    real(r8) :: u(ncan)                       ! Tridiagonal solution
    integer  :: method                        ! Solution type: (1) uses cval (2) uses 0.5(cval0+cval)
    !---------------------------------------------------------------------

    method = 2

    select case (method)

       ! Method 1 uses scalar values at time n+1 (cval)

       case (1)

       ! Build tridiagonal vectors at ic = 1

       ic = 1
       a(ic) = 0._r8
       b(ic) = ga(ic)
       c(ic) = -ga(ic)
       d(ic) = (s(0) + s(ic)) * fac

       ! Build tridiagonal vectors for ic = 2 -> ncan-1

       do ic = 2, ncan-1
          a(ic) = -ga(ic-1)
          b(ic) = ga(ic-1) + ga(ic)
          c(ic) = -ga(ic)
          d(ic) = s(ic) * fac
       end do

       ! Build tridiagonal vectors at ic = ncan

       ic = ncan
       a(ic) = -ga(ic-1)
       b(ic) = ga(ic-1) + ga(ic)
       c(ic) = 0._r8
       d(ic) = cref * ga(ic) + s(ic) * fac

       ! Method 2 uses scalar values at time n (cval0) and n+1 (cval)
       ! ie., 0.5 cval(n) + 0.5 cval(n+1)

       case (2)

       ! Build tridiagonal vectors at ic = 1

       ic = 1
       a(ic) = 0._r8
       b(ic) = 0.5_r8 * ga(ic)
       c(ic) = -0.5_r8 * ga(ic)
       d(ic) = (s(0) + s(ic)) * fac - 0.5_r8 * ga(ic) * cval0(ic) + 0.5_r8 * ga(ic) * cval0(ic+1)

       ! Build tridiagonal vectors for ic = 2 -> ncan-1

       do ic = 2, ncan-1
          a(ic) = -0.5_r8 * ga(ic-1)
          b(ic) = 0.5_r8 * (ga(ic-1) + ga(ic))
          c(ic) = -0.5_r8 * ga(ic)
          d(ic) = s(ic) * fac + 0.5_r8 * ga(ic-1) * cval0(ic-1) &
                - 0.5_r8 * (ga(ic-1) + ga(ic)) * cval0(ic) + 0.5_r8 * ga(ic) * cval0(ic+1)
       end do

       ! Build tridiagonal vectors at ic = ncan

       ic = ncan
       a(ic) = -0.5_r8 * ga(ic-1)
       b(ic) = 0.5_r8 * (ga(ic-1) + ga(ic))
       c(ic) = 0._r8
       d(ic) = cref * ga(ic) + s(ic) * fac &
             + 0.5_r8 * ga(ic-1) * cval0(ic-1) - 0.5_r8 * (ga(ic-1) + ga(ic)) * cval0(ic)

    end select

    ! Solve tridiagonal system for concentration profile

    call tridiag (a, b, c, d, u, ncan)

    ! Copy concentration profile (needed because u is for levels 1 -> ncan
    ! but cval is for levels 0 -> ncan)

    do ic = 1, ncan
       cval(ic) = u(ic)
    end do

  end subroutine ScalarProfileIterative

  !-----------------------------------------------------------------------
  subroutine ScalarProfile (p, niter, soilstate_inst, temperature_inst, &
  energyflux_inst, waterflux_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Form a vertical aerodynamic conductance network that when connected with
    ! the network of leaf-level sources forms a system of equations that is solved
    ! to give the within-canopy vertical scalar profiles of temperature and water vapor,
    ! leaf temperature, and the leaf source fluxes. The boundary condition is the
    ! above-canopy scalar values at the reference height, and the ground
    ! temperature and fluxes are calculated as part of the implicit solution.
    !
    ! !USES:
    use clm_varpar, only : isun, isha, nlevcan, nleaf
    use clm_varctl, only : dtime_sub, leaf_temp_iter
    use clm_time_manager, only : get_nstep
    use WaterVaporMod, only : SatVap, LatVap
    use ColumnType, only : col
    use PatchType, only : patch
    use SoilStateType, only : soilstate_type
    use TemperatureType, only : temperature_type
    use EnergyFluxType, only : energyflux_type
    use WaterFluxType, only : waterflux_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    use LeafTemperatureMod, only : LeafTemperature
    use SoilFluxesMultilayerMod, only : SoilFluxesMultilayer
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p         ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: niter     ! Iteration index
    type(soilstate_type), intent(inout) :: soilstate_inst
    type(temperature_type), intent(inout) :: temperature_inst
    type(energyflux_type), intent(inout) :: energyflux_inst
    type(waterflux_type), intent(inout) :: waterflux_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c                    ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                   ! Canopy layer index
    integer  :: il                   ! Sunlit (1) or shaded (2) leaf index
    integer  :: nstep                ! Current model time step number

    real(r8) :: dtime                ! Model time step (s)
    real(r8) :: lambda               ! Latent heat of vaporization (J/mol)
    real(r8) :: esat                 ! Saturation vapor pressure (Pa)
    real(r8) :: desat                ! Temperature derivative of saturation vapor pressure (Pa/K)

    real(r8) :: qsat0                ! Saturation vapor pressure at ground (mol/mol)
    real(r8) :: dqsat0               ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: gsw                  ! Soil conductance for water vapor (mol H2O/m2/s)
    real(r8) :: gs0                  ! Total soil-to-air conductance for water vapor (mol H2O/m2/s)
    real(r8) :: alpha0               ! Coefficient for ground temperature (dimensionless)
    real(r8) :: beta0                ! Coefficient for ground temperature (K mol/mol)
    real(r8) :: delta0               ! Coefficient for ground temperature (K)
    real(r8) :: c02                  ! Soil heat flux term (W/m2/K)
    real(r8) :: c01                  ! Soil heat flux term (W/m2)
    real(r8) :: den                  ! Intermediate calculation
    real(r8) :: t0                   ! Soil surface temperature (K)
    real(r8) :: sh0                  ! Ground sensible heat flux (W/m2)
    real(r8) :: et0                  ! Ground evaporation flux (mol/m2/s)
    real(r8) :: g0                   ! Soil heat flux (W/m2)

    real(r8) :: pai                          ! Plant area of layer (m2/m2)
    real(r8) :: gleaf_sh(nlevcan,nleaf)      ! Leaf conductance for sensible heat (mol/m2/s)
    real(r8) :: gleaf_et(nlevcan,nleaf)      ! Leaf conductance for water vapor (mol/m2/s)
    real(r8) :: heatcap(nlevcan,nleaf)       ! Heat capacity of leaves (J/m2/K)
    real(r8) :: avail_energy(nlevcan,nleaf)  ! Available energy for leaf (W/m2)
    real(r8) :: qsat                         ! Saturation vapor pressure (mol/mol)
    real(r8) :: dqsat(nlevcan,nleaf)         ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: qsat_term(nlevcan,nleaf)     ! Intermediate calculation for saturation vapor pressure (mol/mol)
    real(r8) :: alpha(nlevcan,nleaf)         ! Coefficient for leaf temperature (dimensionless)
    real(r8) :: beta(nlevcan,nleaf)          ! Coefficient for leaf temperature (K mol/mol)
    real(r8) :: delta(nlevcan,nleaf)         ! Coefficient for leaf temperature (K)

    real(r8) :: rho_dz_over_dt               ! Intermediate calculation for canopy air storage

    real(r8) :: a1(nlevcan)                  ! Coefficient for canopy air temperature
    real(r8) :: b11(nlevcan)                 ! Coefficient for canopy air temperature
    real(r8) :: b12(nlevcan)                 ! Coefficient for canopy air temperature
    real(r8) :: c1(nlevcan)                  ! Coefficient for canopy air temperature
    real(r8) :: d1(nlevcan)                  ! Coefficient for canopy air temperature

    real(r8) :: a2(nlevcan)                  ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: b21(nlevcan)                 ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: b22(nlevcan)                 ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: c2(nlevcan)                  ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: d2(nlevcan)                  ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: ga_prof_ic_minus_one         ! Special case for ic-1 = 0: use gs0 not ga_prof

    real(r8) :: ainv, binv                   ! "a" and "b" elements of 2x2 matrix to invert
    real(r8) :: cinv, dinv                   ! "c" and "d" elements of 2x2 matrix to invert
    real(r8) :: det                          ! Determinant of 2x2 matrix

    real(r8) :: e11(0:nlevcan)               ! Coefficient for canopy air temperature
    real(r8) :: e12(0:nlevcan)               ! Coefficient for canopy air temperature
    real(r8) :: f1(0:nlevcan)                ! Coefficient for canopy air temperature
    real(r8) :: e21(0:nlevcan)               ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: e22(0:nlevcan)               ! Coefficient for canopy air water vapor (mole fraction)
    real(r8) :: f2(0:nlevcan)                ! Coefficient for canopy air water vapor (mole fraction)

    ! These are needed only for conservation checks

    real(r8) :: storage_sh(nlevcan)        ! Heat storage flux in air (W/m2)
    real(r8) :: storage_et(nlevcan)        ! Water vapor storage flux in air (mol/m2/s)
    real(r8) :: stveg(nlevcan)             ! Leaf storage heat flux (W/m2)
    real(r8) :: shsrc(nlevcan)             ! Leaf sensible heat flux (W/m2)
    real(r8) :: etsrc(nlevcan)             ! Leaf water vapor flux (mol/m2/s)
    real(r8) :: stveg_leaf                 ! Leaf storage heat flux from LeafTemperatureMod (W/m2)
    real(r8) :: shsrc_leaf                 ! Leaf sensible heat flux from LeafTemperatureMod (W/m2)
    real(r8) :: etsrc_leaf                 ! Leaf water vapor flux from LeafTemperatureMod (mol/m2/s)
    real(r8) :: error                      ! Energy imbalance (W/m2)
    real(r8) :: sum_src                    ! Sum of source flux over leaf layers
    real(r8) :: sum_storage                ! Sum of storage flux over leaf layers

    integer, parameter :: tridiag = 1      ! 1 = use tridiagonal solver. 0 = use lapack solver
    integer :: neq                         ! LAPACK: number of linear equations to solve
    integer :: nrhs                        ! LAPACK: number of right-hand sides
    integer :: info                        ! LAPACK: error flag
    integer :: meq                         ! Index to linear equation matrix

    real(r8), allocatable :: asolv(:,:)    ! LAPACK: A coefficient matrix
    real(r8), allocatable :: bsolv(:)      ! LAPACK: B matrix
    integer,  allocatable :: ipiv(:)       ! LAPACK: vector of pivot indices

    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    ntop      => mlcanopy_inst%ntop            , &  ! Index for top leaf layer
    ncan      => mlcanopy_inst%ncan            , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai            , &  ! Layer plant area index (m2/m2)
    fwet      => mlcanopy_inst%fwet            , &  ! Fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry            , &  ! Fraction of plant area index that is green and dry
    zw        => mlcanopy_inst%zw              , &  ! Canopy heights at layer interfaces (m)
    tref      => mlcanopy_inst%tref            , &  ! Air temperature at reference height (K)
    thref     => mlcanopy_inst%thref           , &  ! Atmospheric potential temperature (K)
    eref      => mlcanopy_inst%eref            , &  ! Vapor pressure at reference height (Pa)
    pref      => mlcanopy_inst%pref            , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol          , &  ! Molar density at reference height (mol/m3)
    cpair     => mlcanopy_inst%cpair           , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    gbh       => mlcanopy_inst%gbh             , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv             , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs              , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    cpleaf    => mlcanopy_inst%cpleaf          , &  ! Leaf heat capacity (J/m2 leaf/K)
    rnleaf    => mlcanopy_inst%rnleaf          , &  ! Leaf net radiation (W/m2 leaf)
    fracsun   => mlcanopy_inst%fracsun         , &  ! Sunlit fraction of canopy layer
    fracsha   => mlcanopy_inst%fracsha         , &  ! Shaded fraction of canopy layer
    ga_prof   => mlcanopy_inst%ga_prof         , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    tair_old  => mlcanopy_inst%tair_old        , &  ! Air temperature profile for previous timestep (K)
    eair_old  => mlcanopy_inst%eair_old        , &  ! Vapor pressure profile for previous timestep (Pa)
    tveg_old  => mlcanopy_inst%tveg_old        , &  ! Vegetation temperature profile for previous timestep (K)
    rhg       => mlcanopy_inst%rhg             , &  ! Relative humidity of airspace at soil surface (fraction)
    rnsoi     => mlcanopy_inst%rnsoi           , &  ! Net radiation, ground (W/m2)
    z         => col%z                         , &  ! Soil layer depth (m)
    zi        => col%zi                        , &  ! Soil layer depth at layer interface (m)
    snl       => col%snl                       , &  ! Number of snow layers
    thk       => soilstate_inst%thk_col        , &  ! Soil layer thermal conductivity (W/m/K)
    soilresis => soilstate_inst%soilresis_col  , &  ! Soil evaporative resistance (s/m)
    t_soisno  => temperature_inst%t_soisno_col , &  ! Soil temperature (K)
                                                    ! *** Output ***
    tair      => mlcanopy_inst%tair            , &  ! Air temperature profile (K)
    eair      => mlcanopy_inst%eair            , &  ! Vapor pressure profile (Pa)
    tveg      => mlcanopy_inst%tveg            , &  ! Vegetation temperature profile (K)
    shair     => mlcanopy_inst%shair           , &  ! Canopy air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair           , &  ! Canopy air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair           , &  ! Canopy air storage heat flux (W/m2)
                                                    ! *** These are only because leaf and soil fluxes are done here again
    tleaf     => mlcanopy_inst%tleaf           , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf          , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf          , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf          , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf          , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf          , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    shsoi     => mlcanopy_inst%shsoi           , &  ! Sensible heat flux, ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi           , &  ! Latent heat flux, ground (W/m2)
    gsoi      => mlcanopy_inst%gsoi            , &  ! Soil heat flux (W/m2)
    etsoi     => mlcanopy_inst%etsoi           , &  ! Water vapor flux, ground (mol H2O/m2/s)
    tg        => mlcanopy_inst%tg                &  ! Soil surface temperature (K)
    )

    c = patch%column(p)

    ! Get current step (counter)

    nstep = get_nstep()

    ! Time step (s)

    dtime = dtime_sub

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    !---------------------------------------------------------------------
    ! Terms for ground temperature, which is calculated from the energy balance:
    ! Rn0 - H0 - lambda*E0 - G = 0, and is rewritten as:
    ! T0 = alpha0*T(1) + beta0*q(1) + delta0
    !---------------------------------------------------------------------

    call SatVap (tair_old(p,0), esat, desat)           ! Vapor pressure (Pa) at ground temperature
    qsat0 = esat / pref(p) ; dqsat0 = desat / pref(p)  ! Pa -> mol/mol

    gsw = 1._r8 / soilresis(c)                         ! Soil conductance for water vapor: s/m -> m/s
    gsw = gsw * rhomol(p)                              ! m/s -> mol H2O/m2/s
    gs0 = ga_prof(p,0) * gsw / (ga_prof(p,0) + gsw)    ! Total conductance

    c02 = thk(c,snl(c)+1) / (z(c,snl(c)+1)-zi(c,snl(c)))  ! Soil heat flux term (W/m2/K)
    c01 = -c02 * t_soisno(c,snl(c)+1)                     ! Soil heat flux term (W/m2)

    den = cpair(p) * ga_prof(p,0) + lambda * rhg(p) * gs0 * dqsat0 + c02
    alpha0 = cpair(p) * ga_prof(p,0) / den
    beta0 = lambda * gs0 / den
    delta0 = (rnsoi(p) - lambda * rhg(p) * gs0 * (qsat0 - dqsat0 * tair_old(p,0)) - c01) / den

    !---------------------------------------------------------------------
    ! alpha, beta, delta coefficients for leaf temperature:
    !
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
    !
    ! and leaf terms needed for scalar conservations equations. These
    ! are defined for all ncan layers but the leaf terms only exist
    ! for the layers with leaves.
    !---------------------------------------------------------------------

    do ic = 1, ncan(p)

       if (dpai(p,ic) > 0.0_r8) then ! leaf layer

          ! Calculate terms for sunlit and shaded leaves

          do il = 1, nleaf

             ! Leaf conductances (here these are per unit leaf area)

             gleaf_sh(ic,il) = 2._r8 * gbh(p,ic,il)
             gleaf_et(ic,il) = gs(p,ic,il)*gbv(p,ic,il)/(gs(p,ic,il)+gbv(p,ic,il)) * fdry(p,ic) + gbv(p,ic,il) * fwet(p,ic)

             ! Heat capacity of leaves

             heatcap(ic,il) = cpleaf(p,ic)

             ! Available energy: net radiation

             avail_energy(ic,il) = rnleaf(p,ic,il)

             ! Saturation vapor pressure for leaf temperature at time n: Pa -> mol/mol

             call SatVap (tveg_old(p,ic,il), esat, desat)
             qsat = esat / pref(p) ; dqsat(ic,il) = desat / pref(p)

             ! Term for linearized vapor pressure at leaf temperature:
             ! qsat(tveg) = qsat(tveg_old) + dqsat * (tveg - tveg_old)
             ! Here qsat_term contains the terms with tveg_old

             qsat_term(ic,il) = qsat - dqsat(ic,il) * tveg_old(p,ic,il)

             ! alpha, beta, delta coefficients for leaf temperature

             den = heatcap(ic,il) / dtime + gleaf_sh(ic,il) * cpair(p) + gleaf_et(ic,il) * lambda * dqsat(ic,il)
             alpha(ic,il) = gleaf_sh(ic,il) * cpair(p) / den
             beta(ic,il) = gleaf_et(ic,il) * lambda / den
             delta(ic,il) = avail_energy(ic,il) / den &
                          - lambda * gleaf_et(ic,il) * qsat_term(ic,il) / den &
                          + heatcap(ic,il) / dtime * tveg_old(p,ic,il) / den

             ! Now scale flux terms for leaf area

             if (il == isun) pai = fracsun(p,ic) * dpai(p,ic)
             if (il == isha) pai = fracsha(p,ic) * dpai(p,ic)

             gleaf_sh(ic,il) = gleaf_sh(ic,il) * pai
             gleaf_et(ic,il) = gleaf_et(ic,il) * pai
             heatcap(ic,il) = heatcap(ic,il) * pai
             avail_energy(ic,il) = avail_energy(ic,il) * pai

          end do

       else ! non-leaf layer

          ! Zero out terms

          do il = 1, nleaf
             gleaf_sh(ic,il) = 0._r8
             gleaf_et(ic,il) = 0._r8
             heatcap(ic,il) = 0._r8
             avail_energy(ic,il) = 0._r8
             dqsat(ic,il) = 0._r8
             qsat_term(ic,il) = 0._r8
             alpha(ic,il) = 0._r8
             beta(ic,il) = 0._r8
             delta(ic,il) = 0._r8
          end do

       end if

    end do

    !---------------------------------------------------------------------
    ! a,b,c,d coefficients for air temperature:
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    !
    ! a,b,c,d coefficients for water vapor (mole fraction):
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !---------------------------------------------------------------------

    do ic = 1, ncan(p)

       ! Storage term

       rho_dz_over_dt = rhomol(p) * (zw(p,ic) - zw(p,ic-1)) / dtime

       ! a,b,c,d coefficients for air temperature

       a1(ic) = -ga_prof(p,ic-1)
       b11(ic) = rho_dz_over_dt + ga_prof(p,ic-1) + ga_prof(p,ic) &
               + gleaf_sh(ic,isun) * (1._r8 - alpha(ic,isun)) + gleaf_sh(ic,isha) * (1._r8 - alpha(ic,isha))
       b12(ic) = -gleaf_sh(ic,isun) * beta(ic,isun) - gleaf_sh(ic,isha) * beta(ic,isha)
       c1(ic) = -ga_prof(p,ic)
       d1(ic) = rho_dz_over_dt * tair_old(p,ic) + gleaf_sh(ic,isun) * delta(ic,isun) + gleaf_sh(ic,isha) * delta(ic,isha)

       ! Special case for top canopy layer

       if (ic == ncan(p)) then
          c1(ic) = 0._r8
          d1(ic) = d1(ic) + ga_prof(p,ic) * thref(p)
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a1(ic) = 0._r8
          b11(ic) = b11(ic) - ga_prof(p,ic-1) * alpha0
          b12(ic) = b12(ic) - ga_prof(p,ic-1) * beta0
          d1(ic) = d1(ic) + ga_prof(p,ic-1) * delta0
       end if

       ! a,b,c,d coefficients for water vapor (mole fraction)

       if (ic == 1) then
          ga_prof_ic_minus_one = gs0
       else
          ga_prof_ic_minus_one = ga_prof(p,ic-1)
       end if

       a2(ic) = -ga_prof_ic_minus_one
       b21(ic) = -gleaf_et(ic,isun) * dqsat(ic,isun) * alpha(ic,isun) - gleaf_et(ic,isha) * dqsat(ic,isha) * alpha(ic,isha)
       b22(ic) = rho_dz_over_dt + ga_prof_ic_minus_one + ga_prof(p,ic) &
               + gleaf_et(ic,isun) * (1._r8 - dqsat(ic,isun) * beta(ic,isun)) &
               + gleaf_et(ic,isha) * (1._r8 - dqsat(ic,isha) * beta(ic,isha))
       c2(ic) = -ga_prof(p,ic)
       d2(ic) = rho_dz_over_dt * eair_old(p,ic) / pref(p) &
              + gleaf_et(ic,isun) * (dqsat(ic,isun) * delta(ic,isun) + qsat_term(ic,isun)) &
              + gleaf_et(ic,isha) * (dqsat(ic,isha) * delta(ic,isha) + qsat_term(ic,isha)) 

       ! Special case for top canopy layer

       if (ic == ncan(p)) then
          c2(ic) = 0._r8
          d2(ic) = d2(ic) + ga_prof(p,ic) * eref(p) / pref(p)
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a2(ic) = 0._r8
          b21(ic) = b21(ic) - gs0 * rhg(p) * dqsat0 * alpha0
          b22(ic) = b22(ic) - gs0 * rhg(p) * dqsat0 * beta0
          d2(ic) = d2(ic) + gs0 * rhg(p) * (qsat0 + dqsat0 * (delta0 - tair_old(p,0)))
       end if

    end do

    !---------------------------------------------------------------------
    ! Use tridiagonal solver
    !
    ! Solve for air temperature and water vapor (mole fraction):
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! The solution rewrites these equations so that:
    ! T(i) = f1(i) - e11(i)*T(i+1) - e12(i)*q(i+1) 
    ! q(i) = f2(i) - e21(i)*T(i+1) - e22(i)*q(i+1) 
    !
    ! Note that as used here eair = mol/mol
    !---------------------------------------------------------------------

    if (tridiag == 1) then

    e11(0) = 0._r8
    e12(0) = 0._r8
    e21(0) = 0._r8
    e22(0) = 0._r8
    f1(0) = 0._r8
    f2(0) = 0._r8

    do ic = 1, ncan(p)

       ! Elements in the 2x2 matrix to invert: 
       !
       !                      | a b |
       ! B(i) - A(i)*E(i-1) = | c d |
       !
       ! and determinant of the matrix

       ainv = b11(ic) - a1(ic) * e11(ic-1)
       binv = b12(ic) - a1(ic) * e12(ic-1)
       cinv = b21(ic) - a2(ic) * e21(ic-1)
       dinv = b22(ic) - a2(ic) * e22(ic-1)
       det = ainv * dinv - binv * cinv

       ! E(i) = [B(i) - A(i)*E(i-1)]^(-1) * C(i)
       e11(ic) = dinv * c1(ic) / det
       e12(ic) = -binv * c2(ic) / det
       e21(ic) = -cinv * c1(ic) / det
       e22(ic) = ainv * c2(ic) / det

       ! F(i) = [B(i) - A(i)*E(i-1)]^(-1) * [D(i) - A(i)*F(i-1)]
       f1(ic) =  (dinv*(d1(ic) - a1(ic)*f1(ic-1)) - binv*(d2(ic) - a2(ic)*f2(ic-1))) / det
       f2(ic) = (-cinv*(d1(ic) - a1(ic)*f1(ic-1)) + ainv*(d2(ic) - a2(ic)*f2(ic-1))) / det

    end do

    ic = ncan(p)
    tair(p,ic) = f1(ic)
    eair(p,ic) = f2(ic)

    do ic = ncan(p)-1, 1, -1
       tair(p,ic) = f1(ic) - e11(ic)*tair(p,ic+1) - e12(ic)*eair(p,ic+1) 
       eair(p,ic) = f2(ic) - e21(ic)*tair(p,ic+1) - e22(ic)*eair(p,ic+1) 
    end do

    end if

    !---------------------------------------------------------------------
    ! Use LAPACK to solve the linear system of 2*ncan equations with 2*ncan
    ! unknowns for T and q at each level in the canopy.
    !
    ! Solve for air temperature and water vapor (mole fraction) where for each
    ! level i there are two equations for T(i) and q(i):
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! Rewrite these equations in the matrix form A*X = B where for each level i:
    !
    !       a(1)*x(1) + a(2)*x(2) + ... + a(2n-1)*x(2n-1) + a(2n)*x(2n) = b(i)
    ! with:
    !       x(1   ) = T(1)
    !       x(2   ) = q(1)
    !              ...
    !       x(2n-1) = T(n)
    !       x(2n  ) = q(n)
    !
    ! for canopy levels i = 1 (bottom layer) to i = n (top layer)
    !
    ! Solve A*X = B using LAPACK dgesv routine where A is an n-by-n matrix
    ! and X and B are n-by-m matrices.
    !
    ! call dgesv (neq, nrhs, a, lda, ipiv, b, ldb, info)
    !
    ! NEQ  (input) number of linear equations
    ! NRHS (input) number of right-hand sides (number of columns of B)
    ! A    (input/output) On entry, the LDA-by-NEQ coefficient matrix A. On exit, the results from factorization
    ! LDA  (input) leading dimension of A (equal to NEQ)
    ! IPIV (output) An integer vector of the NEQ pivot indices
    ! B    (input/output) On entry, the LDB-by-NRHS matrix B. On exit, the LDB-by-NRHS solution X if INFO = 0
    ! LDB  (input) leading dimension of B (equal to NEQ)
    ! INFO (output) equals 0 if successful. info > 1 if unsucessful and denotes the first j where U(jj)=0
    !---------------------------------------------------------------------

    if (tridiag == 0) then

    neq = ncan(p) * 2
    nrhs = 1

    allocate (asolv(neq,neq))
    allocate (bsolv(neq))
    allocate (ipiv(neq))

    asolv(:,:) = 0._r8
    bsolv(:) = 0._r8
    ipiv(:) = 0

    ! Fill in B matrix

    meq = 0
    do ic = 1, ncan(p)
       meq = meq + 1
       bsolv(meq) = d1(ic) ! first equation is T(i)
       meq = meq + 1
       bsolv(meq) = d2(ic) ! second equation is q(i)
    end do

    ! Fill in  A matrix

    meq = 0
    do ic = 1, ncan(p)

       ! First equation is for T

       meq = meq + 1
       if (ic == 1) then                ! First canopy layer (i = 1)
          asolv(meq,meq+0) = b11(ic)    ! T(i)
          asolv(meq,meq+1) = b12(ic)    ! q(i)
          asolv(meq,meq+2) = c1(ic)     ! T(i+1)
       else if (ic < ncan(p)) then      ! Middle canopy layer (2 <= i < n)
          asolv(meq,meq-2) = a1(ic)     ! T(i-1)
          asolv(meq,meq+0) = b11(ic)    ! T(i)
          asolv(meq,meq+1) = b12(ic)    ! q(i)
          asolv(meq,meq+2) = c1(ic)     ! T(i+1)
       else                             ! Top canopy layer (i = n)
          asolv(meq,meq-2) = a1(ic)     ! T(i-1)
          asolv(meq,meq+0) = b11(ic)    ! T(i)
          asolv(meq,meq+1) = b12(ic)    ! q(i)
       end if

       ! Second equation is for q

       meq = meq + 1
       if (ic == 1) then                ! First canopy layer (i = 1)
          asolv(meq,meq-1) = b21(ic)    ! T(i)
          asolv(meq,meq+0) = b22(ic)    ! q(i)
          asolv(meq,meq+2) = c2(ic)     ! q(i+1)
       else if (ic < ncan(p)) then      ! Middle canopy layer (2 <= i < n)
          asolv(meq,meq-2) = a2(ic)     ! q(i-1)
          asolv(meq,meq-1) = b21(ic)    ! T(i)
          asolv(meq,meq+0) = b22(ic)    ! q(i)
          asolv(meq,meq+2) = c2(ic)     ! q(i+1)
       else                             ! Top canopy layer (i = n)
          asolv(meq,meq-2) = a2(ic)     ! q(i-1)
          asolv(meq,meq-1) = b21(ic)    ! T(i)
          asolv(meq,meq+0) = b22(ic)    ! q(i)
       end if

    end do

    ! Call lapack solver

    call dgesv (neq, nrhs, asolv, neq, ipiv, bsolv, neq, info)
 
    if (info /= 0) then
       call endrun (msg=' ERROR: ScalarProfile: LAPACK solver error')
    end if

    ! Map solution matrix B to T and q

    meq = 0
    do ic = 1, ncan(p)
       meq = meq + 1
       tair(p,ic) = bsolv(meq)  ! T(i)
       meq = meq + 1
       eair(p,ic) = bsolv(meq)  ! q(i)
    end do

    deallocate (asolv)
    deallocate (bsolv)
    deallocate (ipiv)

    end if

    !---------------------------------------------------------------------
    ! Calculate leaf temperature:
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
    !---------------------------------------------------------------------

    do ic = 1, ncan(p)
       tveg(p,ic,isun) = alpha(ic,isun)*tair(p,ic) + beta(ic,isun)*eair(p,ic) + delta(ic,isun)
       tveg(p,ic,isha) = alpha(ic,isha)*tair(p,ic) + beta(ic,isha)*eair(p,ic) + delta(ic,isha)
       if (dpai(p,ic) > 0._r8) then
          pai = fracsun(p,ic) * dpai(p,ic)
          if (pai == 0._r8) tveg(p,ic,isun) = tair(p,ic)
          pai = fracsha(p,ic) * dpai(p,ic)
          if (pai == 0._r8) tveg(p,ic,isha) = tair(p,ic)
       else
          tveg(p,ic,isun) = tair(p,ic)
          tveg(p,ic,isha) = tair(p,ic)
       end if
    end do

    !---------------------------------------------------------------------
    ! Calculate ground fluxes
    !---------------------------------------------------------------------

    ! Ground temperature

    t0 = alpha0 * tair(p,1) + beta0 * eair(p,1) + delta0
    tveg(p,0,isun) = t0
    tveg(p,0,isha) = t0
    tair(p,0) = t0

    ! Water vapor at ground (mol/mol)

    eair(p,0) = rhg(p) * (qsat0 + dqsat0 * (t0 - tair_old(p,0)))

    ! Calculate fluxes

    sh0 = -cpair(p) * (tair(p,1) - t0) * ga_prof(p,0)
    et0 = -(eair(p,1) - eair(p,0)) * gs0
    g0 = c01 + c02 * t0

    !---------------------------------------------------------------------
    ! Convert water vapor from mol/mol to Pa
    !---------------------------------------------------------------------

    do ic = 0, ncan(p)
       eair(p,ic) = eair(p,ic) * pref(p)
    end do

    !---------------------------------------------------------------------
    ! Vertical sensible heat and water vapor fluxes between layers
    !---------------------------------------------------------------------

    do ic = 1, ncan(p)-1
       shair(p,ic) = -cpair(p) * (tair(p,ic+1) - tair(p,ic)) * ga_prof(p,ic)
       etair(p,ic) = -(eair(p,ic+1) - eair(p,ic)) / pref(p) * ga_prof(p,ic)
    end do
    ic = ncan(p)
    shair(p,ic) = -cpair(p) * (thref(p) - tair(p,ic)) * ga_prof(p,ic)
    etair(p,ic) = -(eref(p) - eair(p,ic)) / pref(p) * ga_prof(p,ic)

    !---------------------------------------------------------------------
    ! Canopy air storage flux (W/m2)
    !---------------------------------------------------------------------

    do ic = 1, ncan(p)
       den = (zw(p,ic) - zw(p,ic-1)) / dtime
       storage_sh(ic) = rhomol(p) * cpair(p) * (tair(p,ic) - tair_old(p,ic)) * den
       storage_et(ic) = rhomol(p) * (eair(p,ic) - eair_old(p,ic)) / pref(p) * den
       stair(p,ic) = storage_sh(ic) + storage_et(ic) * lambda
    end do

    !---------------------------------------------------------------------
    ! Canopy conservation checks
    !---------------------------------------------------------------------

    ! Canopy layer fluxes as calculated in this routine. Check energy balance

    do ic = 1, ncan(p)
       shsrc(ic) = 0._r8 ; etsrc(ic) = 0._r8 ; stveg(ic) = 0._r8
       do il = 1, nleaf
          shsrc(ic) = shsrc(ic) + cpair(p) * (tveg(p,ic,il) - tair(p,ic)) * gleaf_sh(ic,il)
          call SatVap (tveg_old(p,ic,il), esat, desat)
          etsrc(ic) = etsrc(ic) + (esat + desat * (tveg(p,ic,il) - tveg_old(p,ic,il)) - eair(p,ic)) / pref(p) &
                    * gleaf_et(ic,il)
          stveg(ic) = stveg(ic) + heatcap(ic,il) * (tveg(p,ic,il) - tveg_old(p,ic,il)) / dtime
       end do

       error = avail_energy(ic,isun) + avail_energy(ic,isha) - shsrc(ic) - lambda * etsrc(ic) - stveg(ic)
       if (abs(error) > 0.001_r8) then
          call endrun (msg=' ERROR: ScalarProfile: Leaf energy balance error')
       end if
    end do

    ! Flux conservation at each layer. Note special case for first canopy layer.

    do ic = 1, ncan(p)

       if (ic == 1) then
          error = storage_sh(ic) - (sh0 + shsrc(ic) - shair(p,ic))
       else
          error = storage_sh(ic) - (shair(p,ic-1) + shsrc(ic) - shair(p,ic))
       end if
       if (abs(error) > 0.001_r8) then
          call endrun (msg=' ERROR: ScalarProfile: Sensible heat layer conservation error')
       end if

       if (ic == 1) then
          error = storage_et(ic) - (et0 + etsrc(ic) - etair(p,ic))
       else
          error = storage_et(ic) - (etair(p,ic-1) + etsrc(ic) - etair(p,ic))
       end if
       error = error * lambda
       if (abs(error) > 0.001_r8) then
          call endrun (msg=' ERROR: ScalarProfile: Latent heat layer conservation error')
       end if

    end do

    ! Flux conservation for canopy sensible heat. This is to check canopy
    ! conservation equation (so the sum is to ntop not ncan).

    sum_src = 0._r8
    sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + shsrc(ic)
       sum_storage = sum_storage + storage_sh(ic)
    end do

    error = (sh0 + sum_src - sum_storage) - shair(p,ntop(p))
    if (abs(error) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Sensible heat canopy conservation error')
    end if

    ! Flux conservation for canopy latent heat. This is to check canopy
    ! conservation equation (so the sum is to ntop not ncan).

    sum_src = 0._r8
    sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + etsrc(ic)
       sum_storage = sum_storage + storage_et(ic)
    end do

    error = (et0 + sum_src - sum_storage) - etair(p,ntop(p))  ! mol H2O/m2/s
    error = error * lambda                                    ! W/m2
    if (abs(error) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Latent heat canopy conservation error')
    end if

    !---------------------------------------------------------------------
    ! Ground temperature energy balance conservation
    !---------------------------------------------------------------------

    error = rnsoi(p) - sh0 - lambda * et0 - g0
    if (abs(error) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Ground temperature energy balance error')
    end if

    !---------------------------------------------------------------------
    ! Use LeafTemperatureMod to calculate leaf fluxes for the current
    ! air temperature and vapor pressure profiles in the canopy. The
    ! flux and leaf temperature calculations there are the same as here, so
    ! the answers are the same in both routines. But remember that here
    ! the fluxes for sunlit/shaded leaves are multiplied by their leaf area.
    !---------------------------------------------------------------------

    if (.not. leaf_temp_iter) then

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0) then

             call LeafTemperature (p, ic, isun, mlcanopy_inst)
             call LeafTemperature (p, ic, isha, mlcanopy_inst)

             shsrc_leaf = (shleaf(p,ic,isun) * fracsun(p,ic) + shleaf(p,ic,isha) * fracsha(p,ic)) * dpai(p,ic)
             etsrc_leaf = (trleaf(p,ic,isun) + evleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &
                        + (trleaf(p,ic,isha) + evleaf(p,ic,isha)) * fracsha(p,ic) * dpai(p,ic)
             stveg_leaf = (stleaf(p,ic,isun) * fracsun(p,ic) + stleaf(p,ic,isha) * fracsha(p,ic)) * dpai(p,ic)

             ! Compare flux calculations between LeafTemperatureMod and here

             if (abs(shsrc(ic)-shsrc_leaf) > 0.001_r8) then
                call endrun (msg=' ERROR: ScalarProfile: Leaf sensible heat flux error')
             end if

             if (abs(lambda*etsrc(ic)-lambda*etsrc_leaf) > 0.001_r8) then
                call endrun (msg=' ERROR: ScalarProfile: Leaf latent heat flux error')
             end if

             if (abs(stveg(ic)-stveg_leaf) > 0.001_r8) then
                call endrun (msg=' ERROR: ScalarProfile: Leaf heat storage error')
             end if

             if (abs(tleaf(p,ic,isun)-tveg(p,ic,isun)) > 1.e-06_r8) then
                call endrun (msg=' ERROR: ScalarProfile: Leaf temperature error (sunlit)')
             end if

             if (abs(tleaf(p,ic,isha)-tveg(p,ic,isha)) > 1.e-06_r8) then
                call endrun (msg=' ERROR: ScalarProfile: Leaf temperature error (shaded)')
             end if

          end if
       end do

    end if

!   go to 100

    !---------------------------------------------------------------------
    ! Calculate soil fluxes
    !---------------------------------------------------------------------

    call SoilFluxesMultilayer (p, soilstate_inst, temperature_inst, &
    energyflux_inst, waterflux_inst, mlcanopy_inst)

    if (abs(shsoi(p)-sh0) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Soil sensible heat flux error')
    end if

    if (abs(lambda*etsoi(p)-lambda*et0) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Soil latent heat flux error')
    end if

    if (abs(gsoi(p)-g0) > 0.001_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Soil heat flux error')
    end if

    if (abs(tg(p)-t0) > 1.e-06_r8) then
       call endrun (msg=' ERROR: ScalarProfile: Soil temperature error')
    end if

100 continue

    end associate
  end subroutine ScalarProfile

  !-----------------------------------------------------------------------
  subroutine LookupPsihatINI
    !
    ! !DESCRIPTION:
    ! Initialize the look-up tables needed to calculate the RSL psihat functions
    !
    ! !USES:
    use clm_varcon, only : nZ, nL, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
    !
    ! !ARGUMENTS:
    implicit none
    !
    !LOCAL VARIABLES
    character(len=256) :: fin          ! File name
    integer :: nin                     ! Fortran unit number
    integer :: ii, jj                  ! Looping indices
    !---------------------------------------------------------------------

    fin = '../canopy/psihatM.dat'
    nin = shr_file_getUnit()
    open (unit=nin, file=trim(fin), action="read", status="old", err=10)
    read (nin,*,err=11) dtLgridM(1,1),(dtLgridM(1,jj), jj = 1,nL)
    do ii = 1, nZ
      read (nin,*,err=12) zdtgridM(ii,1),(psigridM(ii,jj), jj = 1,nL)
    end do
    close (nin)
    call shr_file_freeUnit (nin)
    
    fin = '../canopy/psihatH.dat'
    nin = shr_file_getUnit()
    open (unit=nin, file=trim(fin), action="read", status="old", err=20)
    read (nin,*,err=21) dtLgridH(1,1),(dtLgridH(1,jj), jj = 1,nL)
    do ii = 1,nZ
      read (nin,*,err=22) zdtgridH(ii,1),(psigridH(ii,jj), jj = 1,nL)
    end do
    close (nin)
    call shr_file_freeUnit (nin)

    return

10  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error opening psihatM.dat')

11  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading dtLgridM')

12  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading psigridM')

20  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error opening psihatH.dat')

21  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading dtLgridH')

22  continue
    call endrun (msg=' ERROR: CanopyTurbulenceMod: LookupPsihatINI error reading psigridH')

  end subroutine LookupPsihatINI

  !-----------------------------------------------------------------------
  subroutine ObuFuncCLM (p, iterate, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate the Obukhov length as in CLM
    !
    ! !USES:
    use clm_varcon, only : grav, vkc
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p       ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: iterate ! 1 = iterate for Obukhov length. 0 = use specified Obukhov length
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: iter               ! Iteration index
    real(r8) :: zldis              ! Atmospheric reference height minus displacement height, z-d (m)
    real(r8) :: z0mv               ! Roughness length over vegetation, momentum (m)
    real(r8) :: z0hv               ! Roughness length over vegetation, sensible heat (m)
    real(r8) :: egvf               ! Vegetation fraction of roughness length (dimensionless)
    real(r8) :: lt                 ! lai+sai for egvf (m2/m2)
    real(r8) :: obudif             ! Difference in Obukhov length between iterations (m)
    real(r8) :: obucur             ! Obukhov length from previous iteration (m)
    real(r8) :: zlog_m             ! Momentum      {log[(z-d)/z0m] - psi_m(zeta) + ...}
    real(r8) :: zlog_h             ! Sensible heat {log[(z-d)/z0h] - psi_h(zeta) + ...}
    real(r8) :: ur                 ! Wind speed at reference height (m/s)
    real(r8) :: dth                ! Difference of potential temperature between reference height and surface (K)
    real(r8) :: dqh                ! Difference of humidity between reference height and surface (kg/kg)
    real(r8) :: dthv               ! Difference of virtual potential temperature between reference height and surface (K)
    real(r8) :: um                 ! Wind speed at reference height including stability effect (m/s)
    real(r8) :: thvstar            ! Virtual potential temperature scale (K)
    real(r8) :: zeta               ! Monin-Obukhov stability parameter (z-d)/L
    real(r8) :: wc                 ! Convective velocity (m/s)
    real(r8) :: ram                ! Aerodynamic resistance for momentum (s/m)
    real(r8) :: rah                ! Aerodynamic resistance for sensible heat (s/m)
    real(r8) :: hcdis              ! Canopy height minus displacement height (m)

    integer,  parameter :: itermax = 40          ! Maximum number of iterations
    real(r8), parameter :: obutol = 0.1_r8       ! Convergence tolerance for Obukhov length (m)
    real(r8), parameter :: alpha_aero = 1.0_r8   ! Constant for aerodynamic parameter weighting
    real(r8), parameter :: tlsai_crit = 2.0_r8   ! Critical value of lai+sai for which aerodynamic parameters are maximum
    real(r8), parameter :: zii = 1000.0_r8       ! Convective boundary layer height (m)
    real(r8), parameter :: beta_con = 1.0_r8     ! Coefficient of convective velocity (-)
    !---------------------------------------------------------------------

    associate ( &
                                            ! *** Input ***
    z0mr    => pftcon%z0mr             , &  ! Ratio of momentum roughness length to canopy top height
    displar => pftcon%displar          , &  ! Ratio of displacement height to canopy top height
    ztop    => mlcanopy_inst%ztop      , &  ! Canopy height (m)
    lai     => mlcanopy_inst%lai       , &  ! Leaf area index of canopy (m2/m2)
    sai     => mlcanopy_inst%sai       , &  ! Stem area index of canopy (m2/m2)
    zref    => mlcanopy_inst%zref      , &  ! Reference height (m)
    uref    => mlcanopy_inst%uref      , &  ! Wind speed at reference height (m/s)
    qref    => mlcanopy_inst%qref      , &  ! Specific humidity at reference height (kg/kg)
    thref   => mlcanopy_inst%thref     , &  ! Atmospheric potential temperature (K)
    thvref  => mlcanopy_inst%thvref    , &  ! Atmospheric virtual potential temperature (K)
    rhomol  => mlcanopy_inst%rhomol    , &  ! Molar density at reference height (mol/m3)
    taf     => mlcanopy_inst%taf       , &  ! Air temperature at canopy top (K)
    qaf     => mlcanopy_inst%qaf       , &  ! Specific humidity at canopy top (kg/kg)
                                            ! *** Output ***
    z0mg    => mlcanopy_inst%z0mg      , &  ! Roughness length of ground (m)
    ustar   => mlcanopy_inst%ustar     , &  ! Friction velocity (m/s)
    tstar   => mlcanopy_inst%tstar     , &  ! Temperature scale (K)
    qstar   => mlcanopy_inst%qstar     , &  ! Water vapor scale (kg/kg)
    gah     => mlcanopy_inst%gah       , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu     => mlcanopy_inst%obu       , &  ! Obukhov length (m)
    obu_gah => mlcanopy_inst%obu_gah   , &  ! Obukhov length used for gah (m)
    PrSc    => mlcanopy_inst%PrSc      , &  ! Prandtl (Schmidt) number at canopy top
    zdisp   => mlcanopy_inst%zdisp     , &  ! Displacement height (m)
    uaf     => mlcanopy_inst%uaf         &  ! Wind speed at canopy top (m/s)
    )

    PrSc(p) = 1._r8

    ! Set roughness length of ground

    z0mg(p) = 0.01_r8

    ! Set roughness length of vegetation and displacement height

    zdisp(p) = displar(patch%itype(p)) * ztop(p)
    z0mv = z0mr(patch%itype(p)) * ztop(p)

    ! Modify aerodynamic parameters for sparse/dense canopy

    lt = min(lai(p)+sai(p), tlsai_crit)
    egvf = (1._r8 - alpha_aero * exp(-lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
    zdisp(p) = egvf * zdisp(p)
    z0mv = exp(egvf * log(z0mv) + (1._r8 - egvf) * log(z0mg(p)))
    z0hv = z0mv
    zldis = zref(p) - zdisp(p)

    ! Wind speed at reference height

    ur = max(1.0_r8, uref(p))

    ! Temperature and humidity differences between reference height and canopy

    dth = thref(p) - taf(p)
    dqh = qref(p) - qaf(p)
    dthv = dth * (1._r8 + 0.61_r8 * qref(p)) + 0.61_r8 * thref(p) * dqh
!   dthv = thvref(p) - taf(p)*(1._r8 + 0.61_r8 * qaf(p))                  ! Equivalent equation

    ! Initialize Obukhov length (obu) and wind speed (um)

    if (iterate == 1) then
       call MoninObukIniCLM (ur, thvref(p), dthv, zldis, z0mv, um, obu(p))
    else
       um = ur
    end if

    ! Iterative calculation for Obukhov length

    iter = 0
    obudif = 1.e36_r8

    do while (iter <= itermax .and. obudif > obutol)

       iter = iter + 1
       obucur = obu(p)

       ! Friction velocity and stability dependent log-z functions

       call FrictionVelocityCLM (zref(p), zdisp(p), z0mv, z0hv, obu(p), um, ustar(p), zlog_m, zlog_h)

       ! Above-canopy aerodynamic resistances

       ram = um / (ustar(p)*ustar(p))            ! Resistance: s/m
       rah = zlog_h / (vkc*ustar(p))             ! Resistance: s/m
       gah(p) = (1._r8 / rah) * rhomol(p)        ! Conductance: m/s -> mol/m2/s
       obu_gah(p) = obu(p)

       ! Wind speed at canopy height

       hcdis = ztop(p) - zdisp(p)
       if (zeta < zetam) then             ! very unstable (zeta < zetam)
          uaf(p) = um - ustar(p) / vkc * ( log(zetam*obu(p)/hcdis)      &
                 - StabilityFunc1(zetam) + StabilityFunc1(hcdis/obu(p)) &
                 + 1.14_r8 * ((-zeta)**0.333_r8 - (-zetam)**0.333_r8) )
       else if (zeta < 0._r8) then        ! unstable (zetam <= zeta < 0)
          uaf(p) = um - ustar(p) / vkc * ( log(zldis/hcdis)             &
                 - StabilityFunc1(zeta) + StabilityFunc1(hcdis/obu(p)) )
       else if (zeta <=  1._r8) then      ! stable (0 <= zeta <= 1)
          uaf(p) = um - ustar(p) / vkc * ( log(zldis/hcdis)             &
                 + 5._r8 * zeta - 5._r8 * hcdis/obu(p) )
       else                               ! very stable (zeta > 1)
          uaf(p) = um - ustar(p) / vkc * ( log(obu(p)/hcdis)            &
                 + 5._r8 - 5._r8 * hcdis / obu(p)                       &
                 + 5._r8 * log(zeta) + zeta - 1._r8 )
       end if

       ! Update Obukhov length and wind speed including the stability effect

       tstar(p) = dth * vkc / zlog_h
       qstar(p) = dqh * vkc / zlog_h
       thvstar = tstar(p) * (1._r8 + 0.61_r8 * qref(p)) + 0.61_r8 * thref(p) * qstar(p)
       zeta = zldis * vkc * grav * thvstar / (ustar(p)**2 * thvref(p))

       if (zeta >= 0._r8) then                       ! stable
          zeta = min(zetamaxstable, max(zeta,0.01_r8))
          um = max(ur, 0.1_r8)
       else                                          ! unstable
          zeta = max(-100._r8, min(zeta,-0.01_r8))
          wc = beta_con*(-grav*ustar(p)*thvstar*zii/thvref(p))**0.333_r8
          um = sqrt(ur*ur + wc*wc)
       end if
       obu(p) = zldis / zeta

       ! Check for convergence

       obudif = abs(obu(p)-obucur)

       ! Exit iteration after first pass if no iteration required

       if (iterate == 0) exit

    end do

    end associate
  end subroutine ObuFuncCLM

  !-----------------------------------------------------------------------
  subroutine MoninObukIniCLM (ur, thv, dthv, zldis, z0m, um, obu)
    !
    ! !DESCRIPTION:
    ! CLM initialization of the Obukhov length
    !
    ! !USES:
    use clm_varcon, only : grav
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: ur    ! Wind speed at reference height (m/s)
    real(r8), intent(in)  :: thv   ! Atmospheric virtual potential temperature (K)
    real(r8), intent(in)  :: dthv  ! Difference of virtual potential temperature between reference height and surface (K)
    real(r8), intent(in)  :: zldis ! Atmospheric reference height minus displacement height, z-d (m)
    real(r8), intent(in)  :: z0m   ! Roughness length, momentum (m)
    real(r8), intent(out) :: um    ! Wind speed at reference height including stability effect (m/s)
    real(r8), intent(out) :: obu   ! Obukhov length (m)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ustar              ! Friction velocity (m/s)
    real(r8) :: wc                 ! Convective velocity (m/s)
    real(r8) :: rib                ! Bulk Richardson number
    real(r8) :: zeta               ! Monin-Obukhov stability parameter (z-d)/L
    !---------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar = 0.06_r8
    wc = 0.5_r8

    if (dthv >= 0._r8) then
       um = max(ur, 0.1_r8)
    else
       um = sqrt(ur*ur + wc*wc)
    end if

    rib = grav * zldis * dthv / (thv * um * um)

    if (rib >= 0._r8) then
       ! neutral or stable
       zeta = rib * log(zldis/z0m) / (1._r8 - 5._r8 * min(rib,0.19_r8))
       zeta = min(zetamaxstable, max(zeta,0.01_r8))
    else
       ! unstable
       zeta = rib * log(zldis/z0m)
       zeta = max(-100._r8, min(zeta,-0.01_r8))
    end if

    obu = zldis / zeta

  end subroutine MoninObukIniCLM

  !-----------------------------------------------------------------------
  subroutine FrictionVelocityCLM (zref, displa, z0m, z0h, obu, um, ustar, zlog_m, zlog_h)
    !
    ! !DESCRIPTION:
    ! CLM friction velocity and stability dependent log-z functions
    !
    ! !USES:
    use clm_varcon, only : vkc
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: zref           ! Atmospheric reference height (m)
    real(r8), intent(in)  :: displa         ! Displacement height (m)
    real(r8), intent(in)  :: z0m            ! Roughness length over vegetation, momentum (m)
    real(r8), intent(in)  :: z0h            ! Roughness length over vegetation, sensible heat (m)
    real(r8), intent(in)  :: obu            ! Obukhov length (m)
    real(r8), intent(in)  :: um             ! Wind speed including the stablity effect (m/s)
    real(r8), intent(out) :: ustar          ! Friction velocity (m/s)
    real(r8), intent(out) :: zlog_m         ! Momentum      {log[(z-d)/z0m] - psi_m(zeta) + ...}
    real(r8), intent(out) :: zlog_h         ! Sensible heat {log[(z-d)/z0h] - psi_h(zeta) + ...}
    !
    ! !LOCAL VARIABLES:
    real(r8):: zldis                        ! Atmospheric reference height minus displacement height, z-d (m)
    real(r8):: zeta                         ! Monin-Obukhov stability parameter (z-d)/L
    !---------------------------------------------------------------------

    ! Monin-Obukhov "zeta"

    zldis = zref - displa
    zeta = zldis / obu

    ! Stability functions - psi_m and psi_h

     call GetPsiCLM (zldis, zeta, obu, z0m, z0h, zlog_m, zlog_h)

    ! Friction velocity

    ustar = vkc * um / zlog_m

  end subroutine FrictionVelocityCLM

  !-----------------------------------------------------------------------
  subroutine GetPsiCLM (zldis, zeta, obu, z0m, z0h, zlog_m, zlog_h)
    !
    ! !DESCRIPTION:
    ! CLM surface-layer stability functions
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: zldis           ! Atmospheric reference height minus displacement height, z-d (m)
    real(r8), intent(in)  :: zeta            ! Monin-Obukhov stability parameter (z-d)/L
    real(r8), intent(in)  :: obu             ! Obukhov length (m)
    real(r8), intent(in)  :: z0m             ! Roughness length for momentum (m)
    real(r8), intent(in)  :: z0h             ! Roughness length for sensible heat (m)
    real(r8), intent(out) :: zlog_m          ! Momentum      {log[(z-d)/z0m] - psi_m(zeta) + ...}
    real(r8), intent(out) :: zlog_h          ! Sensible heat {log[(z-d)/z0h] - psi_h(zeta) + ...}
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    ! Momentum

    if (zeta < zetam) then             ! very unstable (zeta < zetam)
       zlog_m = log(zetam*obu/z0m) - StabilityFunc1(zetam) + StabilityFunc1(z0m/obu) &
              + 1.14_r8 * ((-zeta)**0.333_r8 - (-zetam)**0.333_r8)
    else if (zeta < 0._r8) then        ! unstable (zetam <= zeta < 0)
       zlog_m = log(zldis/z0m) - StabilityFunc1(zeta) + StabilityFunc1(z0m/obu)
    else if (zeta <=  1._r8) then      ! stable (0 <= zeta <= 1)
       zlog_m = log(zldis/z0m) + 5._r8*zeta - 5._r8*z0m/obu
    else                               ! very stable (zeta > 1)
       zlog_m = log(obu/z0m) + 5._r8 - 5._r8*z0m/obu + 5._r8*log(zeta) + zeta - 1._r8
    end if
      
    ! Sensible heat and other scalar fluxes

    if (zeta < zetah) then             ! very unstable (zeta < zetah)
       zlog_h = log(zetah*obu/z0h) - StabilityFunc2(zetah) + StabilityFunc2(z0h/obu) &
              + 0.8_r8 * ((-zetah)**(-0.333_r8) - (-zeta)**(-0.333_r8))
    else if (zeta < 0._r8) then        ! unstable (zetah <= zeta < 0)
       zlog_h = log(zldis/z0h) - StabilityFunc2(zeta) + StabilityFunc2(z0h/obu)
    else if (zeta <=  1._r8) then      ! stable (0 <= zeta <= 1)
       zlog_h = log(zldis/z0h) + 5._r8*zeta - 5._r8*z0h/obu
    else                               ! very stable (zeta > 1)
       zlog_h = log(obu/z0h) + 5._r8 - 5._r8*z0h/obu + 5._r8*log(zeta) + zeta - 1._r8
    end if

  end subroutine GetPsiCLM

  !-----------------------------------------------------------------------
  function StabilityFunc1 (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! CLM Monin-Obukhov stability function for momentum
    !
    ! !USES:
    use clm_varcon, only : pi => rpi
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter (z-d)/L
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for momentum
    !---------------------------------------------------------------------

    x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
    psi = 2._r8 * log((1._r8+x)/2._r8) + log((1._r8+x*x)/2._r8) - 2._r8*atan(x) + pi * 0.5_r8

  end function StabilityFunc1

  !-----------------------------------------------------------------------
  function StabilityFunc2 (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! CLM Monin-Obukhov stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter (z-d)/L
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/2
    real(r8) :: psi               ! psi for scalars
    !---------------------------------------------------------------------

    x = sqrt(1._r8 - 16._r8 * zeta)
    psi = 2._r8 * log((1._r8+x)/2._r8)

  end function StabilityFunc2

end module CanopyTurbulenceMod

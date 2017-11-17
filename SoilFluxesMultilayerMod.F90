module SoilFluxesMultilayerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate soil surface temperature and energy balance
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilFluxesMultilayer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilFluxesMultilayer (p, soilstate_inst, temperature_inst, &
  energyflux_inst, waterflux_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Soil surface temperature and energy balance
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use clm_varcon, only : mmh2o, rgasc, grav, hvap
    use clm_varctl, only : iulog
    use WaterVaporMod, only : SatVap, LatVap
    use ColumnType, only : col
    use PatchType, only : patch
    use SoilStateType, only : soilstate_type
    use TemperatureType, only : temperature_type
    use EnergyFluxType, only : energyflux_type
    use WaterFluxType, only : waterflux_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p  ! Patch index for CLM g/l/c/p hierarchy
    type(soilstate_type), intent(inout) :: soilstate_inst
    type(temperature_type), intent(in) :: temperature_inst
    type(energyflux_type), intent(out) :: energyflux_inst
    type(waterflux_type), intent(out) :: waterflux_inst
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c                             ! Column index for CLM g/l/c/p hierarchy
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: gws                           ! Soil conductance for water vapor (mol H2O/m2/s)
    real(r8) :: gw                            ! Total conductance for water vapor (mol H2O/m2/s)
    real(r8) :: esat                          ! Saturation vapor pressure (Pa)
    real(r8) :: desat                         ! Derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                          ! Saturation vapor pressure of air (mol/mol)
    real(r8) :: dqsat                         ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: num1, num2, num3, num4, den   ! Intermediate terms
    real(r8) :: dshsoi                        ! Temperature derivative of sensible heat flux (W/m2/K)
    real(r8) :: dlhsoi                        ! Temperature derivative of latent heat flux (W/m2/K)
    real(r8) :: detsoi                        ! Temperature derivative of evaporation flux (mol H2O/m2/s/K)
    real(r8) :: err                           ! Surface energy imbalance (W/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                                ! *** Input ***
    swsoi          => mlcanopy_inst%swsoi                  , &  ! Absorbed solar radiation, ground (W/m2)
    irsoi          => mlcanopy_inst%irsoi                  , &  ! Absorbed longwave radiation, ground (W/m2)
    tref           => mlcanopy_inst%tref                   , &  ! Air temperature at reference height (K)
    pref           => mlcanopy_inst%pref                   , &  ! Air pressure at reference height (Pa)
    rhomol         => mlcanopy_inst%rhomol                 , &  ! Molar density at reference height (mol/m3)
    cpair          => mlcanopy_inst%cpair                  , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    tair           => mlcanopy_inst%tair                   , &  ! Air temperature profile (K)
    eair           => mlcanopy_inst%eair                   , &  ! Vapor pressure profile (Pa)
    ga_prof        => mlcanopy_inst%ga_prof                , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    tair_old       => mlcanopy_inst%tair_old               , &  ! Air temperature profile for previous timestep (K)
    z              => col%z                                , &  ! Soil layer depth (m)
    zi             => col%zi                               , &  ! Soil layer depth at layer interface (m)
    snl            => col%snl                              , &  ! Number of snow layers
    soilresis      => soilstate_inst%soilresis_col         , &  ! Soil evaporative resistance (s/m)
    thk            => soilstate_inst%thk_col               , &  ! Soil layer thermal conductivity (W/m/K)
    smp_l          => soilstate_inst%smp_l_col             , &  ! Soil layer matric potential (mm)
    t_soisno       => temperature_inst%t_soisno_col        , &  ! Soil temperature (K)
                                                                ! *** Output ***
    rnsoi          => mlcanopy_inst%rnsoi                  , &  ! Net radiation, ground (W/m2)
    shsoi          => mlcanopy_inst%shsoi                  , &  ! Sensible heat flux, ground (W/m2)
    lhsoi          => mlcanopy_inst%lhsoi                  , &  ! Latent heat flux, ground (W/m2)
    gsoi           => mlcanopy_inst%gsoi                   , &  ! Soil heat flux (W/m2)
    etsoi          => mlcanopy_inst%etsoi                  , &  ! Water vapor flux, ground (mol H2O/m2/s)
    tg             => mlcanopy_inst%tg                     , &  ! Soil surface temperature (K)
    eg             => mlcanopy_inst%eg                     , &  ! Soil surface vapor pressure (Pa)
    rhg            => mlcanopy_inst%rhg                    , &  ! Relative humidity of airspace at soil surface (fraction)
    eflx_sh_grnd   => energyflux_inst%eflx_sh_grnd_patch   , &  ! CLM: Sensible heat flux from ground (W/m2)
    eflx_sh_snow   => energyflux_inst%eflx_sh_snow_patch   , &  ! CLM: Sensible heat flux from snow (W/m2)
    eflx_sh_h2osfc => energyflux_inst%eflx_sh_h2osfc_patch , &  ! CLM: Sensible heat flux from surface water (W/m2)
    eflx_sh_soil   => energyflux_inst%eflx_sh_soil_patch   , &  ! CLM: Sensible heat flux from soil (W/m2)
    qflx_evap_soi  => waterflux_inst%qflx_evap_soi_patch   , &  ! CLM: Soil evaporation (kg/m2/s)
    qflx_ev_snow   => waterflux_inst%qflx_ev_snow_patch    , &  ! CLM: Evaporation flux from snow (kg/m2/s)
    qflx_ev_h2osfc => waterflux_inst%qflx_ev_h2osfc_patch  , &  ! CLM: Evaporation flux from h2osfc (kg/m2/s)
    qflx_ev_soil   => waterflux_inst%qflx_ev_soil_patch    , &  ! CLM: Evaporation flux from soil (kg/m2/s)
    cgrnds         => energyflux_inst%cgrnds_patch         , &  ! CLM: Deriv. of soil sensible heat flux wrt soil temp (W/m2/K)
    cgrndl         => energyflux_inst%cgrndl_patch         , &  ! CLM: Deriv. of soil evaporation flux wrt soil temp (kg/m2/s/K)
    cgrnd          => energyflux_inst%cgrnd_patch            &  ! CLM: Deriv. of soil energy flux wrt to soil temp (W/m2/K)
    )

    c = patch%column(p)

    ! Current ground temperature

    tg(p) = tair_old(p,0)

    ! Net radiation

    rnsoi(p) = swsoi(p,ivis) + swsoi(p,inir) + irsoi(p)

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    ! Relative humidity in soil airspace

    rhg(p) = exp(grav * mmh2o * smp_l(c,1)*1.e-03_r8 / (rgasc * t_soisno(c,1)))

    ! Soil conductance to water vapour diffusion

    gws = 1._r8 / soilresis(c)                       ! s/m -> m/s
    gws = gws * rhomol(p)                            ! m/s -> mol H2O/m2/s
    gw = ga_prof(p,0) * gws / (ga_prof(p,0) + gws)   ! total conductance

    ! Saturation vapor pressure at ground temperature (Pa -> mol/mol)

    call SatVap (tg(p), esat, desat)
    qsat = esat / pref(p) ; dqsat = desat / pref(p)

    ! Calculate soil surface temperature

    num1 = cpair(p) * ga_prof(p,0)
    num2 = lambda * gw
    num3 = thk(c,snl(c)+1) / (z(c,snl(c)+1)-zi(c,snl(c)))
    num4 = rnsoi(p) - num2 * rhg(p) * (qsat - dqsat * tg(p)) + num3 * t_soisno(c,snl(c)+1)
    den = num1 + num2 * dqsat * rhg(p) + num3
    tg(p) = (num1*tair(p,1) + num2*eair(p,1)/pref(p) + num4) / den

    ! Sensible heat flux

    shsoi(p) = cpair(p) * (tg(p) - tair(p,1)) * ga_prof(p,0)

    ! Latent heat flux - remember that tair_old(p,0) is tg(p) before the update for time n+1

    eg(p) = rhg(p) * (esat + desat * (tg(p) - tair_old(p,0)))
    lhsoi(p) = lambda / pref(p) * (eg(p) - eair(p,1)) * gw

    ! Soil heat flux

    gsoi(p) = thk(c,snl(c)+1) * (tg(p) - t_soisno(c,snl(c)+1)) / (z(c,snl(c)+1)-zi(c,snl(c)))
       
    ! Error check

    err = rnsoi(p) - shsoi(p) - lhsoi(p) - gsoi(p)
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: SoilFluxesMultilayerMod: energy balance error')
    end if

    ! Water vapor flux: W/m2 -> mol H2O/m2/s

    etsoi(p) = lhsoi(p) / lambda

    ! Needed for CLM soil temperature - Snow fluxes must be set equal to
    ! fluxes from the soil even if there is no snow

    eflx_sh_grnd(p) = shsoi(p)                ! Sensibe heat flux from ground (W/m2)
    eflx_sh_snow(p) = shsoi(p)                ! Sensible heat flux from snow (W/m2)
    eflx_sh_h2osfc(p) = 0._r8                 ! Sensible heat flux from surface water (W/m2)
    eflx_sh_soil(p) = shsoi(p)                ! Sensible heat flux from soil (W/m2)

    qflx_evap_soi(p) = etsoi(p) * mmh2o       ! Soil evaporation (kg/m2/s)
    qflx_ev_snow(p) = etsoi(p) * mmh2o        ! Evaporation flux from snow (kg/m2/s)
    qflx_ev_h2osfc(p) = 0._r8                 ! Evaporation flux from h2osfc (kg/m2/s)
    qflx_ev_soil(p) = etsoi(p) * mmh2o        ! Evaporation flux from soil (kg/m2/s)

    dshsoi = cpair(p) * ga_prof(p,0)
    dlhsoi  = lambda / pref(p) * gw * rhg(p) * desat
    detsoi = dlhsoi / lambda

    cgrnds(p) = dshsoi                        ! Temperature derivative of soil sensible heat flux (W/m2/K)
    cgrndl(p) = detsoi * mmh2o                ! Temperature deriative of soil evaporation flux (kg/m2/s/K)
    cgrnd(p) = cgrnds(p) + cgrndl(p)*hvap     ! Temperature derivative of soil energy flux (W/m2/K)

    end associate
  end subroutine SoilFluxesMultilayer

end module SoilFluxesMultilayerMod

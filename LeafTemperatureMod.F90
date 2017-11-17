module LeafTemperatureMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf temperature and energy fluxes (excluding evaporation of intercepted water)
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LeafTemperature       ! Leaf temperature and energy fluxes
  public :: LeafHeatCapacity      ! Leaf heat capacity
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafTemperature (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf temperature and energy fluxes
    !
    ! !USES:
    use clm_varctl, only : iulog, dtime_sub
    use WaterVaporMod, only : SatVap, LatVap
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p            ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic           ! Aboveground layer index
    integer, intent(in) :: il           ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime                   ! Model time step (s)
    real(r8) :: esat                    ! Saturation vapor pressure (Pa)
    real(r8) :: desat                   ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                    ! Saturation vapor pressure of air (mol/mol)
    real(r8) :: dqsat                   ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: lambda                  ! Latent heat of vaporization (J/mol)
    real(r8) :: gleaf                   ! Leaf conductance for transpiration (mol H2O/m2 leaf/s)
    real(r8) :: gw                      ! Total conductance including evaporation (mol H2O/m2 leaf/s)
    real(r8) :: num1, num2, num3, den   ! Intermediate calculation
    real(r8) :: err                     ! Energy balance error (W/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                ! *** Input ***
    dpai      => mlcanopy_inst%dpai        , &  ! Layer plant area index (m2/m2)
    tref      => mlcanopy_inst%tref        , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref        , &  ! Air pressure at reference height (Pa)
    cpair     => mlcanopy_inst%cpair       , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    tair      => mlcanopy_inst%tair        , &  ! Air temperature profile (K)
    eair      => mlcanopy_inst%eair        , &  ! Vapor pressure profile (Pa)
    gbh       => mlcanopy_inst%gbh         , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv         , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rnleaf    => mlcanopy_inst%rnleaf      , &  ! Leaf net radiation (W/m2 leaf)
    tleaf_old => mlcanopy_inst%tleaf_old   , &  ! Leaf temperature for previous timestep (K)
    cpleaf    => mlcanopy_inst%cpleaf      , &  ! Leaf heat capacity (J/m2 leaf/K)
    fwet      => mlcanopy_inst%fwet        , &  ! Fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry        , &  ! Fraction of plant area index that is green and dry
                                                ! *** Output ***
    tleaf     => mlcanopy_inst%tleaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    evleaf    => mlcanopy_inst%evleaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    trleaf    => mlcanopy_inst%trleaf        &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    ! Time step (s)

    dtime = dtime_sub

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       ! Saturation vapor pressure (Pa -> mol/mol)

       call SatVap (tleaf_old(p,ic,il), esat, desat)
       qsat = esat / pref(p) ; dqsat = desat / pref(p)

       ! Latent heat of vaporization

       lambda = LatVap(tref(p))

       ! Leaf conductance for transpiration

       gleaf = gs(p,ic,il) * gbv(p,ic,il) / (gs(p,ic,il) + gbv(p,ic,il))

       ! Total conductance (including evaporation)

       gw = gleaf * fdry(p,ic) + gbv(p,ic,il) * fwet(p,ic)

       ! Linearized leaf temperature calculation that balances the energy budget

       num1 = 2._r8 * cpair(p) * gbh(p,ic,il)
       num2 = lambda * gw 
       num3 = rnleaf(p,ic,il) - lambda * gw * (qsat - dqsat * tleaf_old(p,ic,il)) &
            + cpleaf(p,ic) / dtime * tleaf_old(p,ic,il)
       den = cpleaf(p,ic) / dtime + num1 + num2 * dqsat
       tleaf(p,ic,il) = (num1 * tair(p,ic) + num2 * eair(p,ic) / pref(p) + num3) / den

       ! Storage flux

       stleaf(p,ic,il) = (tleaf(p,ic,il) - tleaf_old(p,ic,il)) * cpleaf(p,ic) / dtime

       ! Sensible heat flux

       shleaf(p,ic,il) = 2._r8 * cpair(p) * (tleaf(p,ic,il) - tair(p,ic)) * gbh(p,ic,il)

       ! Transpiration and evaporation water fluxes: mol H2O/m2/s

       num1 = qsat + dqsat * (tleaf(p,ic,il) - tleaf_old(p,ic,il)) - eair(p,ic) / pref(p)
       trleaf(p,ic,il) = gleaf * fdry(p,ic) * num1
       evleaf(p,ic,il) = gbv(p,ic,il) * fwet(p,ic) * num1

       ! Latent heat flux

       lhleaf(p,ic,il) = (trleaf(p,ic,il) + evleaf(p,ic,il)) * lambda

       ! Error check

       err = rnleaf(p,ic,il) - shleaf(p,ic,il) - lhleaf(p,ic,il) - stleaf(p,ic,il)
       if (abs(err) > 1.e-03_r8) then
!         write (iulog,*) 'LeafTemperatureMod error:'
!         write (iulog,*) rnleaf(p,ic,il), shleaf(p,ic,il), lhleaf(p,ic,il), stleaf(p,ic,il), err
!         write (iulog,*) tleaf(p,ic,il), tleaf_old(p,ic,il)
          call endrun (msg=' ERROR: LeafTemperatureMod: energy balance error')
       end if

    else ! non-leaf layer

       tleaf(p,ic,il) = tair(p,ic)
       stleaf(p,ic,il) = 0._r8
       shleaf(p,ic,il) = 0._r8
       lhleaf(p,ic,il) = 0._r8
       evleaf(p,ic,il) = 0._r8
       trleaf(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine LeafTemperature

  !-----------------------------------------------------------------------
  subroutine LeafHeatCapacity (p, ic, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf heat capacity
    !
    ! !USES:
    use clm_varcon, only : cpliq, cpbio
    use clm_varctl, only : no_storage
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                 ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic                ! Aboveground layer index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: lma                          ! Leaf carbon mass per area (kg C / m2 leaf)
    real(r8) :: dry_weight                   ! Leaf dry mass per area (kg DM / m2 leaf)
    real(r8) :: fresh_weight                 ! Leaf fresh mass per area (kg FM / m2 leaf)
    real(r8) :: leaf_water                   ! Leaf water (kg H2O / m2 leaf)
    real(r8), parameter :: fcarbon = 0.5_r8  ! Fraction of dry biomass that is carbon
    real(r8), parameter :: fwater = 0.7_r8   ! Fraction of fresh biomass that is water
    !---------------------------------------------------------------------

    associate ( &
                                             ! *** Input ***
    slatop    => pftcon%slatop          , &  ! Specific leaf area at top of canopy (m2/gC)
    dpai      => mlcanopy_inst%dpai     , &  ! Layer plant area index (m2/m2)
                                             ! *** Output ***
    cpleaf    => mlcanopy_inst%cpleaf     &  ! Leaf heat capacity (J/m2 leaf/K)
    )

    ! Leaf heat capacity - need to convert specific leaf area (m2/gC) to
    ! leaf mass per area (kgC/m2) and then convert to dry weight (assume
    ! carbon is 50% of dry biomass). Then need to convert dry biomass to
    ! fresh biomass (assume 70% of fresh biomass is water). Then remember
    ! that 70% of fresh biomass is water when calculating heat capacity.

    if (dpai(p,ic) > 0._r8) then ! leaf layer
       lma = 1._r8 / slatop(patch%itype(p)) * 0.001_r8                ! m2 / g C -> kg C / m2
       dry_weight = lma / fcarbon                                     ! kg C / m2 -> kg DM / m2
       fresh_weight = dry_weight / (1._r8 - fwater)                   ! kg DM / m2 -> kg FM / m2
       leaf_water = fwater * fresh_weight                             ! Leaf water (kg H2O / m2 leaf)
       cpleaf(p,ic) = cpbio * dry_weight + cpliq * leaf_water         ! Heat capacity (J/K/m2 leaf) 
    else ! non-leaf layer
       cpleaf(p,ic) = 0._r8
    end if

    ! Use very low leaf heat capacity to reduce storage term

    if (no_storage) cpleaf(p,ic) = 1._r8

    end associate
  end subroutine LeafHeatCapacity

end module LeafTemperatureMod

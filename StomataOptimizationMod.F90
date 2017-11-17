module StomataOptimizationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
  ! using water-use efficiency optimization and cavitation check
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: StomataOptimization      ! Photosynthesis and stomatal conductance
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: StomataFluxes           ! Leaf calculations for a specificed stomatal conductance
  private :: StomataEfficiency       ! Water-use efficiency and cavitation checks for maximum gs
  private :: LeafWaterPotential      ! Leaf water potential from transpiration rate
  private :: LeafTranspiration       ! Leaf transpiration flux
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine StomataOptimization (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Photosynthesis and stomatal conductance with optimization
    !
    ! !USES:
    use clm_varctl, only : iulog
    use MathToolsMod, only : zbrent
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic             ! Aboveground layer index
    integer, intent(in) :: il             ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gs1, gs2                  ! Initial guess for gs (mol H2O/m2/s)
    real(r8) :: check1, check2            ! Water-use efficiency and cavitation check for gs1 and gs2
    real(r8), parameter :: tol = 0.004_r8 ! gs is updated to accuracy tol (mol H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
    dpai        => mlcanopy_inst%dpai    , &  ! Layer plant area index (m2/m2)
    lwp         => mlcanopy_inst%lwp     , &  ! Leaf water potential of canopy layer (MPa)
    psil        => mlcanopy_inst%psil    , &  ! Leaf water potential (MPa)
    gs          => mlcanopy_inst%gs        &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    ! Initialize leaf water potential (for sunlit or shaded leaf) to the layer
    ! value of the previous time step

    psil(p,ic,il) = lwp(p,ic)

    ! Low and high initial estimates for gs (mol H2O/m2/s)

    gs1 = 0.002_r8
    gs2 = 2._r8

    ! Calculate gs

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       ! Check for minimum stomatal conductance linked to low light or drought stress
       ! based on the water-use efficiency and cavitation checks for gs1 and gs2

       check1 = StomataEfficiency (p, ic, il, mlcanopy_inst, gs1)
       check2 = StomataEfficiency (p, ic, il, mlcanopy_inst, gs2)

       if (check1 * check2 < 0._r8) then

          ! Calculate gs using the function StomataEfficiency to iterate gs
          ! to an accuracy of tol (mol H2O/m2/s)

          gs(p,ic,il) = zbrent ('StomataOptimization', p, ic, il, mlcanopy_inst, StomataEfficiency, gs1, gs2, tol)

       else

          ! Low light or drought stress. Set gs to minimum conductance

          gs(p,ic,il) = 0.002_r8

       end if

    else ! non-leaf layer

       gs(p,ic,il) = 0._r8

    end if

    ! Leaf fluxes and leaf water potential for this gs

    call StomataFluxes (p, ic, il, mlcanopy_inst, gs(p,ic,il), psil(p,ic,il))

    end associate
  end subroutine StomataOptimization

  !-----------------------------------------------------------------------
  function StomataEfficiency (p, ic, il, mlcanopy_inst, gs_val) result(val)
    !
    ! !DESCRIPTION:
    ! Stomata water-use efficiency check and cavitation check to determine maximum gs. 
    ! For the stomatal conductance gs_val, calculate photosynthesis and leaf
    ! water potential for an increase in stomatal conductance equal to "delta".
    ! The returned value is positive if this increase produces a change in
    ! photosynthesis > iota*vpd*delta or if the leaf water potential is > minlwp.
    ! The returned value is negative if the increase produces a change in
    ! photosynthesis < iota*vpd*delta or if the leaf water potential is < minlwp. 
    !
    ! !USES:
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic       ! Aboveground layer index
    integer, intent(in)  :: il       ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: gs_val   ! Value for gs to use in calculations
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: delta                ! Small difference for gs (mol H2O/m2/s)
    real(r8) :: leafwp               ! Current leaf water potential (MPa)
    real(r8) :: gs2                  ! Lower value for gs (mol H2O/m2/s)
    real(r8) :: an2                  ! Leaf photosynthesis at gs2 (umol CO2/m2/s)
    real(r8) :: gs1                  ! Higher value for gs (mol H2O/m2/s)
    real(r8) :: an1                  ! Leaf photosynthesis at gs1 (umol CO2/m2/s)
    real(r8) :: wue                  ! Water-use efficiency check
    real(r8) :: minpsi               ! Cavitation check
    real(r8) :: val                  ! Returned minimum of the two checks
    !---------------------------------------------------------------------

    associate ( &
    minlwp    => pftcon%minlwp       , &  ! Minimum leaf water potential (MPa)
    iota      => pftcon%iota         , &  ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    psil      => mlcanopy_inst%psil  , &  ! Leaf water potential (MPa)
    an        => mlcanopy_inst%an    , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    vpd       => mlcanopy_inst%vpd   , &  ! Leaf vapor pressure deficit (Pa)
    pref      => mlcanopy_inst%pref    &  ! Air pressure at reference height (Pa)
    )

    ! Specify "delta" as a small difference in gs (mol H2O/m2/s)

    delta = 0.001_r8

    ! Photosynthesis at lower gs (gs_val - delta)

    leafwp = psil(p,ic,il)
    gs2 = gs_val - delta
    call StomataFluxes (p, ic, il, mlcanopy_inst, gs2, leafwp)
    an2 = an(p,ic,il)

    ! Photosynthesis at higher gs (gs_val)

    leafwp = psil(p,ic,il)
    gs1 = gs_val
    call StomataFluxes (p, ic, il, mlcanopy_inst, gs1, leafwp)
    an1 = an(p,ic,il)

    ! Efficiency check: wue < 0 when d(An) / d(gs) < iota * vpd

    wue = (an1 - an2) - iota(patch%itype(p)) * delta * (vpd(p,ic,il) / pref(p))

    ! Cavitation check: minpsi < 0 when leafwp < minlwp

    minpsi = leafwp - minlwp(patch%itype(p))

    ! Return the minimum of the two checks

    val = min(wue, minpsi)

    end associate
  end function StomataEfficiency

  !-----------------------------------------------------------------------
  subroutine StomataFluxes (p, ic, il, mlcanopy_inst, gs_val, leafwp)
    !
    ! !DESCRIPTION:
    ! Calculate leaf temperature, energy fluxes, photosynthesis, and leaf
    ! water potential for a specificed stomatal conductance (gs_val)
    !
    ! !USES:
    use clm_varctl, only : leaf_temp_iter
    use LeafTemperatureMod, only : LeafTemperature
    use LeafPhotosynthesisMod, only : LeafPhotosynthesis
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)     :: p         ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)     :: ic        ! Aboveground layer index
    integer, intent(in)     :: il        ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in)    :: gs_val    ! Value for gs to use in calculations
    real(r8), intent(inout) :: leafwp    ! Leaf water potential (MPa)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    associate (gs => mlcanopy_inst%gs)   ! Leaf stomatal conductance (mol H2O/m2 leaf/s)

    ! Use specified gs (gs_val)

    gs(p,ic,il) = gs_val

    ! Leaf temperature and energy fluxes

    if (leaf_temp_iter) then
       call LeafTemperature (p, ic, il, mlcanopy_inst)
    end if

    ! Leaf photosynthesis

    call LeafPhotosynthesis (p, ic, il, mlcanopy_inst)

    ! Leaf transpiration

    if (.not. leaf_temp_iter) then
       call LeafTranspiration (p, ic, il, mlcanopy_inst)
    end if

    ! Leaf water potential

    call LeafWaterPotential (p, ic, il, mlcanopy_inst, leafwp)

    end associate
  end subroutine StomataFluxes

  !-----------------------------------------------------------------------
  subroutine LeafWaterPotential (p, ic, il, mlcanopy_inst, leafwp)
    !
    ! !DESCRIPTION:
    ! Calculate leaf water potential for the current transpiration rate
    !
    ! !USES:
    use clm_varcon, only : denh2o, grav
    use clm_varctl, only : dtime_sub
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)     :: p              ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)     :: ic             ! Aboveground layer index
    integer, intent(in)     :: il             ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(inout) :: leafwp         ! Leaf water potential (MPa)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dtime                         ! Model time step (s)
    real(r8) :: y0                            ! Leaf water potential at beginning of timestep (MPa)
    real(r8) :: dy                            ! Change in leaf water potential (MPa)
    real(r8) :: a, b                          ! Intermediate calculation
    real(r8) :: head = denh2o*grav*1.e-06_r8  ! Head of pressure  (MPa/m)
    !---------------------------------------------------------------------

    associate ( &
    capac       => pftcon%capac          , &  ! Plant capacitance (mmol H2O/m2 leaf area/MPa)
    dpai        => mlcanopy_inst%dpai    , &  ! Layer plant area index (m2/m2)
    zs          => mlcanopy_inst%zs      , &  ! Canopy height for scalar concentration and source (m)
    psis        => mlcanopy_inst%psis    , &  ! Weighted soil water potential (MPa)
    lsc         => mlcanopy_inst%lsc     , &  ! Leaf-specific conductance for canopy layer (mmol H2O/m2 leaf/s/MPa)
    trleaf      => mlcanopy_inst%trleaf    &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    ! Get step size

    dtime = dtime_sub

    ! Change in leaf water potential is: dy / dt = (a - y) / b. The integrated change 
    ! over a full model timestep is: dy = (a - y0) * (1 - exp(-dt/b))

    if (dpai(p,ic) > 0._r8) then ! leaf layer
       y0 = leafwp
       a = psis(p) - head *  zs(p,ic) - 1000._r8 * trleaf(p,ic,il) / lsc(p,ic)
       b = capac(patch%itype(p)) / lsc(p,ic)
       dy = (a - y0) * (1._r8 - exp(-dtime/b))
       leafwp = y0 + dy
    else ! non-leaf layer
       leafwp = 0._r8
    end if

    end associate
  end subroutine LeafWaterPotential

  !-----------------------------------------------------------------------
  subroutine LeafTranspiration (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate leaf water transpiration flux
    !
    ! !USES:
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
    real(r8) :: esat                    ! Saturation vapor pressure (Pa)
    real(r8) :: desat                   ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: lambda                  ! Latent heat of vaporization (J/mol)
    real(r8) :: gleaf                   ! Leaf conductance for transpiration (mol H2O/m2 leaf/s)
    !---------------------------------------------------------------------

    associate ( &
                                                ! *** Input ***
    dpai      => mlcanopy_inst%dpai        , &  ! Layer plant area index (m2/m2)
    tref      => mlcanopy_inst%tref        , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref        , &  ! Air pressure at reference height (Pa)
    eair      => mlcanopy_inst%eair        , &  ! Vapor pressure profile (Pa)
    fdry      => mlcanopy_inst%fdry        , &  ! Fraction of plant area index that is green and dry
    tleaf     => mlcanopy_inst%tleaf       , &  ! Leaf temperature (K)
    gs        => mlcanopy_inst%gs          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv         , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
                                                ! *** Output ***
    trleaf    => mlcanopy_inst%trleaf        &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    )

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       ! Saturation vapor pressure

       call SatVap (tleaf(p,ic,il), esat, desat)

       ! Latent heat of vaporization

       lambda = LatVap(tref(p))

       ! Leaf conductance for transpiration

       gleaf = gs(p,ic,il) * gbv(p,ic,il) / (gs(p,ic,il) + gbv(p,ic,il))

       ! Transpiration flux: mol H2O/m2/s

       trleaf(p,ic,il) = gleaf * fdry(p,ic) * (esat - eair(p,ic)) / pref(p)

    else

       trleaf(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine LeafTranspiration

end module StomataOptimizationMod

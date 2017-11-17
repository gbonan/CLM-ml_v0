module PlantHydraulicsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate plant hydraulics
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilResistance         ! Calculate soil resistance and water uptake
  public :: PlantResistance        ! Calculate whole-plant resistance
  !-----------------------------------------------------------------------

contains

  subroutine PlantResistance (num_exposedvegp, filter_exposedvegp, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate whole-plant leaf-specific conductance (soil-to-leaf)
    ! 
    ! !USES:
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_exposedvegp       ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:) ! CLM patch filter for non-snow-covered vegetation
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                          ! Filter index
    integer  :: p                          ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                         ! Aboveground layer index
    real(r8) :: rplant                     ! Aboveground plant hydraulic resistance (MPa.s.m2/mmol H2O)
    !---------------------------------------------------------------------

    associate ( &
                                           ! *** Input ***
    gplant    => pftcon%gplant       , &   ! Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/MPa)
    ncan      => mlcanopy_inst%ncan  , &   ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai  , &   ! Layer plant area index (m2/m2)
    zs        => mlcanopy_inst%zs    , &   ! Canopy height for scalar concentration and source (m)
    rsoil     => mlcanopy_inst%rsoil , &   ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
                                           ! *** Output ***
    lsc       => mlcanopy_inst%lsc     &   ! Leaf-specific conductance of canopy layer (mmol H2O/m2 leaf/s/MPa)
    )

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then ! leaf layer

             ! Aboveground plant stem resistance, xylem-to-leaf (MPa.s.m2/mmol H2O)

!            rplant = zs(p,ic) / gplant(patch%itype(p))       ! gplant is conductivity (mmol/m/s/MPa)
             rplant = 1._r8 / gplant(patch%itype(p))          ! gplant is conductance (mmol/m2/s/MPa)

             ! Leaf specific conductance, soil-to-leaf (mmol H2O/m2/s/MPa)

             lsc(p,ic) = 1._r8 / (rsoil(p) + rplant)

          else ! non-leaf layer

             lsc(p,ic) = 0._r8

          end if

       end do
    end do

    end associate
  end subroutine PlantResistance

  !-----------------------------------------------------------------------
  subroutine SoilResistance (num_exposedvegp, filter_exposedvegp, &
  soilstate_inst, waterstate_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate soil hydraulic resistance and water uptake from each soil layer
    ! 
    ! !USES:
    use clm_varpar, only : nlevsoi
    use clm_varcon, only : pi => rpi, denh2o, grav, mmh2o
    use ColumnType, only : col
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use SoilStateType, only : soilstate_type
    use WaterStateType, only : waterstate_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_exposedvegp       ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:) ! CLM patch filter for non-snow-covered vegetation
    type(soilstate_type), intent(in) :: soilstate_inst
    type(waterstate_type), intent(in) :: waterstate_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: j                                ! Soil layer index
    real(r8) :: s                                ! Soil layer water content relative to saturation (fraction)
    real(r8) :: root_cross_sec_area              ! Root cross-sectional area (m2 root)
    real(r8) :: root_biomass_density             ! Root biomass density (g biomass / m3 soil) 
    real(r8) :: root_length_density              ! Root length density (m root / m3 soil) 
    real(r8) :: root_dist                        ! Mean distance between roots (m)
    real(r8) :: hk                               ! Hydraulic conductivity (mm/s -> mmol/m/s/MPa)
    real(r8) :: soilr1                           ! Soil-to-root resistance (MPa.s.m2/mmol H2O)
    real(r8) :: soilr2                           ! Root-to-stem resistance (MPa.s.m2/mmol H2O)
    real(r8) :: soilr                            ! Belowground resistance (MPa.s.m2/mmol H2O) 
    real(r8) :: smp_mpa(nlevsoi)                 ! Soil matric potential (MPa)
    real(r8) :: evap(nlevsoi)                    ! Maximum transpiration (mmol H2O/m2/s)
    real(r8) :: totevap                          ! Total maximum transpiration (mmol H2O/m2/s)

    real(r8), parameter :: head = denh2o*grav*1.e-06_r8  ! Head of pressure  (MPa/m)
    !---------------------------------------------------------------------

    associate ( &
                                                        ! *** Input ***
    minlwp       => pftcon%minlwp                  , &  ! Minimum leaf water potential (MPa)
    root_radius  => pftcon%root_radius             , &  ! Fine root radius (m)
    root_density => pftcon%root_density            , &  ! Fine root density (g biomass / m3 root)
    root_resist  => pftcon%root_resist             , &  ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
    dz           => col%dz                         , &  ! Soil layer thickness (m)
    watsat       => soilstate_inst%watsat_col      , &  ! Soil layer volumetric water content at saturation (porosity)
    hksat        => soilstate_inst%hksat_col       , &  ! Soil layer hydraulic conductivity at saturation (mm H2O/s)
    bsw          => soilstate_inst%bsw_col         , &  ! Soil layer Clapp and Hornberger "b" parameter
    smp_l        => soilstate_inst%smp_l_col       , &  ! Soil layer matric potential (mm)
    rootfr       => soilstate_inst%rootfr_patch    , &  ! Fraction of roots in each soil layer
    h2osoi_vol   => waterstate_inst%h2osoi_vol_col , &  ! Soil layer volumetric water content (m3/m3)
    h2osoi_ice   => waterstate_inst%h2osoi_ice_col , &  ! Soil layer ice lens (kg/m2)
    lai          => mlcanopy_inst%lai              , &  ! Leaf area index of canopy (m2/m2)
    root_biomass => mlcanopy_inst%root_biomass     , &  ! Fine root biomass (g biomass / m2)
                                                        ! *** Output ***
    psis         => mlcanopy_inst%psis             , &  ! Weighted soil water potential (MPa)
    rsoil        => mlcanopy_inst%rsoil            , &  ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    soil_et_loss => mlcanopy_inst%soil_et_loss       &  ! Fraction of total transpiration from each soil layer (-)
    )

    ! Soil and root resistances for each layer

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       c = patch%column(p)
       root_cross_sec_area = pi * root_radius(patch%itype(p))**2

       rsoil(p) = 0._r8
       do j = 1, nlevsoi

          ! Hydraulic conductivity and matric potential for each layer

          s = max(min(h2osoi_vol(c,j)/watsat(c,j), 1._r8), 0.01_r8)
          hk = hksat(c,j) * s**(2._r8 * bsw(c,j) + 3._r8)           ! mm/s
          hk = hk * 1.e-03_r8 / head                                ! mm/s -> m/s -> m2/s/MPa
          hk = hk * denh2o / mmh2o * 1000._r8                       ! m2/s/MPa -> mmol/m/s/MPa
          smp_mpa(j) = smp_l(c,j) * 1.e-03_r8 * head                ! mm -> m -> MPa

          ! Root biomass density: g biomass / m3 soil

          root_biomass_density = root_biomass(p) * rootfr(p,j) / dz(c,j)
          root_biomass_density = max(root_biomass_density, 1.e-10_r8)

          ! Root length density: m root per m3 soil

          root_length_density = root_biomass_density / (root_density(patch%itype(p)) * root_cross_sec_area)

          ! Distance between roots: m

          root_dist = sqrt (1._r8  / (root_length_density * pi))

          ! Soil-to-root resistance (MPa.s.m2/mmol H2O)

          soilr1 = log(root_dist/root_radius(patch%itype(p))) / (2._r8 * pi * root_length_density * dz(c,j) * hk)

          ! Root-to-stem resistance (MPa.s.m2/mmol H2O)

          soilr2 = root_resist(patch%itype(p)) / (root_biomass_density * dz(c,j))

          ! Belowground resistance (MPa.s.m2/mmol H2O) 

          soilr = soilr1 + soilr2

          ! Total belowground resistance. First sum the conductances (1/soilr)
          ! for each soil layer and then convert back to a resistance after the
          ! summation.

          rsoil(p) = rsoil(p) + 1._r8 / soilr

          ! Maximum transpiration for each layer (mmol H2O/m2/s). No negative
          ! transpiration and no transpiration from frozen soil.

          evap(j) = (smp_mpa(j) - minlwp(patch%itype(p))) / soilr
          evap(j) = max (evap(j), 0._r8)
          if (h2osoi_ice(c,j) > 0._r8) evap(j) = 0._r8

       end do

       ! Belowground resistance: resistance = 1 / conductance

       rsoil(p) = lai(p) / rsoil(p)

       ! Weighted soil water potential (MPa) and fractional uptake from soil layers

       totevap = sum(evap)
       psis(p) = 0._r8
       soil_et_loss(p,:) = 0._r8

       do j = 1, nlevsoi
          psis(p) = psis(p) + smp_mpa(j) * evap(j)
          if (totevap > 0._r8) then
             soil_et_loss(p,j) = evap(j) / totevap
          else
             soil_et_loss(p,j) = 1._r8 / nlevsoi
          end if
       end do

       if (totevap > 0._r8) then
          psis(p) = psis(p) / totevap
       else
          psis(p) = minlwp(patch%itype(p))
       end if

    end do

    end associate
  end subroutine SoilResistance

end module PlantHydraulicsMod

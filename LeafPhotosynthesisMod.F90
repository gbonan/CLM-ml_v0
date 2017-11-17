module LeafPhotosynthesisMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate leaf photosynthesis and stomatal conductance
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  !
  ! Leaf-level photosynthesis variables
  !
  ! *** Input parameters ***
  !
  private
  real(r8) :: kc25        ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
  real(r8) :: ko25        ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
  real(r8) :: cp25        ! CO2 compensation point at 25C (umol/mol)

  real(r8) :: vcmaxha     ! Activation energy for vcmax (J/mol)
  real(r8) :: jmaxha      ! Activation energy for jmax (J/mol)
  real(r8) :: rdha        ! Activation energy for rd (J/mol)
  real(r8) :: kcha        ! Activation energy for kc (J/mol)
  real(r8) :: koha        ! Activation energy for ko (J/mol)
  real(r8) :: cpha        ! Activation energy for cp (J/mol)

  real(r8) :: vcmaxhd     ! Deactivation energy for vcmax (J/mol)
  real(r8) :: jmaxhd      ! Deactivation energy for jmax (J/mol)
  real(r8) :: rdhd        ! Deactivation energy for rd (J/mol)

  real(r8) :: vcmaxse     ! Entropy term for vcmax (J/mol/K)
  real(r8) :: jmaxse      ! Entropy term for jmax (J/mol/K)
  real(r8) :: rdse        ! Entropy term for rd (J/mol/K)

  real(r8) :: vcmaxc      ! Scaling factor for high temperature inhibition (25 C = 1.0)
  real(r8) :: jmaxc       ! Scaling factor for high temperature inhibition (25 C = 1.0)
  real(r8) :: rdc         ! Scaling factor for high temperature inhibition (25 C = 1.0)

  real(r8) :: qe_c4       ! Quantum yield, used only for C4 (mol CO2 / mol photons)
  real(r8) :: phi_psii    ! Quantum yield of PS II
  real(r8) :: theta_j     ! Empirical curvature parameter for electron transport rate

  real(r8) :: vpd_min     ! Minimum vapor pressure deficit for Medlyn stomatal conductance (Pa)
  !
  ! *** Calculated variables ***
  !
  real(r8) :: vcmax       ! Maximum carboxylation rate (umol/m2/s)
  real(r8) :: jmax        ! Maximum electron transport rate (umol/m2/s)
  real(r8) :: je          ! Electron transport rate (umol/m2/s)
  real(r8) :: kc          ! Michaelis-Menten constant for CO2 (umol/mol)
  real(r8) :: ko          ! Michaelis-Menten constant for O2 (mmol/mol)
  real(r8) :: cp          ! CO2 compensation point (umol/mol)
  real(r8) :: rdleaf      ! Leaf respiration rate (umol CO2/m2 leaf/s)
  real(r8) :: kp          ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
  real(r8) :: g0          ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
  real(r8) :: g1          ! Ball-Berry slope of conductance-photosynthesis relationship

  real(r8) :: ceair       ! Vapor pressure of air, constrained (Pa)
  real(r8) :: esat        ! Saturation vapor pressure (Pa)
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: PhotosynthesisParam ! Leaf-level parameters for photosynthesis model
  public :: LeafPhotosynthesis  ! Leaf photosynthesis and stomatal conductance
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CiFunc             ! Function used to iterate Ci calculation
  private :: CiFuncGs           ! Function to evaluate Ci, with gs specified as input
  private :: ft                 ! Photosynthesis temperature response
  private :: fth                ! Photosynthesis temperature inhibition
  private :: fth25              ! Scaling factor for photosynthesis temperature inhibition
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine PhotosynthesisParam (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf-level parameters for photosynthesis model
    !
    ! !USES:
    use clm_varctl, only : use_acclim
    use clm_varcon, only : tfrz
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p  ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: sco           ! Relative specificity of rubisco
    !---------------------------------------------------------------------

    associate ( &
    o2ref       => mlcanopy_inst%o2ref       , &  ! Atmospheric O2 at reference height (mmol/mol)
    tacclim     => mlcanopy_inst%tacclim       &  ! Average air temperature for acclimation (K)
    )

    !---------------------------------------------------------------------
    ! kc, ko, cp at 25C: Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    ! Derive sco from cp with o2=0.209 mol/mol and re-calculate cp to allow
    ! variation in o2
    !---------------------------------------------------------------------

    kc25 = 404.9_r8                 ! umol/mol
    ko25 = 278.4_r8                 ! mmol/mol
    cp25 = 42.75_r8                 ! umol/mol

    sco = 0.5_r8 * 0.209_r8 / (cp25 * 1.e-06_r8) ! cp25 (umol/mol) -> (mol/mol)
    cp25 = 0.5_r8 * o2ref(p) / sco * 1000._r8    ! O2 is mmol/mol. Multiply by 1000 for umol/mol

    !---------------------------------------------------------------------
    ! Activation energy:
    ! Bernacchi et al (2001) Plant, Cell Environment 24:253-259
    ! Bernacchi et al (2003) Plant, Cell Environment 26:1419-1430
    !
    ! Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    !---------------------------------------------------------------------

    kcha    = 79430._r8
    koha    = 36380._r8
    cpha    = 37830._r8
    vcmaxha = 65330._r8
    jmaxha  = 43540._r8
    rdha    = 46390._r8

    if (use_acclim) then
       vcmaxha = 72000._r8
       jmaxha  = 50000._r8
    end if

    !---------------------------------------------------------------------
    ! High temperature deactivation: 
    ! Leuning (2002) Plant, Cell Environment 25:1205-1210
    ! The factor "c" scales the deactivation to a value of 1.0 at 25C
    !
    ! Acclimation from: Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190
    !---------------------------------------------------------------------

    vcmaxhd = 150000._r8
    jmaxhd  = 150000._r8
    rdhd    = 150000._r8

    if (use_acclim) then
       vcmaxhd = 200000._r8
       jmaxhd  = 200000._r8
    end if

    vcmaxse = 490._r8
    jmaxse  = 490._r8
    rdse    = 490._r8

    if (use_acclim) then
       vcmaxse = 668.39_r8 - 1.07_r8 * min(max((tacclim(p)-tfrz),11._r8),35._r8)
       jmaxse  = 659.70_r8 - 0.75_r8 * min(max((tacclim(p)-tfrz),11._r8),35._r8)
    end if

    vcmaxc = fth25 (vcmaxhd, vcmaxse)
    jmaxc  = fth25 (jmaxhd, jmaxse)
    rdc    = fth25 (rdhd, rdse)

    !---------------------------------------------------------------------
    ! Miscellaneous parameters
    !---------------------------------------------------------------------

    qe_c4 = 0.05_r8
    phi_psii = 0.70_r8
!   phi_psii = 0.85_r8
    theta_j = 0.90_r8

    vpd_min = 100._r8 ! Minimum vapor pressure deficit for Medlyn stomatal conductance (Pa)

    end associate
  end subroutine PhotosynthesisParam 

  !-----------------------------------------------------------------------
  subroutine LeafPhotosynthesis (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use clm_varctl, only : gstyp, iulog
    use WaterVaporMod, only : SatVap
    use MathToolsMod, only : hybrid, quadratic
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic       ! Aboveground layer index
    integer, intent(in) :: il       ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: qabs                ! PAR utilized by PS II (umol photons/m2/s)
    real(r8) :: desat               ! Derivative of saturation vapor pressure (Pa/K)
    real(r8) :: gs_err              ! gs for error check
    real(r8) :: an_err              ! An for error check
    real(r8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(r8) :: r1,r2               ! Roots of quadratic equation
    real(r8) :: ci0, ci1            ! Initial estimates for Ci
    real(r8), parameter :: tol = 0.1_r8 ! Convergence tolerance for Ci (mol/mol)
    !---------------------------------------------------------------------

    associate ( &
                                             ! *** Input ***
    c3psn     => pftcon%c3psn           , &  ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
    g0opt     => pftcon%g0opt           , &  ! Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
    g1opt     => pftcon%g1opt           , &  ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    dpai      => mlcanopy_inst%dpai     , &  ! Layer plant area index (m2/m2)
    btran     => mlcanopy_inst%btran    , &  ! Ball-Berry soil wetness factor (-)
    vcmax25   => mlcanopy_inst%vcmax25  , &  ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    jmax25    => mlcanopy_inst%jmax25   , &  ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    rd25      => mlcanopy_inst%rd25     , &  ! Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
    kp25      => mlcanopy_inst%kp25     , &  ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)
    eair      => mlcanopy_inst%eair     , &  ! Vapor pressure profile (Pa)
    cair      => mlcanopy_inst%cair     , &  ! Atmospheric CO2 profile (umol/mol)
    tleaf     => mlcanopy_inst%tleaf    , &  ! Leaf temperature (K)
    gbv       => mlcanopy_inst%gbv      , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc      , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    apar      => mlcanopy_inst%apar     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
                                             ! *** Output ***
    rd        => mlcanopy_inst%rd       , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci        => mlcanopy_inst%ci       , &  ! Leaf intercellular CO2 (umol/mol)
    hs        => mlcanopy_inst%hs       , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd       => mlcanopy_inst%vpd      , &  ! Leaf vapor pressure deficit (Pa)
                                             ! *** Output from calls to CiFunc or CiFuncGs ***
    ac        => mlcanopy_inst%ac       , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj        => mlcanopy_inst%aj       , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap        => mlcanopy_inst%ap       , &  ! Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag        => mlcanopy_inst%ag       , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an        => mlcanopy_inst%an       , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs        => mlcanopy_inst%cs       , &  ! Leaf surface CO2 (umol/mol)
    gs        => mlcanopy_inst%gs         &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       !------------------------------------------------------------------
       ! Adjust photosynthetic and conductance parameters for temperature
       ! and soil water
       !------------------------------------------------------------------

       ! C3 temperature response

       kc     = kc25          * ft(tleaf(p,ic,il), kcha)
       ko     = ko25          * ft(tleaf(p,ic,il), koha)
       cp     = cp25          * ft(tleaf(p,ic,il), cpha)
       vcmax  = vcmax25(p,ic) * ft(tleaf(p,ic,il), vcmaxha) * fth(tleaf(p,ic,il), vcmaxhd, vcmaxse, vcmaxc) 
       jmax   = jmax25(p,ic)  * ft(tleaf(p,ic,il), jmaxha)  * fth(tleaf(p,ic,il), jmaxhd, jmaxse, jmaxc)
       rdleaf = rd25(p,ic)    * ft(tleaf(p,ic,il), rdha)    * fth(tleaf(p,ic,il), rdhd, rdse, rdc)

       ! C4 temperature response

       if (nint(c3psn(patch%itype(p))) == 0) then
          vcmax  = vcmax25(p,ic) * 2.0**( (tleaf(p,ic,il)-(tfrz+25._r8)) / 10._r8 ) 
          vcmax  = vcmax / ( 1._r8 + exp( 0.2_r8*((tfrz+15._r8)-tleaf(p,ic,il)) ) )
          vcmax  = vcmax / ( 1._r8 + exp( 0.3_r8*(tleaf(p,ic,il)-(tfrz+40._r8)) ) ) 
          rdleaf = rd25(p,ic) * 2.0**( (tleaf(p,ic,il)-(tfrz+25._r8)) / 10._r8 ) 
          rdleaf = rdleaf / ( 1._r8 + exp( 1.3_r8*(tleaf(p,ic,il)-(tfrz+55._r8)) ) )
       end if
       kp = kp25(p,ic) * 2.0**( (tleaf(p,ic,il)-(tfrz+25._r8)) / 10._r8 ) 

       ! Soil water

       if (gstyp == 0) then
          vcmax = vcmax * btran(p)
          g0 = g0opt(patch%itype(p))
       end if

       if (gstyp == 1) then
          vcmax = vcmax * btran(p)
          g0 = max( g0opt(patch%itype(p)) * btran(p), 1.e-06_r8 )
       end if

       g1 = g1opt(patch%itype(p))

       ! Save leaf respiration

       rd(p,ic,il) = rdleaf

       !------------------------------------------------------------------
       ! Saturation vapor pressure at leaf temperature
       !------------------------------------------------------------------

       call SatVap (tleaf(p,ic,il), esat, desat)

       ! Constrain eair >= 0.05*esat[tleaf] so that solution does not blow up. This 
       ! ensures that hs does not go to zero. Also eair <= esat[tleaf] so that hs <= 1. 

!      ceair = min( max(eair(p,ic), 0.05_r8*esat), esat )
       ceair = min( max(eair(p,ic), 0.20_r8*esat), esat )

       !------------------------------------------------------------------
       ! Electron transport rate for C3 plants
       !------------------------------------------------------------------

       qabs = 0.5_r8 * phi_psii * apar(p,ic,il)
       aquad = theta_j
       bquad = -(qabs + jmax)
       cquad = qabs * jmax
       call quadratic (aquad, bquad, cquad, r1, r2)
       je = min(r1,r2)

       !------------------------------------------------------------------
       ! Ci calculation
       !------------------------------------------------------------------

       if (gstyp <= 1) then

          ! Initial estimates for Ci

          if (nint(c3psn(patch%itype(p))) == 1) then
             ci0 = 0.7_r8 * cair(p,ic)
          else
             ci0 = 0.4_r8 * cair(p,ic)
          end if
          ci1 = ci0 * 0.99_r8

          ! Solve for Ci: Use CiFunc to iterate photosynthesis calculations
          ! until the change in Ci is < tol. Ci has units umol/mol

          ci(p,ic,il) = hybrid ('LeafPhotosynthesis', p, ic, il, mlcanopy_inst, CiFunc, ci0, ci1, tol)

       else if (gstyp == 2) then

          ! Calculate photosynthesis for a specified stomatal conductance

          ci(p,ic,il) = CiFuncGs (p, ic, il, mlcanopy_inst)

       end if

       !------------------------------------------------------------------
       ! Make sure iterative solution is correct
       !------------------------------------------------------------------

       if (gs(p,ic,il) < 0._r8) then
          call endrun (msg=' ERROR: LeafPhotosynthesisMod: negative stomatal conductance')
       end if

       ! Compare with Ball-Berry model: gs = g1 * An * hs/cs + g0
       ! Use hs calculated with ceair

       if (gstyp == 1) then
          hs(p,ic,il) = (gbv(p,ic,il)*ceair + gs(p,ic,il)*esat) / ((gbv(p,ic,il)+gs(p,ic,il))*esat)
          gs_err = g1*max(an(p,ic,il), 0._r8)*hs(p,ic,il)/cs(p,ic,il) + g0
          if (abs(gs(p,ic,il)-gs_err)*1.e06_r8 > 1.e-04_r8) then
             call endrun (msg=' ERROR: LeafPhotosynthesisMod: failed Ball-Berry error check')
          end if
       end if

       ! Compare with Medlyn model: gs = g0 + 1.6 * (1 + g1 / sqrt(Ds)) * An / cs
       ! Use Ds calculated with ceair and also (esat - ceair) > vpd_min

       if (gstyp == 0) then
          ceair = min(ceair, esat - vpd_min)
          hs(p,ic,il) = (gbv(p,ic,il)*ceair + gs(p,ic,il)*esat) / ((gbv(p,ic,il)+gs(p,ic,il))*esat)
          vpd(p,ic,il) = esat - hs(p,ic,il) * esat
          gs_err = g0 + 1.6_r8 * (1._r8 + g1 / sqrt(vpd(p,ic,il)*0.001_r8)) * max(an(p,ic,il),0._r8)/cs(p,ic,il)
          if (abs(gs(p,ic,il)-gs_err)*1.e06_r8 > 1._r8) then
             call endrun (msg=' ERROR: LeafPhotosynthesisMod: failed Medlyn error check')
          end if
       end if

       ! Compare with diffusion equation: An = (ca - ci) * gleaf

       an_err = (cair(p,ic) - ci(p,ic,il)) / (1._r8 / gbc(p,ic,il) + 1.6_r8 / gs(p,ic,il))
       if (an(p,ic,il) > 0._r8 .and. abs(an(p,ic,il)-an_err) > 0.01_r8) then
          call endrun (msg=' ERROR: LeafPhotosynthesisMod: failed diffusion error check')
       end if

       !------------------------------------------------------------------
       ! Relative humidity and vapor pressure at leaf surface
       !------------------------------------------------------------------

       hs(p,ic,il) = (gbv(p,ic,il)*eair(p,ic) + gs(p,ic,il)*esat) / ((gbv(p,ic,il)+gs(p,ic,il))*esat)
       vpd(p,ic,il) = max(esat - hs(p,ic,il)*esat, 0.1_r8)

    else ! non-leaf layer

       rd(p,ic,il) = 0._r8
       if (gstyp <= 1) then
          ci(p,ic,il) = CiFunc (p, ic, il, mlcanopy_inst, 0._r8)
       else if (gstyp == 2) then
          ci(p,ic,il) = CiFuncGs (p, ic, il, mlcanopy_inst)
       end if
       hs(p,ic,il) = 0._r8
       vpd(p,ic,il) = 0._r8

   end if

    end associate
  end subroutine LeafPhotosynthesis

  !-----------------------------------------------------------------------
  function CiFunc (p, ic, il, mlcanopy_inst, ci_val) result(ci_dif)
    !
    ! !DESCRIPTION:
    ! Calculate leaf photosynthesis and stomatal conductance for a specified Ci
    ! (ci_val).  Then calculate a new Ci from the diffusion equation. This
    ! function equals zero when Ci has converged to the value that satisfies
    ! the metabolic, stomatal constraint, and diffusion equations.
    !
    ! !USES:
    use clm_varctl, only : use_colim, gstyp
    use MathToolsMod, only : quadratic
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p       ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic      ! Aboveground layer index
    integer, intent(in)  :: il      ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: ci_val  ! Input value for Ci (umol/mol)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(r8) :: r1,r2               ! Roots of quadratic equation
    real(r8) :: ai                  ! Intermediate co-limited photosynthesis (umol CO2/m2/s)
    real(r8) :: gleaf               ! Leaf CO2 conductance (mol CO2/m2/s)
    real(r8) :: cinew               ! New value for Ci
    real(r8) :: ci_dif              ! Difference in Ci
    real(r8) :: term                ! Term for Medlyn stomatal conductance
    real(r8) :: vpd_term            ! Vapor pressure deficit for Medlyn stomatal conductance (kPa)
    !---------------------------------------------------------------------

    associate ( &
                                          ! *** Input ***
    c3psn  => pftcon%c3psn           , &  ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant 
    dpai   => mlcanopy_inst%dpai     , &  ! Layer plant area index (m2/m2)
    o2ref  => mlcanopy_inst%o2ref    , &  ! Atmospheric O2 at reference height (mmol/mol)
    cair   => mlcanopy_inst%cair     , &  ! Atmospheric CO2 profile (umol/mol)
    gbv    => mlcanopy_inst%gbv      , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gbc    => mlcanopy_inst%gbc      , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    apar   => mlcanopy_inst%apar     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
                                          ! *** Output ***
    ac     => mlcanopy_inst%ac       , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj     => mlcanopy_inst%aj       , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap     => mlcanopy_inst%ap       , &  ! Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag     => mlcanopy_inst%ag       , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an     => mlcanopy_inst%an       , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs     => mlcanopy_inst%cs       , &  ! Leaf surface CO2 (umol/mol)
    gs     => mlcanopy_inst%gs         &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       !------------------------------------------------------------------
       ! Metabolic (demand-based) photosynthetic rate
       !------------------------------------------------------------------

       if (nint(c3psn(patch%itype(p))) == 1) then

          ! C3: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax * max(ci_val-cp, 0._r8) / (ci_val + kc*(1._r8+o2ref(p)/ko))
 
          ! C3: RuBP-limited photosynthesis
          aj(p,ic,il) = je * max(ci_val-cp, 0._r8) / (4._r8*ci_val + 8._r8*cp)

          ! C3: Product-limited photosynthesis
          ap(p,ic,il) = 0._r8

       else

          ! C4: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax
 
          ! C4: RuBP-limited photosynthesis
          aj(p,ic,il) = qe_c4 * apar(p,ic,il)

          ! C4: PEP carboxylase-limited (CO2-limited)
          ap(p,ic,il) = kp * max(ci_val, 0._r8)

       end if

       ! Net photosynthesis as the minimum or co-limited rate

       if (use_colim) then

          ! First co-limit ac and aj

          if (nint(c3psn(patch%itype(p))) == 1) then
             aquad = 0.98_r8
          else
             aquad = 0.80_r8
          end if
          bquad = -(ac(p,ic,il) + aj(p,ic,il))
          cquad = ac(p,ic,il) * aj(p,ic,il)
          call quadratic (aquad, bquad, cquad, r1, r2)
          ai = min(r1,r2)

          ! Now co-limit using ap, but only for C4 plants

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = ai
          else
             aquad = 0.95_r8
             bquad = -(ai + ap(p,ic,il))
             cquad = ai * ap(p,ic,il)
             call quadratic (aquad, bquad, cquad, r1, r2)
             ag(p,ic,il) = min(r1,r2)
          end if

       else if (.not. use_colim) then

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il))
          else
             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il),ap(p,ic,il))
          end if

       end if

       ! Prevent photosynthesis from ever being negative

       ac(p,ic,il) = max(ac(p,ic,il), 0._r8)
       aj(p,ic,il) = max(aj(p,ic,il), 0._r8)
       ap(p,ic,il) = max(ap(p,ic,il), 0._r8)
       ag(p,ic,il) = max(ag(p,ic,il), 0._r8)

       ! Net photosynthesis

       an(p,ic,il) = ag(p,ic,il) - rdleaf

       !------------------------------------------------------------------
       ! CO2 at leaf surface
       !------------------------------------------------------------------

       cs(p,ic,il) = cair(p,ic) - an(p,ic,il) / gbc(p,ic,il)
       cs(p,ic,il) = max(cs(p,ic,il), 1._r8)

       !------------------------------------------------------------------
       ! Stomatal constraint function
       !------------------------------------------------------------------

       ! Ball-Berry stomatal conductance
       ! Quadratic gs calculation given An. Valid for An >= 0. With An <= 0, gs = g0
 
       if (gstyp == 1) then
          if (an(p,ic,il) > 0._r8) then
             aquad = cs(p,ic,il)
             bquad = cs(p,ic,il)*(gbv(p,ic,il) - g0) - g1*an(p,ic,il)
             cquad = -gbv(p,ic,il) * (cs(p,ic,il)*g0 + g1*an(p,ic,il)*ceair/esat)
             call quadratic (aquad, bquad, cquad, r1, r2)
             gs(p,ic,il) = max(r1,r2)
          else
             gs(p,ic,il) = g0
          end if
       end if

       ! Medlyn stomatal conductance
       ! Quadratic gs calculation given An. Valid for An >= 0. With An <= 0, gs = g0.
       ! Note that vapor pressure deficit is limited to be > vpd_min

       if (gstyp == 0) then
          if (an(p,ic,il) > 0._r8) then
             vpd_term = max((esat - ceair), vpd_min) * 0.001_r8
             term = 1.6_r8 * an(p,ic,il) / cs(p,ic,il)
             aquad = 1._r8
             bquad = -(2._r8 * (g0 + term) + (g1 * term)**2 / (gbv(p,ic,il) * vpd_term))
             cquad = g0 * g0 + (2._r8 * g0 + term * (1._r8 - g1 * g1 / vpd_term)) * term
             call quadratic (aquad, bquad, cquad, r1, r2)
             gs(p,ic,il) = max(r1,r2)
          else
             gs(p,ic,il) = g0
          end if
       end if

       !------------------------------------------------------------------
       ! Diffusion (supply-based) photosynthetic rate - Calculate Ci
       ! from the diffusion rate
       !------------------------------------------------------------------

       gleaf = 1._r8 / (1._r8/gbc(p,ic,il) + 1.6_r8/gs(p,ic,il))
       cinew = cair(p,ic) - an(p,ic,il) / gleaf

       !------------------------------------------------------------------
       ! CiFunc is the difference between the current Ci and the new Ci
       !------------------------------------------------------------------

       ci_dif = cinew - ci_val
       if (an(p,ic,il) < 0._r8) ci_dif = 0._r8

    else ! non-leaf layer

       ac(p,ic,il) = 0._r8
       aj(p,ic,il) = 0._r8
       ap(p,ic,il) = 0._r8
       ag(p,ic,il) = 0._r8
       an(p,ic,il) = 0._r8
       cs(p,ic,il) = 0._r8
       gs(p,ic,il) = 0._r8
       ci_dif = 0._r8

    end if

    end associate
  end function CiFunc

  !-----------------------------------------------------------------------
  function CiFuncGs (p, ic, il, mlcanopy_inst) result(ci_val)
    !
    ! !DESCRIPTION:
    ! Calculate leaf photosynthesis for a specified stomatal conductance.
    ! Then calculate Ci from the diffusion equation. 
    !
    ! !USES:
    use clm_varctl, only : use_colim
    use MathToolsMod, only : quadratic
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic       ! Aboveground layer index
    integer, intent(in) :: il       ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: gleaf               ! Leaf CO2 conductance (mol CO2/m2/s)
    real(r8) :: a0,e0,d0            ! Terms for quadratic photosynthesis calculation
    real(r8) :: aquad,bquad,cquad   ! Terms for quadratic equations
    real(r8) :: r1,r2               ! Roots of quadratic equation
    real(r8) :: ai                  ! Intermediate co-limited photosynthesis (umol CO2/m2/s)
    real(r8) :: ci_val              ! Calculated value for Ci (umol/mol)
    !---------------------------------------------------------------------

    associate ( &
                                          ! *** Input ***
    c3psn  => pftcon%c3psn           , &  ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant 
    dpai   => mlcanopy_inst%dpai     , &  ! Layer plant area index (m2/m2)
    o2ref  => mlcanopy_inst%o2ref    , &  ! Atmospheric O2 at reference height (mmol/mol)
    cair   => mlcanopy_inst%cair     , &  ! Atmospheric CO2 profile (umol/mol)
    gbc    => mlcanopy_inst%gbc      , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    gs     => mlcanopy_inst%gs       , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    apar   => mlcanopy_inst%apar     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
                                          ! *** Output ***
    ac     => mlcanopy_inst%ac       , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj     => mlcanopy_inst%aj       , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap     => mlcanopy_inst%ap       , &  ! Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag     => mlcanopy_inst%ag       , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an     => mlcanopy_inst%an       , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    cs     => mlcanopy_inst%cs         &  ! Leaf surface CO2 (umol/mol)
    )

    !---------------------------------------------------------------------
    ! Calculate leaf photosynthesis for a specified stomatal conductance.
    ! Then calculate Ci from the diffusion equation. 
    !
    ! This routine uses a quadratic equation to solve for net photosynthesis (An).
    ! A general equation for C3 photosynthesis is:
    !
    !      a (Ci - Cp)
    ! An = ----------- - Rd
    !        e Ci + d
    !
    ! where:
    !
    ! An = Net leaf photosynthesis (umol CO2/m2/s)
    ! Rd = Leaf respiration (umol CO2/m2/s)
    ! Ci = Intercellular CO2 concentration (umol/mol)
    ! Cp = CO2 compensation point (umol/mol)
    ! 
    ! Rubisco-limited photosynthesis (Ac)
    ! a  = Vcmax
    ! e  = 1
    ! d  = Kc (1 + Oi/Ko)
    !
    ! RuBP-limited photosynthesis (Aj)
    ! a = J
    ! e = 4
    ! d = 8 Cp
    !
    ! where:
    !
    ! Vcmax = Maximum carboxylation rate (umol/m2/s)
    ! Kc    = Michaelis-Menten constant for CO2 (umol/mol)
    ! Ko    = Michaelis-Menten constant for O2 (mmol/mol)
    ! Oi    = Intercellular O2 concentration (mmol/mol)
    ! J     = Electron transport rate (umol/m2/s)
    !
    ! Ci is calculated from the diffusion equation:
    !
    !                   1.4   1.6
    ! An = (Ca - Ci) / (--- + ---)
    !                   gb    gs
    !
    !            1.4   1.6
    ! Ci = Ca - (--- + ---) An
    !            gb    gs
    !
    ! where:
    ! 
    ! Ca  = Atmospheric CO2 concentration (umol/mol)
    ! gb  = Leaf boundary layer conductance (mol H2O/m2/s)
    ! gs  = Leaf stomatal conductance (mol H2O/m2/s)
    ! 1.4 = Corrects gb for the diffusivity of CO2 compared with H2O
    ! 1.6 = Corrects gs for the diffusivity of CO2 compared with H2O
    !
    ! The resulting quadratic equation is: a An**2 + b An + c = 0
    !
    ! A similar approach is used for C4 photosynthesis
    !---------------------------------------------------------------------

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       ! Leaf conductance: gbc has units mol CO2/m2/s, gs has units mol H2O/m2/s,
       ! gleaf has units mol CO2/m2/s

       gleaf = 1._r8 / (1._r8/gbc(p,ic,il) + 1.6_r8/gs(p,ic,il))

       !------------------------------------------------------------------
       ! Gross assimilation rates
       !------------------------------------------------------------------

       if (nint(c3psn(patch%itype(p))) == 1) then

          ! C3: Rubisco-limited photosynthesis

          a0 = vcmax
          e0 = 1._r8
          d0 = kc * (1._r8 + o2ref(p) / ko)

          aquad = e0 / gleaf
          bquad = -(e0*cair(p,ic) + d0) - (a0 - e0*rdleaf) / gleaf
          cquad = a0 * (cair(p,ic) - cp) - rdleaf * (e0*cair(p,ic) + d0)

          call quadratic (aquad, bquad, cquad, r1, r2)
          ac(p,ic,il) = min(r1,r2) + rdleaf
 
          ! C3: RuBP-limited photosynthesis

          a0 = je
          e0 = 4._r8
          d0 = 8._r8 * cp

          aquad = e0 / gleaf
          bquad = -(e0*cair(p,ic) + d0) - (a0 - e0*rdleaf) / gleaf
          cquad = a0 * (cair(p,ic) - cp) - rdleaf * (e0*cair(p,ic) + d0)

          call quadratic (aquad, bquad, cquad, r1, r2)
          aj(p,ic,il) = min(r1,r2) + rdleaf

          ! C3: Product-limited photosynthesis

          ap(p,ic,il) = 0._r8

       else

          ! C4: Rubisco-limited photosynthesis
          ac(p,ic,il) = vcmax
 
          ! C4: RuBP-limited photosynthesis
          aj(p,ic,il) = qe_c4 * apar(p,ic,il)

          ! C4: PEP carboxylase-limited (CO2-limited)
          ap(p,ic,il) = kp * (cair(p,ic) * gleaf + rdleaf) / (gleaf + kp)

       end if

       !------------------------------------------------------------------
       ! Net assimilation as the minimum or co-limited rate
       !------------------------------------------------------------------

       if (use_colim) then

          ! First co-limit ac and aj

          if (nint(c3psn(patch%itype(p))) == 1) then
             aquad = 0.98_r8
          else
             aquad = 0.80_r8
          end if
          bquad = -(ac(p,ic,il) + aj(p,ic,il))
          cquad = ac(p,ic,il) * aj(p,ic,il)
          call quadratic (aquad, bquad, cquad, r1, r2)
          ai = min(r1,r2)

          ! Now co-limit using ap, but only for C4 plants

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = ai
          else
             aquad = 0.95_r8
             bquad = -(ai + ap(p,ic,il))
             cquad = ai * ap(p,ic,il)
             call quadratic (aquad, bquad, cquad, r1, r2)
             ag(p,ic,il) = min(r1,r2)
          end if

       else if (.not. use_colim) then

          if (nint(c3psn(patch%itype(p))) == 1) then
             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il))
          else
             ag(p,ic,il) = min(ac(p,ic,il),aj(p,ic,il),ap(p,ic,il))
          end if

       end if
 
       an(p,ic,il) = ag(p,ic,il) - rdleaf

       !------------------------------------------------------------------
       ! Leaf surface CO2
       !------------------------------------------------------------------

       cs(p,ic,il) = cair(p,ic) - an(p,ic,il) / gbc(p,ic,il)

       !------------------------------------------------------------------
       ! Intercelluar CO2
       !------------------------------------------------------------------

       ci_val = cair(p,ic) - an(p,ic,il) / gleaf

    else

       ac(p,ic,il) = 0._r8
       aj(p,ic,il) = 0._r8
       ap(p,ic,il) = 0._r8
       ag(p,ic,il) = 0._r8
       an(p,ic,il) = 0._r8
       cs(p,ic,il) = 0._r8
       ci_val = 0._r8

    end if

    end associate
  end function CiFuncGs

  !-----------------------------------------------------------------------
  subroutine C13Fractionation (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! 13C fractionation for photosynthesis
    !
    ! !USES:
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p        ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in) :: ic       ! Aboveground layer index
    integer, intent(in) :: il       ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    associate ( &
                                             ! *** Input ***
    c3psn     => pftcon%c3psn           , &  ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
    dpai      => mlcanopy_inst%dpai     , &  ! Layer plant area index (m2/m2)
    apar      => mlcanopy_inst%apar     , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    ci        => mlcanopy_inst%ci       , &  ! Leaf intercellular CO2 (umol/mol)
    cair      => mlcanopy_inst%cair     , &  ! Atmospheric CO2 profile (umol/mol)
                                             ! *** Output ***
    alphapsn => mlcanopy_inst%alphapsn    &  ! Leaf 13C fractionation factor for photosynthesis (-)
    )

    if (dpai(p,ic) > 0._r8) then ! leaf layer

       if (apar(p,ic,il) > 0._r8) then ! day
          if (nint(c3psn(patch%itype(p))) == 1) then ! C3
             alphapsn(p,ic,il) = 1._r8 + (4.4_r8 + 22.6_r8 * ci(p,ic,il) / cair(p,ic)) / 1000._r8
          else ! C4
             alphapsn(p,ic,il) = 1._r8 + 4.4_r8 / 1000._r8
          end if
       else ! night
          alphapsn(p,ic,il) = 1._r8
       end if

    else ! non-leaf layer

       alphapsn(p,ic,il) = 0._r8

    end if

    end associate
  end subroutine C13Fractionation

  !-----------------------------------------------------------------------
  function ft (tl, ha) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature response
    !
    ! !USES:
    use clm_varcon, only : rgasc, tfrz
    !
    ! !ARGUMENTS:
    implicit none
    real(r8) :: tl                   ! Leaf temperature (K)
    real(r8) :: ha                   ! Activation energy (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = exp( ha / (rgasc * (tfrz+25._r8)) * (1._r8 - (tfrz+25._r8) / tl) )

  end function ft

  !-----------------------------------------------------------------------
  function fth (tl, hd, se, c) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature inhibition
    !
    ! !USES:
    use clm_varcon, only : rgasc
    !
    ! !ARGUMENTS:
    implicit none
    real(r8) :: tl                   ! Leaf temperature (K)
    real(r8) :: hd                   ! Deactivation energy (J/mol)
    real(r8) :: se                   ! Entropy term (J/mol/K)
    real(r8) :: c                    ! Scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = c / ( 1._r8 + exp( (-hd + se*tl) / (rgasc*tl) ) )

  end function fth

  !-----------------------------------------------------------------------
  function fth25 (hd, se) result(ans)
    !
    ! !DESCRIPTION:
    ! Scaling factor for photosynthesis temperature inhibition
    !
    ! !USES:
    use clm_varcon, only : rgasc, tfrz
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: hd       ! Deactivation energy (J/mol)
    real(r8), intent(in) :: se       ! Entropy term (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd + se * (tfrz+25._r8)) / (rgasc * (tfrz+25._r8)) )

  end function fth25

end module LeafPhotosynthesisMod

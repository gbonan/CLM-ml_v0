module CanopyNitrogenProfileMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy profile of nitrogen and photosynthetic capacity
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyNitrogenProfile
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine CanopyNitrogenProfile (num_exposedvegp, filter_exposedvegp, &
  mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy profile of nitrogen and photosynthetic capacity
    !
    ! !USES:
    use clm_varctl, only : use_clm45kn
    use clm_varctl, only : use_acclim
    use clm_varcon, only : tfrz
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
    integer  :: f               ! Filter index
    integer  :: p               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic              ! Aboveground layer index
    real(r8) :: vcmax25top      ! Canopy top - Maximum carboxylation rate at 25C (umol/m2/s)
    real(r8) :: jmax25top       ! Canopy top - C3: Maximum electron transport rate at 25C (umol/m2/s)
    real(r8) :: rd25top         ! Canopy top - Leaf respiration rate at 25C (umol CO2/m2/s)
    real(r8) :: kp25top         ! Canopy top - C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
    real(r8) :: kn              ! Leaf nitrogen decay coefficient
    real(r8) :: nscale          ! Nitrogen scaling coefficient
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    vcmaxpft    => pftcon%vcmaxpft             , &  ! Maximum carboxylation rate at 25C (umol/m2/s)
    c3psn       => pftcon%c3psn                , &  ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
    tacclim     => mlcanopy_inst%tacclim       , &  ! Average air temperature for acclimation (K)
    ncan        => mlcanopy_inst%ncan          , &  ! Number of aboveground layers
    dpai        => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    sumpai      => mlcanopy_inst%sumpai        , &  ! Cumulative plant area index (m2/m2)
                                                    ! *** Output ***
    vcmax25     => mlcanopy_inst%vcmax25       , &  ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    jmax25      => mlcanopy_inst%jmax25        , &  ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    rd25        => mlcanopy_inst%rd25          , &  ! Leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
    kp25        => mlcanopy_inst%kp25            &  ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)
    )

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! vcmax and other parameters (at 25C and top of canopy). jmax acclimation
       ! from Kattge and Knorr (2007) Plant, Cell Environment 30:1176-1190

       vcmax25top = vcmaxpft(patch%itype(p))

       if (nint(c3psn(patch%itype(p))) == 1) then
          if (use_acclim) then
             jmax25top = (2.59_r8 - 0.035_r8*min(max((tacclim(p)-tfrz),11._r8),35._r8)) * vcmax25top
          else
             jmax25top = 1.67_r8 * vcmax25top
          end if
          rd25top = 0.015_r8 * vcmax25top
          kp25top = 0._r8
       else
          jmax25top = 0._r8
          rd25top = 0.025_r8 * vcmax25top
          kp25top = 0.02_r8 * vcmax25top
       end if

       ! Leaf nitrogen decay coefficient

       if (use_clm45kn) then
          kn = 0.3_r8
       else
          kn = exp(0.00963_r8 * vcmax25top - 2.43_r8)
       end if

       ! Layer values

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then ! leaf layer
             nscale = exp (-kn * sumpai(p,ic))
             vcmax25(p,ic) = vcmax25top * nscale
             jmax25(p,ic) = jmax25top * nscale
             rd25(p,ic) = rd25top * nscale
             kp25(p,ic) = kp25top * nscale
          else ! non-leaf leaf layer
             vcmax25(p,ic) = 0._r8
             jmax25(p,ic) = 0._r8
             rd25(p,ic) = 0._r8
             kp25(p,ic) = 0._r8
          end if
       end do

    end do

    end associate
  end subroutine CanopyNitrogenProfile

end module CanopyNitrogenProfileMod

module LeafFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LeafFluxes         ! Leaf flux calculations
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: TleafFunc         ! Function to evaluate leaf temperature
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafFluxes (p, ic, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
    !
    ! !USES:
    use clm_varctl, only : gstyp, leaf_temp_iter
    use StomataOptimizationMod, only : StomataOptimization
    use MathToolsMod, only : hybrid
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
    real(r8) :: t0, t1                  ! Initial estimates for leaf temperature (K)
    real(r8), parameter :: tol=0.1_r8   ! Accuracy tolerance for tleaf (K)
    real(r8) :: dummy                   ! Dummy argument, not needed
    !---------------------------------------------------------------------

    associate ( &
    dpai      => mlcanopy_inst%dpai  , &  ! Layer plant area index (m2/m2)
    tleaf     => mlcanopy_inst%tleaf , &  ! Leaf temperature (K)
    tair      => mlcanopy_inst%tair    &  ! Air temperature profile (K)
    )

    if (dpai(p,ic) > 0._r8) then

       ! Calculate fluxes for leaf layers using TleafFunc for Ball-Berry style
       ! stomatal model or StomataOptimization for water-use efficiency
       ! optimization model. These routines calculate leaf temperature and
       ! stomatal conductance simultaneously. If leaf_temp_iter = false, the leaf
       ! temperature calculation is turned off and the stomatal conductance
       ! routines use leaf temperature from the previous sub-timestep.

       if (gstyp <= 1) then

          ! Initial estimates for leaf temperature

          t0 = tleaf(p,ic,il) - 1._r8
          t1 = tleaf(p,ic,il) + 1._r8

          ! Solve for tleaf: Use TleafFunc to iterate leaf temperature, energy fluxes,
          ! photosynthesis and stomatal conductance. This temperature is refined to an
          ! accuracy of tol. Do not use the returned value (dummy), and instead use
          ! the tleaf calculated in the final call to TleafFunc.

          dummy = hybrid ('LeafFluxes', p, ic, il, mlcanopy_inst, TleafFunc, t0, t1, tol)

       else if (gstyp == 2) then

          ! Iterate leaf temperature, flux calculations, and stomatal conductance
          ! using water-use efficiency optimization and cavitation check

          call StomataOptimization (p, ic, il, mlcanopy_inst)

       end if

    else ! non-leaf layer

       ! Zero out fluxes

       if (gstyp <= 1) then
          dummy = TleafFunc (p, ic, il, mlcanopy_inst, tair(p,ic))
       else if (gstyp == 2) then
          call StomataOptimization (p, ic, il, mlcanopy_inst)
       end if

    end if

    end associate
  end subroutine LeafFluxes

  !-----------------------------------------------------------------------
  function TleafFunc (p, ic, il, mlcanopy_inst, tleaf_val) result(tleaf_dif)
    !
    ! !DESCRIPTION:
    ! Calculate leaf fluxes for an input leaf temperature (tleaf_val) and
    ! compare the new temperature to the prior temperature. This function
    ! equals zero when tleaf does not change between iterations.
    !
    ! !USES:
    use clm_varctl, only : iulog, leaf_temp_iter
    use LeafPhotosynthesisMod, only : LeafPhotosynthesis
    use LeafTemperatureMod, only : LeafTemperature
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p           ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic          ! Aboveground layer index
    integer, intent(in)  :: il          ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: tleaf_val   ! Input value for tleaf
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: tleaf_dif               ! Difference in tleaf
    !---------------------------------------------------------------------

    associate (tleaf => mlcanopy_inst%tleaf)  ! Leaf temperature (K)

    if (tleaf_val < 0._r8) then
       call endrun (msg=' ERROR: LeafFluxesMod: TleafFunc error')
    end if

    tleaf(p,ic,il) = tleaf_val

    ! Leaf photosynthesis and stomatal conductance

    call LeafPhotosynthesis (p, ic, il, mlcanopy_inst)

    ! Leaf temperature and energy fluxes

    if (leaf_temp_iter) then
       call LeafTemperature (p, ic, il, mlcanopy_inst)
    end if

    ! Compare with prior value for leaf temperature

    tleaf_dif = tleaf(p,ic,il) - tleaf_val

    end associate
  end function TleafFunc

end module LeafFluxesMod

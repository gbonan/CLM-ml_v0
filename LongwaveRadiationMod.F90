module LongwaveRadiationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate longwave radiation transfer through canopy
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LongwaveRadiation
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LongwaveRadiation (bounds, num_exposedvegp, filter_exposedvegp, &
  mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Longwave radiation transfer through canopy using Norman (1979)
    !
    ! !USES:
    use clm_varpar, only : isun, isha, nlevcan
    use clm_varcon, only : sb
    use clm_varctl, only : iulog
    use decompMod, only : bounds_type
    use pftconMod, only : pftcon
    use PatchType, only : patch
    use MathToolsMod, only : tridiag
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp        ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)  ! CLM patch filter for non-snow-covered vegetation
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                             ! Filter index
    integer  :: p                             ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                            ! Aboveground layer index
    integer  :: icm1                          ! Layer below ic (ic-1)
    real(r8) :: emg                           ! Ground (soil) emissivity
    real(r8) :: sumabs                        ! Absorbed radiation for energy conservation check
    real(r8) :: error                         ! Error check
    real(r8) :: omega                         ! Leaf scattering coefficient
    real(r8) :: rho                           ! Leaf reflectance
    real(r8) :: tau                           ! Leaf transmittance
    real(r8) :: trand                         ! Term for longwave radiation transmitted by layer
    real(r8) :: refld                         ! Term for longwave radiation reflected by layer
    real(r8) :: ir_source_sun                 ! Longwave radiation emitted by sunlit leaf (W/m2)
    real(r8) :: ir_source_sha                 ! Longwave radiation emitted by shaded leaf (W/m2)
    real(r8) :: ir_source(nlevcan)            ! Longwave radiation emitted by leaf layer (W/m2)
    integer  :: m                             ! Index to the tridiagonal matrix
    real(r8) :: aic, bic                      ! Intermediate terms for tridiagonal matrix
    real(r8) :: eic, fic                      ! Intermediate terms for tridiagonal matrix
    integer, parameter :: neq = (nlevcan+1)*2 ! Number of tridiagonal equations to solve
    real(r8) :: atri(neq), btri(neq)          ! Entries in tridiagonal matrix
    real(r8) :: ctri(neq), dtri(neq)          ! Entries in tridiagonal matrix
    real(r8) :: utri(neq)                     ! Tridiagonal solution
    real(r8) :: irabs                         ! Absorbed longwave flux (W/m2 ground)
    real(r8) :: irup(bounds%begp:bounds%endp,0:nlevcan) ! Upward longwave flux above canopy layer (W/m2 ground)
    real(r8) :: irdn(bounds%begp:bounds%endp,0:nlevcan) ! Downward longwave flux onto canopy layer (W/m2 ground)
    !---------------------------------------------------------------------

    associate ( &
                                                 ! *** Input ***
    emleaf   => pftcon%emleaf               , &  ! Leaf emissivity (-)
    tg       => mlcanopy_inst%tg            , &  ! Soil surface temperature (K)
    ncan     => mlcanopy_inst%ncan          , &  ! Number of aboveground layers
    nbot     => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop     => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    dpai     => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    irsky    => mlcanopy_inst%irsky         , &  ! Atmospheric longwave radiation (W/m2)
    fracsun  => mlcanopy_inst%fracsun       , &  ! Sunlit fraction of canopy layer
    fracsha  => mlcanopy_inst%fracsha       , &  ! Shaded fraction of canopy layer
    tleaf    => mlcanopy_inst%tleaf         , &  ! Leaf temperature (K)
    td       => mlcanopy_inst%td            , &  ! Exponential transmittance of diffuse radiation through a single leaf layer
                                                 ! *** Output ***
    irleaf   => mlcanopy_inst%irleaf        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
    irveg    => mlcanopy_inst%irveg         , &  ! Absorbed longwave radiation, vegetation (W/m2)
    ircan    => mlcanopy_inst%ircan         , &  ! Upward longwave radiation above canopy (W/m2)
    irsoi    => mlcanopy_inst%irsoi           &  ! Absorbed longwave radiation, ground (W/m2)
    )

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Zero out radiative fluxes for all layers

       irup(p,0) = 0._r8
       irdn(p,0) = 0._r8

       do ic = 1, ncan(p)
          irup(p,ic) = 0.0;
          irdn(p,ic) = 0.0;
          irleaf(p,ic) = 0.0;
       end do

       ! Ground (soil) emissivity

       emg = 0.96_r8

       !------------------------------------------------------------------
       ! Leaf scattering coefficient and terms for longwave radiation reflected
       ! and transmitted by a layer
       !------------------------------------------------------------------

       omega = 1._r8 - emleaf(patch%itype(p))

       ! Intercepted radiation is reflected

       rho = omega 
       tau = 0._r8

       ! Intercepted radiation is both reflected and transmitted

!      rho = omega * 0.5_r8
!      tau = omega * 0.5_r8

       !------------------------------------------------------------------
       ! Emitted longwave radiation is weighted average of sunlit and shaded leaves
       !------------------------------------------------------------------

       do ic = nbot(p), ntop(p)
          ir_source_sun = emleaf(patch%itype(p)) * sb * tleaf(p,ic,isun)**4
          ir_source_sha = emleaf(patch%itype(p)) * sb * tleaf(p,ic,isha)**4
          ir_source(ic) = (ir_source_sun * fracsun(p,ic) + ir_source_sha * fracsha(p,ic)) * (1._r8 - td(p,ic))
       end do

       !------------------------------------------------------------------
       ! Set up and solve tridiagonal system of equations for upward and downward fluxes
       !------------------------------------------------------------------

       ! There are two equations for each leaf layer and the soil. The first
       ! equation is the upward flux and the second equation is the downward flux.

       m = 0

       ! Soil: upward flux

       m = m + 1
       atri(m) = 0._r8
       btri(m) = 1._r8
       ctri(m) = -(1._r8 - emg)
       dtri(m) = emg * sb * tg(p)**4

       ! Soil: downward flux

       refld = (1._r8 - td(p,nbot(p))) * rho
       trand = (1._r8 - td(p,nbot(p))) * tau + td(p,nbot(p))
       aic = refld - trand * trand / refld
       bic = trand / refld

       m = m + 1
       atri(m) = -aic
       btri(m) = 1._r8
       ctri(m) = -bic
       dtri(m) = (1._r8 - bic) * ir_source(nbot(p))

       ! Leaf layers, excluding top layer

       do ic = nbot(p), ntop(p)-1

          ! Upward flux

          refld = (1._r8 - td(p,ic)) * rho
          trand = (1._r8 - td(p,ic)) * tau + td(p,ic)
          fic = refld - trand * trand / refld
          eic = trand / refld

          m = m + 1
          atri(m) = -eic
          btri(m) = 1._r8
          ctri(m) = -fic
          dtri(m) = (1._r8 - eic) * ir_source(ic)

          ! Downward flux

          refld = (1._r8 - td(p,ic+1)) * rho
          trand = (1._r8 - td(p,ic+1)) * tau + td(p,ic+1)
          aic = refld - trand * trand / refld
          bic = trand / refld

          m = m + 1
          atri(m) = -aic
          btri(m) = 1._r8
          ctri(m) = -bic
          dtri(m) = (1._r8 - bic) * ir_source(ic+1)

       end do

       ! Top canopy layer: upward flux

       ic = ntop(p)
       refld = (1._r8 - td(p,ic)) * rho
       trand = (1._r8 - td(p,ic)) * tau + td(p,ic)
       fic = refld - trand * trand / refld
       eic = trand / refld

       m = m + 1
       atri(m) = -eic
       btri(m) = 1._r8
       ctri(m) = -fic
       dtri(m) = (1._r8 - eic) * ir_source(ic)

       ! Top canopy layer: downward flux

       m = m + 1
       atri(m) = 0._r8
       btri(m) = 1._r8
       ctri(m) = 0._r8
       dtri(m) = irsky(p)

       ! Solve tridiagonal system of equations for upward and downward fluxes

       call tridiag (atri, btri, ctri, dtri, utri, m)

       ! Now copy the solution (utri) to the upward (irup) and downward (irdn)
       ! fluxes for each layer
       ! irup =  Upward longwave flux above layer
       ! irdn =  Downward longwave flux onto layer

       m = 0

       ! Soil fluxes

       m = m + 1
       irup(p,0) = utri(m)
       m = m + 1
       irdn(p,0) = utri(m)

       ! Leaf layer fluxes

       do ic = nbot(p), ntop(p)
          m = m + 1
          irup(p,ic) = utri(m)
          m = m + 1
          irdn(p,ic) = utri(m)
       end do

       !------------------------------------------------------------------
       ! Compute fluxes
       !------------------------------------------------------------------

       ! Absorbed longwave radiation for ground (soil)

       irsoi(p) = irdn(p,0) - irup(p,0)

       ! Leaf layer fluxes

       irveg(p) = 0._r8

       do ic = nbot(p), ntop(p)

          ! Absorbed longwave radiation for layer. Note special case for first
          ! leaf layer, where the upward flux from below is from the ground.
          ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
          ! canopy layers without leaves.

          if (ic == nbot(p)) then
             icm1 = 0
          else
             icm1 = ic - 1
          end if
          irabs = emleaf(patch%itype(p)) * (irdn(p,ic)+irup(p,icm1)) * (1._r8 - td(p,ic)) - 2._r8 * ir_source(ic)
          irleaf(p,ic) = irabs / dpai(p,ic)

          ! Sum longwave radiation absorbed by vegetation

          irveg(p) = irveg(p) + irabs

       end do

       ! Canopy emitted longwave radiation

       ircan(p) = irup(p,ntop(p))

       !------------------------------------------------------------------
       ! Conservation check
       !------------------------------------------------------------------

       ! Total radiation balance: absorbed = incoming - outgoing

       sumabs = irsky(p) - ircan(p)
       error = sumabs - (irveg(p) + irsoi(p))
       if (abs(error) > 1.e-03) then
          call endrun (msg='ERROR: LongwaveRadiationMod: total longwave conservation error')
       end if

    end do

    end associate
  end subroutine LongwaveRadiation

end module LongwaveRadiationMod

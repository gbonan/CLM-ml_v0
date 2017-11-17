module CanopyFluxesMultilayerType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Multilayer canopy module data structure
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevgrnd, numrad, nlevcan, nleaf
  use clm_varcon, only : ispval, spval, nan => spval
  use abortutils, only : endrun
  use decompMod , only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC DATA TYPES:

  type, public :: mlcanopy_type

    ! Vegetation input variables

    real(r8), pointer :: ztop(:)            ! Canopy height (m)
    real(r8), pointer :: lai(:)             ! Leaf area index of canopy (m2/m2)
    real(r8), pointer :: sai(:)             ! Stem area index of canopy (m2/m2)
    real(r8), pointer :: root_biomass(:)    ! Fine root biomass (g biomass / m2)
    integer , pointer :: ncan(:)            ! Number of aboveground layers
    integer , pointer :: nbot(:)            ! Index for bottom leaf layer
    integer , pointer :: ntop(:)            ! Index for top leaf layer
    real(r8), pointer :: dlai(:,:)          ! Layer leaf area index (m2/m2)
    real(r8), pointer :: dsai(:,:)          ! Layer stem area index (m2/m2)
    real(r8), pointer :: dpai(:,:)          ! Layer plant area index (m2/m2)
    real(r8), pointer :: sumpai(:,:)        ! Cumulative plant area index (m2/m2) [for nlevcan layers]
    real(r8), pointer :: zs(:,:)            ! Canopy height for scalar concentration and source (m)
    real(r8), pointer :: zw(:,:)            ! Canopy heights at layer interfaces (m)

    ! Atmospheric input variables

    real(r8), pointer :: zref(:)            ! Reference height (m)
    real(r8), pointer :: zref_old(:)        ! Reference height for previous timestep (m)
    real(r8), pointer :: tref(:)            ! Air temperature at reference height (K)
    real(r8), pointer :: uref(:)            ! Wind speed at reference height (m/s)
    real(r8), pointer :: rhref(:)           ! Relative humidity at reference height (%)
    real(r8), pointer :: pref(:)            ! Air pressure at reference height (Pa)
    real(r8), pointer :: co2ref(:)          ! Atmospheric CO2 at reference height (umol/mol)
    real(r8), pointer :: o2ref(:)           ! Atmospheric O2 at reference height (mmol/mol)
    real(r8), pointer :: solar_zen(:)       ! Solar zenith angle (radians)
    real(r8), pointer :: swskyb(:,:)        ! Atmospheric direct beam solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swskyd(:,:)        ! Atmospheric diffuse solar radiation (W/m2) [for numrad wavebands]
    real(r8), pointer :: irsky(:)           ! Atmospheric longwave radiation (W/m2)
    real(r8), pointer :: qflx_rain(:)       ! Rainfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: qflx_snow(:)       ! Snowfall (mm H2O/s = kg H2O/m2/s)
    real(r8), pointer :: tacclim(:)         ! Average air temperature for acclimation (K)

    ! Additional derived input variables

    real(r8), pointer :: eref(:)            ! Vapor pressure at reference height (Pa)
    real(r8), pointer :: qref(:)            ! Specific humidity at reference height (kg/kg)
    real(r8), pointer :: rhoair(:)          ! Air density at reference height (kg/m3)
    real(r8), pointer :: rhomol(:)          ! Molar density at reference height (mol/m3)
    real(r8), pointer :: mmair(:)           ! Molecular mass of air at reference height (kg/mol)
    real(r8), pointer :: cpair(:)           ! Specific heat of air at constant pressure, at reference height (J/mol/K)

    ! Canopy layer variables [for nlevcan layers]

    real(r8), pointer :: vcmax25(:,:)       ! Leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: jmax25(:,:)        ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: kp25(:,:)          ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)
    real(r8), pointer :: rd25(:,:)          ! Leaf respiration rate at 25C for canopy layer (umol/m2/s)
    real(r8), pointer :: wind(:,:)          ! Wind speed profile (m/s)
    real(r8), pointer :: wind_most(:,:)     ! Wind speed profile from MOST (m/s)
    real(r8), pointer :: tair(:,:)          ! Air temperature profile (K)
    real(r8), pointer :: tair_most(:,:)     ! Air temperature profile from MOST (K)
    real(r8), pointer :: eair(:,:)          ! Vapor pressure profile (Pa)
    real(r8), pointer :: cair(:,:)          ! Atmospheric CO2 profile (umol/mol)
    real(r8), pointer :: tveg(:,:,:)        ! Vegetation temperature profile (K)
    real(r8), pointer :: tair_old(:,:)      ! Air temperature profile for previous timestep (K)
    real(r8), pointer :: eair_old(:,:)      ! Vapor pressure profile for previous timestep (Pa)
    real(r8), pointer :: cair_old(:,:)      ! Atmospheric CO2 profile for previous timestep (umol/mol)
    real(r8), pointer :: tveg_old(:,:,:)    ! Vegetation temperature profile for previous timestep (K)
    real(r8), pointer :: fracsun(:,:)       ! Sunlit fraction of canopy layer
    real(r8), pointer :: fracsha(:,:)       ! Shaded fraction of canopy layer
    real(r8), pointer :: irleaf(:,:)        ! Leaf absorbed longwave radiation for canopy layer(W/m2 leaf)
    real(r8), pointer :: lwp(:,:)           ! Leaf water potential of canopy layer (MPa)
    real(r8), pointer :: lsc(:,:)           ! Leaf-specific conductance of canopy layer (mmol H2O/m2 leaf/s/MPa)
    real(r8), pointer :: h2ocan(:,:)        ! Canopy layer intercepted water (kg H2O/m2)
    real(r8), pointer :: fwet(:,:)          ! Fraction of plant area index that is wet
    real(r8), pointer :: fdry(:,:)          ! Fraction of plant area index that is green and dry
    real(r8), pointer :: shair(:,:)         ! Canopy air sensible heat flux (W/m2)
    real(r8), pointer :: etair(:,:)         ! Canopy air water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: stair(:,:)         ! Canopy air storage heat flux (W/m2)
    real(r8), pointer :: sw_prof(:,:,:)     ! Canopy layer absorbed solar radiation (W/m2)
    real(r8), pointer :: ir_prof(:,:)       ! Canopy layer absorbed longwave radiation (W/m2)
    real(r8), pointer :: rn_prof(:,:)       ! Canopy layer net radiation (W/m2)
    real(r8), pointer :: st_prof(:,:)       ! Canopy layer storage heat flux (W/m2)
    real(r8), pointer :: sh_prof(:,:)       ! Canopy layer sensible heat flux (W/m2)
    real(r8), pointer :: lh_prof(:,:)       ! Canopy layer latent heat flux (W/m2)
    real(r8), pointer :: et_prof(:,:)       ! Canopy layer water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: fc_prof(:,:)       ! Canopy layer CO2 flux (umol CO2/m2/s)
    real(r8), pointer :: ga_prof(:,:)       ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)

    ! Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]

    real(r8), pointer :: tleaf(:,:,:)       ! Leaf temperature (K)
    real(r8), pointer :: tleaf_old(:,:,:)   ! Leaf temperature for previous timestep (K)
    real(r8), pointer :: rnleaf(:,:,:)      ! Leaf net radiation (W/m2 leaf)
    real(r8), pointer :: stleaf(:,:,:)      ! Leaf storage heat flux (W/m2 leaf)
    real(r8), pointer :: shleaf(:,:,:)      ! Leaf sensible heat flux (W/m2 leaf)
    real(r8), pointer :: lhleaf(:,:,:)      ! Leaf latent heat flux (W/m2 leaf)
    real(r8), pointer :: swleaf(:,:,:,:)    ! Leaf absorbed solar radiation (W/m2 leaf) [for numrad wavebands]
    real(r8), pointer :: trleaf(:,:,:)      ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: evleaf(:,:,:)      ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    real(r8), pointer :: psil(:,:,:)        ! Leaf water potential (MPa)

    real(r8), pointer :: gbh(:,:,:)         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    real(r8), pointer :: gbv(:,:,:)         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    real(r8), pointer :: gbc(:,:,:)         ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)

    real(r8), pointer :: apar(:,:,:)        ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    real(r8), pointer :: ac(:,:,:)          ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: aj(:,:,:)          ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ap(:,:,:)          ! Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: ag(:,:,:)          ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: an(:,:,:)          ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    real(r8), pointer :: rd(:,:,:)          ! Leaf respiration rate (umol CO2/m2 leaf/s)
    real(r8), pointer :: ci(:,:,:)          ! Leaf intercellular CO2 (umol/mol)
    real(r8), pointer :: cs(:,:,:)          ! Leaf surface CO2 (umol/mol)
    real(r8), pointer :: hs(:,:,:)          ! Leaf fractional humidity at leaf surface (-)
    real(r8), pointer :: vpd(:,:,:)         ! Leaf vapor pressure deficit (Pa)
    real(r8), pointer :: gs(:,:,:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    real(r8), pointer :: alphapsn(:,:,:)    ! Leaf 13C fractionation factor for photosynthesis (-)

    real(r8), pointer :: cpleaf(:,:)        ! Leaf heat capacity (J/m2 leaf/K)

    ! Canopy variables (fluxes are per m2 ground area)

    real(r8), pointer :: swveg(:,:)         ! Absorbed solar radiation, vegetation (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsun(:,:)      ! Absorbed solar radiation, sunlit canopy (W/m2) [for numrad wavebands]
    real(r8), pointer :: swvegsha(:,:)      ! Absorbed solar radiation, shaded canopy (W/m2) [for numrad wavebands]

    real(r8), pointer :: irveg(:)           ! Absorbed longwave radiation, vegetation (W/m2)
    real(r8), pointer :: irvegsun(:)        ! Absorbed longwave radiation, sunlit canopy (W/m2)
    real(r8), pointer :: irvegsha(:)        ! Absorbed longwave radiation, shaded canopy (W/m2)

    real(r8), pointer :: shveg(:)           ! Sensible heat flux, vegetation (W/m2)
    real(r8), pointer :: shvegsun(:)        ! Sensible heat flux, sunlit canopy (W/m2)
    real(r8), pointer :: shvegsha(:)        ! Sensible heat flux, shaded canopy (W/m2)

    real(r8), pointer :: lhveg(:)           ! Latent heat flux, vegetation (W/m2)
    real(r8), pointer :: lhvegsun(:)        ! Latent heat flux, sunlit canopy (W/m2)
    real(r8), pointer :: lhvegsha(:)        ! Latent heat flux, shaded canopy (W/m2)

    real(r8), pointer :: etveg(:)           ! Water vapor flux, vegetation (mol H2O/m2/s)
    real(r8), pointer :: etvegsun(:)        ! Water vapor flux, sunlit canopy (mol H2O/m2/s)
    real(r8), pointer :: etvegsha(:)        ! Water vapor flux, shaded canopy (mol H2O/m2/s)

    real(r8), pointer :: gppveg(:)          ! Gross primary production (umol CO2/m2/s)
    real(r8), pointer :: gppvegsun(:)       ! Gross primary production, sunlit canopy (umol CO2/m2/s)
    real(r8), pointer :: gppvegsha(:)       ! Gross primary production, shaded canopy (umol CO2/m2/s)

    real(r8), pointer :: albcan(:,:)        ! Albedo above canopy [for numrad wavebands]
    real(r8), pointer :: ircan(:)           ! Upward longwave radiation above canopy (W/m2)
    real(r8), pointer :: rnet(:)            ! Net radiation (W/m2)
    real(r8), pointer :: stflx(:)           ! Canopy storage heat flux (W/m2)
    real(r8), pointer :: shflx(:)           ! Sensible heat flux (W/m2)
    real(r8), pointer :: lhflx(:)           ! Latent heat flux (W/m2)
    real(r8), pointer :: etflx(:)           ! Water vapor flux (mol H2O/m2/s)
    real(r8), pointer :: fracminlwp(:)      ! Fraction of canopy with lwp < minlwp

    real(r8), pointer :: ustar(:)           ! Friction velocity (m/s)
    real(r8), pointer :: uforc(:)           ! Wind speed at reference height including stability effect (m/s)
    real(r8), pointer :: uaf(:)             ! Wind speed at canopy top (m/s)
    real(r8), pointer :: taf(:)             ! Air temperature at canopy top (K)
    real(r8), pointer :: qaf(:)             ! Specific humidity at canopy top (kg/kg)
    real(r8), pointer :: eaf(:)             ! Vapor pressure at canopy top (Pa)
    real(r8), pointer :: obu(:)             ! Obukhov length (m)
    real(r8), pointer :: obu_gah(:)         ! Obukhov length used for gah (m)
    real(r8), pointer :: obuold(:)          ! Obukhov length from previous iteration
    integer,  pointer :: nmozsgn(:)         ! Number of times stability changes sign during iteration

    real(r8), pointer :: z0mg(:)            ! Roughness length of ground (m)
    real(r8), pointer :: thref(:)           ! Atmospheric potential temperature (K)
    real(r8), pointer :: thvref(:)          ! Atmospheric virtual potential temperature (K)
    real(r8), pointer :: gah(:)             ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    real(r8), pointer :: PrSc(:)            ! Prandtl (Schmidt) number at canopy top
    real(r8), pointer :: Lc(:)              ! Canopy density length scale (m)
    real(r8), pointer :: zdisp(:)           ! Displacement height (m)
    real(r8), pointer :: tstar(:)           ! Temperature scale (K)
    real(r8), pointer :: qstar(:)           ! Water vapor scale (kg/kg)

    real(r8), pointer :: td(:,:)            ! Exponential transmittance of diffuse radiation through a single leaf layer

    ! Soil energy balance

    real(r8), pointer :: rnsoi(:)           ! Net radiation, ground (W/m2)
    real(r8), pointer :: shsoi(:)           ! Sensible heat flux, ground (W/m2)
    real(r8), pointer :: lhsoi(:)           ! Latent heat flux, ground (W/m2)
    real(r8), pointer :: gsoi(:)            ! Soil heat flux (W/m2)
    real(r8), pointer :: swsoi(:,:)         ! Absorbed solar radiation, ground (W/m2) [for numrad wavebands]
    real(r8), pointer :: irsoi(:)           ! Absorbed longwave radiation, ground (W/m2)
    real(r8), pointer :: etsoi(:)           ! Water vapor flux, ground (mol H2O/m2/s)
    real(r8), pointer :: tg(:)              ! Soil surface temperature (K)

    ! Soil moisture variables

    real(r8), pointer :: btran(:)           ! Ball-Berry soil wetness factor (-)
    real(r8), pointer :: psis(:)            ! Weighted soil water potential (MPa)
    real(r8), pointer :: rsoil(:)           ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    real(r8), pointer :: soil_et_loss(:,:)  ! Fraction of total transpiration from each soil layer (-)
    real(r8), pointer :: eg(:)              ! Soil surface vapor pressure (Pa)
    real(r8), pointer :: rhg(:)             ! Relative humidity of airspace at soil surface (fraction)

    ! Water flux variables

    real(r8), pointer :: qflx_prec_intr(:)  ! Intercepted precipitation (kg H2O/m2/s)

  contains

    procedure, public  :: Init              ! CLM initialization of data type
    procedure, private :: InitAllocate      ! CLM initialization: allocate module data structure
    procedure, private :: InitHistory       ! CLM initialization: setup history file variables
    procedure, private :: InitCold          ! CLM initialization: cold-start initialization
    procedure, public  :: Restart           ! CLM restart file

  end type mlcanopy_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start
    !
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate module data structure
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp      ! Beginning patch index for CLM g/l/c/p hierarchy
    integer :: endp      ! Ending patch index for CLM g/l/c/p hierarchy
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate (this%ztop           (begp:endp))                            ; this%ztop           (:)       = nan
    allocate (this%lai            (begp:endp))                            ; this%lai            (:)       = nan
    allocate (this%sai            (begp:endp))                            ; this%sai            (:)       = nan
    allocate (this%root_biomass   (begp:endp))                            ; this%root_biomass   (:)       = nan
    allocate (this%ncan           (begp:endp))                            ; this%ncan           (:)       = ispval
    allocate (this%nbot           (begp:endp))                            ; this%nbot           (:)       = ispval
    allocate (this%ntop           (begp:endp))                            ; this%ntop           (:)       = ispval
    allocate (this%dlai           (begp:endp,1:nlevcan))                  ; this%dlai           (:,:)     = nan
    allocate (this%dsai           (begp:endp,1:nlevcan))                  ; this%dsai           (:,:)     = nan
    allocate (this%dpai           (begp:endp,1:nlevcan))                  ; this%dpai           (:,:)     = nan
    allocate (this%sumpai         (begp:endp,1:nlevcan))                  ; this%sumpai         (:,:)     = nan
    allocate (this%zs             (begp:endp,0:nlevcan))                  ; this%zs             (:,:)     = nan
    allocate (this%zw             (begp:endp,0:nlevcan))                  ; this%zw             (:,:)     = nan
    allocate (this%zref           (begp:endp))                            ; this%zref           (:)       = nan
    allocate (this%zref_old       (begp:endp))                            ; this%zref_old       (:)       = nan
    allocate (this%tref           (begp:endp))                            ; this%tref           (:)       = nan
    allocate (this%uref           (begp:endp))                            ; this%uref           (:)       = nan
    allocate (this%rhref          (begp:endp))                            ; this%rhref          (:)       = nan
    allocate (this%pref           (begp:endp))                            ; this%pref           (:)       = nan
    allocate (this%co2ref         (begp:endp))                            ; this%co2ref         (:)       = nan
    allocate (this%o2ref          (begp:endp))                            ; this%o2ref          (:)       = nan
    allocate (this%solar_zen      (begp:endp))                            ; this%solar_zen      (:)       = nan
    allocate (this%swskyb         (begp:endp,1:numrad))                   ; this%swskyb         (:,:)     = nan
    allocate (this%swskyd         (begp:endp,1:numrad))                   ; this%swskyd         (:,:)     = nan
    allocate (this%irsky          (begp:endp))                            ; this%irsky          (:)       = nan
    allocate (this%qflx_rain      (begp:endp))                            ; this%qflx_rain      (:)       = nan
    allocate (this%qflx_snow      (begp:endp))                            ; this%qflx_snow      (:)       = nan
    allocate (this%tacclim        (begp:endp))                            ; this%tacclim        (:)       = nan
    allocate (this%eref           (begp:endp))                            ; this%eref           (:)       = nan
    allocate (this%qref           (begp:endp))                            ; this%qref           (:)       = nan
    allocate (this%rhoair         (begp:endp))                            ; this%rhoair         (:)       = nan
    allocate (this%rhomol         (begp:endp))                            ; this%rhomol         (:)       = nan
    allocate (this%mmair          (begp:endp))                            ; this%mmair          (:)       = nan
    allocate (this%cpair          (begp:endp))                            ; this%cpair          (:)       = nan
    allocate (this%vcmax25        (begp:endp,1:nlevcan))                  ; this%vcmax25        (:,:)     = nan
    allocate (this%jmax25         (begp:endp,1:nlevcan))                  ; this%jmax25         (:,:)     = nan
    allocate (this%kp25           (begp:endp,1:nlevcan))                  ; this%kp25           (:,:)     = nan
    allocate (this%rd25           (begp:endp,1:nlevcan))                  ; this%rd25           (:,:)     = nan
    allocate (this%wind           (begp:endp,0:nlevcan))                  ; this%wind           (:,:)     = nan
    allocate (this%wind_most      (begp:endp,0:nlevcan))                  ; this%wind_most      (:,:)     = nan
    allocate (this%tair           (begp:endp,0:nlevcan))                  ; this%tair           (:,:)     = nan
    allocate (this%tair_most      (begp:endp,0:nlevcan))                  ; this%tair_most      (:,:)     = nan
    allocate (this%eair           (begp:endp,0:nlevcan))                  ; this%eair           (:,:)     = nan
    allocate (this%cair           (begp:endp,0:nlevcan))                  ; this%cair           (:,:)     = nan
    allocate (this%tveg           (begp:endp,0:nlevcan,1:nleaf))          ; this%tveg           (:,:,:)   = nan
    allocate (this%tair_old       (begp:endp,0:nlevcan))                  ; this%tair_old       (:,:)     = nan
    allocate (this%eair_old       (begp:endp,0:nlevcan))                  ; this%eair_old       (:,:)     = nan
    allocate (this%cair_old       (begp:endp,0:nlevcan))                  ; this%cair_old       (:,:)     = nan
    allocate (this%tveg_old       (begp:endp,0:nlevcan,1:nleaf))          ; this%tveg_old       (:,:,:)   = nan
    allocate (this%fracsun        (begp:endp,1:nlevcan))                  ; this%fracsun        (:,:)     = nan
    allocate (this%fracsha        (begp:endp,1:nlevcan))                  ; this%fracsha        (:,:)     = nan
    allocate (this%irleaf         (begp:endp,1:nlevcan))                  ; this%irleaf         (:,:)     = nan
    allocate (this%lwp            (begp:endp,1:nlevcan))                  ; this%lwp            (:,:)     = nan
    allocate (this%lsc            (begp:endp,1:nlevcan))                  ; this%lsc            (:,:)     = nan
    allocate (this%h2ocan         (begp:endp,1:nlevcan))                  ; this%h2ocan         (:,:)     = nan
    allocate (this%fwet           (begp:endp,1:nlevcan))                  ; this%fwet           (:,:)     = nan
    allocate (this%fdry           (begp:endp,1:nlevcan))                  ; this%fdry           (:,:)     = nan
    allocate (this%shair          (begp:endp,1:nlevcan))                  ; this%shair          (:,:)     = nan
    allocate (this%etair          (begp:endp,1:nlevcan))                  ; this%etair          (:,:)     = nan
    allocate (this%stair          (begp:endp,1:nlevcan))                  ; this%stair          (:,:)     = nan
    allocate (this%sw_prof        (begp:endp,0:nlevcan,1:numrad))         ; this%sw_prof        (:,:,:)   = nan
    allocate (this%ir_prof        (begp:endp,0:nlevcan))                  ; this%ir_prof        (:,:)     = nan
    allocate (this%rn_prof        (begp:endp,0:nlevcan))                  ; this%rn_prof        (:,:)     = nan
    allocate (this%st_prof        (begp:endp,0:nlevcan))                  ; this%st_prof        (:,:)     = nan
    allocate (this%sh_prof        (begp:endp,0:nlevcan))                  ; this%sh_prof        (:,:)     = nan
    allocate (this%lh_prof        (begp:endp,0:nlevcan))                  ; this%lh_prof        (:,:)     = nan
    allocate (this%et_prof        (begp:endp,0:nlevcan))                  ; this%et_prof        (:,:)     = nan
    allocate (this%fc_prof        (begp:endp,0:nlevcan))                  ; this%fc_prof        (:,:)     = nan
    allocate (this%ga_prof        (begp:endp,0:nlevcan))                  ; this%ga_prof        (:,:)     = nan
    allocate (this%tleaf          (begp:endp,1:nlevcan,1:nleaf))          ; this%tleaf          (:,:,:)   = nan
    allocate (this%tleaf_old      (begp:endp,1:nlevcan,1:nleaf))          ; this%tleaf_old      (:,:,:)   = nan
    allocate (this%rnleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%rnleaf         (:,:,:)   = nan
    allocate (this%stleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%stleaf         (:,:,:)   = nan
    allocate (this%shleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%shleaf         (:,:,:)   = nan
    allocate (this%lhleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%lhleaf         (:,:,:)   = nan
    allocate (this%swleaf         (begp:endp,1:nlevcan,1:nleaf,1:numrad)) ; this%swleaf         (:,:,:,:) = nan
    allocate (this%trleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%trleaf         (:,:,:)   = nan
    allocate (this%evleaf         (begp:endp,1:nlevcan,1:nleaf))          ; this%evleaf         (:,:,:)   = nan
    allocate (this%psil           (begp:endp,1:nlevcan,1:nleaf))          ; this%psil           (:,:,:)   = nan
    allocate (this%gbh            (begp:endp,1:nlevcan,1:nleaf))          ; this%gbh            (:,:,:)   = nan
    allocate (this%gbv            (begp:endp,1:nlevcan,1:nleaf))          ; this%gbv            (:,:,:)   = nan
    allocate (this%gbc            (begp:endp,1:nlevcan,1:nleaf))          ; this%gbc            (:,:,:)   = nan
    allocate (this%apar           (begp:endp,1:nlevcan,1:nleaf))          ; this%apar           (:,:,:)   = nan
    allocate (this%ac             (begp:endp,1:nlevcan,1:nleaf))          ; this%ac             (:,:,:)   = nan
    allocate (this%aj             (begp:endp,1:nlevcan,1:nleaf))          ; this%aj             (:,:,:)   = nan
    allocate (this%ap             (begp:endp,1:nlevcan,1:nleaf))          ; this%ap             (:,:,:)   = nan
    allocate (this%ag             (begp:endp,1:nlevcan,1:nleaf))          ; this%ag             (:,:,:)   = nan
    allocate (this%an             (begp:endp,1:nlevcan,1:nleaf))          ; this%an             (:,:,:)   = nan
    allocate (this%rd             (begp:endp,1:nlevcan,1:nleaf))          ; this%rd             (:,:,:)   = nan
    allocate (this%ci             (begp:endp,1:nlevcan,1:nleaf))          ; this%ci             (:,:,:)   = nan
    allocate (this%cs             (begp:endp,1:nlevcan,1:nleaf))          ; this%cs             (:,:,:)   = nan
    allocate (this%hs             (begp:endp,1:nlevcan,1:nleaf))          ; this%hs             (:,:,:)   = nan
    allocate (this%vpd            (begp:endp,1:nlevcan,1:nleaf))          ; this%vpd            (:,:,:)   = nan
    allocate (this%gs             (begp:endp,1:nlevcan,1:nleaf))          ; this%gs             (:,:,:)   = nan
    allocate (this%alphapsn       (begp:endp,1:nlevcan,1:nleaf))          ; this%alphapsn       (:,:,:)   = nan
    allocate (this%cpleaf         (begp:endp,1:nlevcan))                  ; this%cpleaf         (:,:)     = nan
    allocate (this%swveg          (begp:endp,1:numrad))                   ; this%swveg          (:,:)     = nan
    allocate (this%swvegsun       (begp:endp,1:numrad))                   ; this%swvegsun       (:,:)     = nan
    allocate (this%swvegsha       (begp:endp,1:numrad))                   ; this%swvegsha       (:,:)     = nan
    allocate (this%irveg          (begp:endp))                            ; this%irveg          (:)       = nan
    allocate (this%irvegsun       (begp:endp))                            ; this%irvegsun       (:)       = nan
    allocate (this%irvegsha       (begp:endp))                            ; this%irvegsha       (:)       = nan
    allocate (this%shveg          (begp:endp))                            ; this%shveg          (:)       = nan
    allocate (this%shvegsun       (begp:endp))                            ; this%shvegsun       (:)       = nan
    allocate (this%shvegsha       (begp:endp))                            ; this%shvegsha       (:)       = nan
    allocate (this%lhveg          (begp:endp))                            ; this%lhveg          (:)       = nan
    allocate (this%lhvegsun       (begp:endp))                            ; this%lhvegsun       (:)       = nan
    allocate (this%lhvegsha       (begp:endp))                            ; this%lhvegsha       (:)       = nan
    allocate (this%etveg          (begp:endp))                            ; this%etveg          (:)       = nan
    allocate (this%etvegsun       (begp:endp))                            ; this%etvegsun       (:)       = nan
    allocate (this%etvegsha       (begp:endp))                            ; this%etvegsha       (:)       = nan
    allocate (this%gppveg         (begp:endp))                            ; this%gppveg         (:)       = nan
    allocate (this%gppvegsun      (begp:endp))                            ; this%gppvegsun      (:)       = nan
    allocate (this%gppvegsha      (begp:endp))                            ; this%gppvegsha      (:)       = nan
    allocate (this%albcan         (begp:endp,1:numrad))                   ; this%albcan         (:,:)     = nan
    allocate (this%ircan          (begp:endp))                            ; this%ircan          (:)       = nan
    allocate (this%rnet           (begp:endp))                            ; this%rnet           (:)       = nan
    allocate (this%stflx          (begp:endp))                            ; this%stflx          (:)       = nan
    allocate (this%shflx          (begp:endp))                            ; this%shflx          (:)       = nan
    allocate (this%lhflx          (begp:endp))                            ; this%lhflx          (:)       = nan
    allocate (this%etflx          (begp:endp))                            ; this%etflx          (:)       = nan
    allocate (this%fracminlwp     (begp:endp))                            ; this%fracminlwp     (:)       = nan
    allocate (this%ustar          (begp:endp))                            ; this%ustar          (:)       = nan
    allocate (this%uforc          (begp:endp))                            ; this%uforc          (:)       = nan
    allocate (this%uaf            (begp:endp))                            ; this%uaf            (:)       = nan
    allocate (this%taf            (begp:endp))                            ; this%taf            (:)       = nan
    allocate (this%qaf            (begp:endp))                            ; this%qaf            (:)       = nan
    allocate (this%eaf            (begp:endp))                            ; this%eaf            (:)       = nan
    allocate (this%obu            (begp:endp))                            ; this%obu            (:)       = nan
    allocate (this%obu_gah        (begp:endp))                            ; this%obu_gah        (:)       = nan
    allocate (this%obuold         (begp:endp))                            ; this%obuold         (:)       = nan
    allocate (this%nmozsgn        (begp:endp))                            ; this%nmozsgn        (:)       = ispval
    allocate (this%z0mg           (begp:endp))                            ; this%z0mg           (:)       = nan
    allocate (this%thref          (begp:endp))                            ; this%thref          (:)       = nan
    allocate (this%thvref         (begp:endp))                            ; this%thvref         (:)       = nan
    allocate (this%gah            (begp:endp))                            ; this%gah            (:)       = nan
    allocate (this%PrSc           (begp:endp))                            ; this%PrSc           (:)       = nan
    allocate (this%Lc             (begp:endp))                            ; this%Lc             (:)       = nan
    allocate (this%zdisp          (begp:endp))                            ; this%zdisp          (:)       = nan
    allocate (this%tstar          (begp:endp))                            ; this%tstar          (:)       = nan
    allocate (this%qstar          (begp:endp))                            ; this%qstar          (:)       = nan
    allocate (this%td             (begp:endp,1:nlevcan))                  ; this%td             (:,:)     = nan
    allocate (this%rnsoi          (begp:endp))                            ; this%rnsoi          (:)       = nan
    allocate (this%shsoi          (begp:endp))                            ; this%shsoi          (:)       = nan
    allocate (this%lhsoi          (begp:endp))                            ; this%lhsoi          (:)       = nan
    allocate (this%gsoi           (begp:endp))                            ; this%gsoi           (:)       = nan
    allocate (this%swsoi          (begp:endp,1:numrad))                   ; this%swsoi          (:,:)     = nan
    allocate (this%irsoi          (begp:endp))                            ; this%irsoi          (:)       = nan
    allocate (this%etsoi          (begp:endp))                            ; this%etsoi          (:)       = nan
    allocate (this%tg             (begp:endp))                            ; this%tg             (:)       = nan
    allocate (this%btran          (begp:endp))                            ; this%btran          (:)       = nan
    allocate (this%psis           (begp:endp))                            ; this%psis           (:)       = nan
    allocate (this%rsoil          (begp:endp))                            ; this%rsoil          (:)       = nan
    allocate (this%soil_et_loss   (begp:endp,1:nlevgrnd))                 ; this%soil_et_loss   (:,:)     = nan
    allocate (this%eg             (begp:endp))                            ; this%eg             (:)       = nan
    allocate (this%rhg            (begp:endp))                            ; this%rhg            (:)       = nan
    allocate (this%qflx_prec_intr (begp:endp))                            ; this%qflx_prec_intr (:)       = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory (this, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files
    !
    ! !USES:
    use histFileMod, only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%gppveg(begp:endp) = spval
    call hist_addfld1d (fname='GPPVEG', units='umol/m2s', &
         avgflag='A', long_name='Gross primary production', &
         ptr_patch=this%gppveg, set_lake=0._r8, set_urb=0._r8, default='inactive')

    this%lwp(begp:endp,1:nlevcan) = spval
    call hist_addfld2d (fname='LWP', units='MPa', type2d='nlevcan', &
         avgflag='A', long_name='Leaf water potential of canopy layer', &
         ptr_patch=this%lwp, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold (this, bounds)
    !
    ! !DESCRIPTION:
    ! Cold-start initialization for multilayer canopy
    !
    ! !USES:
    use clm_varpar, only : nlevcan
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: iv               ! Canopy leaf layer index
    !---------------------------------------------------------------------

    ! Initialize leaf water potential and intercepted water

    do p = bounds%begp, bounds%endp
       do iv = 1, nlevcan
          this%lwp(p,iv) = -0.1_r8
          this%h2ocan(p,iv) = 0._r8
       end do
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod, only : restartvar
    !
    ! !ARGUMENTS:
    class(mlcanopy_type) :: this
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !---------------------------------------------------------------------

    ! Example for 2-d patch variable

    call restartvar(ncid=ncid, flag=flag, varname='lwp', xtype=ncd_double,  &
       dim1name='pft', dim2name='levcan', switchdim=.true., &
       long_name='leaf water potential of canopy layer', units='MPa', &
       interpinic_flag='interp', readvar=readvar, data=this%lwp)

  end subroutine Restart

end module CanopyFluxesMultilayerType

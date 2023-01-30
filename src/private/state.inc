integer :: eos_handle, net_handle
integer :: species, num_reactions

type (Net_Info) :: netinfo_target
type (Net_Info), pointer :: netinfo =>null()
type (Net_General_Info), pointer :: g =>null()

real(dp), pointer :: rate_factors(:) =>null() ! (num_reactions)
integer, pointer :: net_reaction_ptr(:) =>null() 

integer, pointer :: reaction_id(:)
integer, dimension(:), pointer :: net_iso =>null(), chem_id =>null()
         
real(dp), pointer :: ending_x(:) =>null() ! (species)
integer :: nfcn    ! number of function evaluations
integer :: njac    ! number of jacobian evaluations
integer :: nstep   ! number of computed steps
integer :: naccpt  ! number of accepted steps
integer :: nrejct  ! number of rejected steps

real(dp), dimension(:), pointer :: &
   rate_raw =>null(), rate_raw_dT =>null(), rate_raw_dRho =>null(), &
   rate_screened =>null(), rate_screened_dT =>null(), rate_screened_dRho =>null()

integer :: neut_id, prot_id

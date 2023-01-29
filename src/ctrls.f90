module ctrls
   use const_def
   use net_def 
   use net_lib

   implicit none


   type rate_t
         integer :: eos_handle, net_handle
         integer :: species, num_reactions

         type (Net_Info) :: netinfo_target
         type (Net_Info), pointer :: netinfo =>null()
         type (Net_General_Info), pointer :: g =>null()

         real(dp), pointer :: rate_factors(:) =>null() ! (num_reactions)
         integer, pointer :: net_reaction_ptr(:) =>null() 

         integer, pointer :: reaction_id(:)
         integer, dimension(:), pointer :: net_iso =>null(), chem_id =>null()

         real(dp), pointer :: burn_work_array(:) =>null(), net_work_array(:) =>null()
                  
         real(dp), pointer :: ending_x(:) =>null() ! (species)
         integer :: nfcn    ! number of function evaluations
         integer :: njac    ! number of jacobian evaluations
         integer :: nstep   ! number of computed steps
         integer :: naccpt  ! number of accepted steps
         integer :: nrejct  ! number of rejected steps

         real(dp), dimension(:), pointer :: &
            rate_raw =>null(), rate_raw_dT =>null(), rate_raw_dRho =>null(), &
            rate_screened =>null(), rate_screened_dT =>null(), rate_screened_dRho =>null()

         integer :: burn_lwork, net_lwork
         integer :: screening_opt

         integer :: neut_id, prot_id

   end type rate_t

   type bbq_t
      character(len=strlen) :: net_name,screening_mode,iso_list_filename
      integer ::  max_steps
      real(dp) :: weak_rate_factor, eps,odescal,stptry
      character(len=strlen) :: inlist
      logical :: use_input_file,use_random_sampling,use_profile
      logical :: just_write_isos,write_iso_list

      type(rate_t) :: state
   end type bbq_t


   type profile_t
      character(len=strlen) :: input_filename,output_filename,input_composition_filename
      logical :: reflective_boundaries,write_comp_every_loop
      integer :: num_loops
   end type profile_t

   type random_t
      real(dp) :: log_time_min, log_time_max
      real(dp) :: log_temp_min, log_temp_max
      real(dp) :: log_rho_min, log_rho_max
      real(dp) :: log_xa_min, log_xa_max
      real(dp) :: neut_prot_limit_frac
      character(len=strlen) :: output_starting_filename,output_ending_filename
      integer :: num_samples, seed
   end type random_t

   type sample_t
      character(len=strlen) :: input_filename,output_filename
      logical :: uniform_composition
   end type sample_t
   

   contains


   subroutine load_bbq_inputs(inlist, options, ierr)
      character(len=*), intent(in) :: inlist
      type(bbq_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit


      character(len=strlen) :: net_name, screening_mode,iso_list_filename
      integer ::  max_steps
      real(dp) :: weak_rate_factor, eps,odescal,stptry
      logical :: use_input_file,use_random_sampling,use_profile
      logical :: just_write_isos,write_iso_list
      integer  :: flush_freq
   
      namelist /bbq/ net_name, screening_mode, weak_rate_factor, &
                     eps, odescal, stptry, max_steps, &
                     use_input_file, use_random_sampling, use_profile,&
                     iso_list_filename, just_write_isos, write_iso_list

      include '../defaults/bbq.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=bbq)
      close(unit)

      options% inlist = inlist
      options% net_name = net_name
      options% screening_mode = screening_mode
      options% iso_list_filename = iso_list_filename
      options% max_steps = max_steps
      options% weak_rate_factor = weak_rate_factor
      options% eps = eps
      options% odescal = odescal
      options% stptry = stptry
      options% use_input_file = use_input_file
      options% use_random_sampling = use_random_sampling
      options% use_profile = use_profile
      options% just_write_isos = just_write_isos
      options% write_iso_list = write_iso_list

   end subroutine load_bbq_inputs


   subroutine load_profile_inputs(inlist, options, ierr)
      character(len=*), intent(in) :: inlist
      type(profile_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit

      character(len=strlen) :: input_filename,output_filename,input_composition_filename
      logical :: reflective_boundaries,write_comp_every_loop
      integer :: num_loops

      namelist /profile/ input_filename,output_filename, reflective_boundaries, num_loops,&
                        input_composition_filename,write_comp_every_loop


      include '../defaults/profile.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=profile)
      close(unit)


      options% input_filename = input_filename
      options% output_filename = output_filename
      options% reflective_boundaries = reflective_boundaries
      options% num_loops = num_loops
      options% input_composition_filename = input_composition_filename
      options% write_comp_every_loop = write_comp_every_loop


   end subroutine load_profile_inputs

   subroutine load_random_inputs(inlist, options, ierr)
      character(len=*), intent(in) :: inlist
      type(random_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit

      real(dp) :: log_time_min, log_time_max
      real(dp) :: log_temp_min, log_temp_max
      real(dp) :: log_rho_min, log_rho_max
      real(dp) :: log_xa_min, log_xa_max
      real(dp) :: neut_prot_limit_frac
      character(len=strlen) :: output_starting_filename,output_ending_filename
      integer :: num_samples, seed
   

      namelist /random/ log_time_min, log_time_max, log_temp_min, log_temp_max, &
                        log_rho_min,  log_rho_max, log_xa_min, log_xa_max, &
                        neut_prot_limit_frac, &
                        output_starting_filename, output_ending_filename,&
                        num_samples, seed


      include '../defaults/random.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=random)
      close(unit)


      options% log_time_min = log_time_min
      options% log_time_max = log_time_max
      options% log_temp_min = log_temp_min
      options% log_temp_max = log_temp_max
      options% log_rho_min = log_rho_min
      options% log_rho_max = log_rho_max
      options% log_xa_min = log_xa_min
      options% log_xa_max = log_xa_max
      options% neut_prot_limit_frac = neut_prot_limit_frac
      options% output_starting_filename = output_starting_filename
      options% output_ending_filename = output_ending_filename
      options% num_samples = num_samples
      options% seed = seed

   end subroutine load_random_inputs


   subroutine load_sampling_inputs(inlist, options, ierr)
      character(len=*), intent(in) :: inlist
      type(sample_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit

      character(len=strlen) :: input_filename,output_filename
      logical :: uniform_composition
   
      namelist /sampling/ input_filename, output_filename, uniform_composition


      include '../defaults/sampling.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=sampling)
      close(unit)


      options% input_filename = input_filename
      options% output_filename = output_filename
      options% uniform_composition = uniform_composition


   end subroutine load_sampling_inputs


end module ctrls

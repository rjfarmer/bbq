module ctrls
   use const_def
   use net_def 
   use net_lib
   use rates_def

   implicit none

   integer, parameter :: max_num_special_rate_factors = 20


   type state_t
      include 'private/state.inc'
   end type state_t

   type nuclear_t
      include 'private/nuclear.inc'

      integer :: screening_opt
   end type nuclear_t

   type bbq_t
      include 'private/bbq.inc'

      type(state_t) :: state
      type(nuclear_t) :: nuclear
   end type bbq_t


   type profile_t
      include 'private/profile.inc'
   end type profile_t

   type random_t
      include 'private/random.inc'
   end type random_t

   type sample_t
      include 'private/sample.inc'
   end type sample_t


   contains


   subroutine load_bbq_inputs(inlist_in, options, ierr)
      character(len=*), intent(in) :: inlist_in
      type(bbq_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit

      include 'private/bbq.inc'
   
      namelist /bbq/ net_name, &
                     eps, odescal, stptry, max_steps, &
                     use_input_file, use_random_sampling, use_profile,&
                     iso_list_filename, just_write_isos, write_iso_list

      include '../defaults/bbq.defaults'

      open(newunit=unit,file=inlist_in,status='old',action='read')
      read(unit,nml=bbq)
      close(unit)

      options% inlist = inlist_in
      options% net_name = net_name
      options% iso_list_filename = iso_list_filename
      options% max_steps = max_steps
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

      include 'private/profile.inc'

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

      include 'private/random.inc'
   

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

      include 'private/sample.inc'
   
      namelist /sampling/ input_filename, output_filename, uniform_composition


      include '../defaults/sampling.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=sampling)
      close(unit)


      options% input_filename = input_filename
      options% output_filename = output_filename
      options% uniform_composition = uniform_composition


   end subroutine load_sampling_inputs


   subroutine load_nuclear_inputs(inlist, options, ierr)
      character(len=*), intent(in) :: inlist
      type(nuclear_t), intent(inout) :: options
      integer, intent(out) :: ierr
      integer :: unit

      include 'private/nuclear.inc'
   
      namelist /nuclear/  &
            num_special_rate_factors,&
            special_rate_factor,&
            reaction_for_special_factor,&
            filename_of_special_rate,&
            weaklib_blend_hi_Z,&
            T9_weaklib_full_off,&
            T9_weaklib_full_on, &
            T9_weaklib_full_off_hi_Z,&
            T9_weaklib_full_on_hi_Z,&
            use_suzuki_weak_rates,&
            use_3a_fl87,&
            use_special_weak_rates,&
            special_weak_states_file,&
            special_weak_transitions_file,&
            ion_coulomb_corrections,&
            electron_coulomb_corrections,&
            net_reaction_filename,&
            jina_reaclib_filename,&
            jina_reaclib_min_T9,&
            rate_tables_dir,&
            rate_cache_suffix,&
            screening_mode,&
            net_logTcut_lo,&
            net_logTcut_lim,&
            max_logT_for_net,&
            weak_rate_factor,  &
            fe56ec_fake_factor,&
            min_T_for_fe56ec_fake_factor,&
            warn_rates_for_high_temp,&
            max_safe_logT_for_rates


      include '../defaults/nuclear.defaults'

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=nuclear)
      close(unit)


      options% num_special_rate_factors = num_special_rate_factors
      options% special_rate_factor = special_rate_factor
      options% reaction_for_special_factor = reaction_for_special_factor
      options% filename_of_special_rate = filename_of_special_rate
      options% weaklib_blend_hi_Z = weaklib_blend_hi_Z
      options% T9_weaklib_full_off = T9_weaklib_full_off
      options% T9_weaklib_full_on = T9_weaklib_full_on
      options% T9_weaklib_full_off_hi_Z = T9_weaklib_full_off_hi_Z
      options% T9_weaklib_full_on_hi_Z = T9_weaklib_full_on_hi_Z
      options% use_suzuki_weak_rates = use_suzuki_weak_rates
      options% use_3a_fl87 = use_3a_fl87
      options% use_special_weak_rates = use_special_weak_rates
      options% special_weak_states_file = special_weak_states_file
      options% special_weak_transitions_file = special_weak_transitions_file
      options% ion_coulomb_corrections = ion_coulomb_corrections
      options% electron_coulomb_corrections = electron_coulomb_corrections
      options% net_reaction_filename = net_reaction_filename
      options% jina_reaclib_filename = jina_reaclib_filename
      options% jina_reaclib_min_T9 = jina_reaclib_min_T9
      options% rate_tables_dir = rate_tables_dir
      options% rate_cache_suffix = rate_cache_suffix
      options% screening_mode = screening_mode
      options% net_logTcut_lo = net_logTcut_lo
      options% net_logTcut_lim = net_logTcut_lim
      options% max_logT_for_net = max_logT_for_net
      options% weak_rate_factor = weak_rate_factor
      options% fe56ec_fake_factor = fe56ec_fake_factor
      options% min_T_for_fe56ec_fake_factor = min_T_for_fe56ec_fake_factor
      options% warn_rates_for_high_temp = warn_rates_for_high_temp
      options% max_safe_logT_for_rates = max_safe_logT_for_rates

   end subroutine load_nuclear_inputs



end module ctrls

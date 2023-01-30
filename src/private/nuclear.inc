integer :: num_special_rate_factors
real(dp) :: special_rate_factor(max_num_special_rate_factors)
character(len=maxlen_reaction_Name) :: &
    reaction_for_special_factor(max_num_special_rate_factors)
character(len=strlen) :: filename_of_special_rate(max_num_special_rate_factors)
integer :: weaklib_blend_hi_Z
real(dp) :: &
    T9_weaklib_full_off, T9_weaklib_full_on, &
    T9_weaklib_full_off_hi_Z, T9_weaklib_full_on_hi_Z

logical :: use_suzuki_weak_rates, use_3a_fl87

logical :: use_special_weak_rates
character (len=1000) :: special_weak_states_file, special_weak_transitions_file
character (len=strlen) :: ion_coulomb_corrections, electron_coulomb_corrections

character (len=strlen) :: net_reaction_filename, jina_reaclib_filename
real(dp) :: jina_reaclib_min_T9

character (len=strlen) :: rate_tables_dir, rate_cache_suffix

character (len=strlen) :: screening_mode
real(dp) :: net_logTcut_lo, net_logTcut_lim

real(dp) :: max_logT_for_net

real(dp) :: weak_rate_factor,  &
   fe56ec_fake_factor, min_T_for_fe56ec_fake_factor

logical :: warn_rates_for_high_temp 
real(dp) :: max_safe_logT_for_rates 

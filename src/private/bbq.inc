character(len=strlen) :: net_name,iso_list_filename
integer ::  max_steps
real(dp) :: eps,odescal,stptry
character(len=strlen) :: inlist
logical :: use_input_file,use_random_sampling,use_profile, use_hydrostatic
logical :: just_write_isos,write_iso_list

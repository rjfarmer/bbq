program main
   use utils_lib
   use bbq_lib
   use sampler_lib
   use profile_lib
   use random_lib
   use ctrls


   type(bbq_t) :: bbq_in


   call net_setup(bbq_in)


   if(bbq_in% just_write_isos .or. bbq_in% write_iso_list) then
      call write_isos(bbq_in, bbq_in% iso_list_filename)
      if(bbq_in% just_write_isos) return
   end if

   if(bbq_in% use_input_file) then
      call run_sampler_from_file(bbq_in)
   else if(bbq_in% use_random_sampling) then
      call run_random(bbq_in)
   else if (bbq_in% use_profile) then
      call run_profile(bbq_in)
   else
      write(*,*) "Must select one mode to run in"
      call mesa_error(__FILE__,__LINE__)
   end if

end program main

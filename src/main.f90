program main
    use utils_lib
    use bbq_lib
    use sampler_lib
    use profile_lib


    call net_setup()


    if(just_write_isos) then
        call write_isos(iso_list_filename)
        return
    end if

    if(use_input_file) then
        call run_sampler_from_file(inlist_fname)
    else if(use_random_sampling) then
        call run_sampler_random(inlist_fname)
    else if (use_profile) then
        call run_profile(inlist_fname)
    else
        write(*,*) "Must select one mode to run in"
        call mesa_error(__FILE__,__LINE__)
    end if

end program main
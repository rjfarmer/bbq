Embedding
=========

This is a guide for how to embed ``bbq`` in another code. To run ``bbq`` 
yourself, without useing the prewitten drivers, look at the ``src/main.f90`` file.

At the most basic you would do::

        use bbq_lib
        implicit none

        type(bbq_t) :: bbq_in
        type(inputs_t) :: in
        type(outputs_t) :: out


        in% logT = my_logT
        in% logRho = my_logRho
        in% time = my_integration_time
        in% xa = my_starting_xa

        ! This should be done once at the start and handles all setup and loading inlists
        ! After this call bbq_in is read only so could be shared between threads
        call net_setup(bbq_in)

        call do_burn(in, out, bbq_in,  ierr)
        if(ierr/=0) return ! Non-zero signifies some error

        ! Output is stored in the out variable





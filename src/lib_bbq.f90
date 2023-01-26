module bbq_lib
   use chem_lib
   use net_def
   use net_lib
   use rates_lib, only: rates_init
   use rates_def
   use const_lib
   use utils_lib
   use math_lib
   use mtx_def
   use num_lib, only: solver_option
   use mtx_lib, only: decsol_option
   use interp_1d_lib, only: interp_pm
   use interp_1d_def, only: pm_work_size

   implicit none


   character(len=strlen) :: net_name='',screening='',iso_list_filename=''
   integer ::  max_steps
   real(dp) :: weak_rate_factor, eps,odescal,stptry
   character(len=strlen) :: inlist_fname = 'inlist'
   logical :: use_input_file=.false.,use_random_sampling=.false.,use_profile=.false.
   logical :: just_write_isos,write_iso_list

   namelist /bbq/ net_name, screening_mode, weak_rate_factor, &
                  eps, odescal, stptry, max_steps, &
                  use_input_file, use_random_sampling, use_profile,&
                  iso_list_filename, just_write_isos, write_iso_list



   real(dp), dimension(:), pointer :: &
      xin_log =>null(), d_eps_nuc_dx =>null(), dxdt =>null(), &
      d_dxdt_dRho =>null(), d_dxdt_dT,y_log =>null()
   real(dp), pointer :: d_dxdt_dx(:, :)   =>null()

   integer :: eos_handle,net_handle , ind
   integer :: species, num_reactions

   type (Net_Info), target :: netinfo_target
   type (Net_Info), pointer :: netinfo =>null()
   type (Net_General_Info), pointer :: g =>null()

   real(dp), pointer :: rate_factors(:) =>null() ! (num_reactions)
   integer, pointer :: net_reaction_ptr(:) =>null() 

   character (len=64) :: cache_suffix ='0'

   integer, pointer :: reaction_id(:)
   integer, dimension(:), pointer :: net_iso =>null(), chem_id =>null()

   integer,parameter :: num_times=1

   character(len=*),parameter :: which_solver = 'rodas4_solver', &
                              small_mtx_decsol = 'lapack',&
                              large_mtx_decsol = 'bcyclic_dble'

   real(dp), pointer :: burn_work_array(:) =>null(), net_work_array(:) =>null()
         
   real(dp), pointer :: ending_x(:) =>null() ! (species)
   integer :: nfcn    ! number of function evaluations
   integer :: njac    ! number of jacobian evaluations
   integer :: nstep   ! number of computed steps
   integer :: naccpt  ! number of accepted steps
   integer :: nrejct  ! number of rejected steps
   integer :: max_order_used
   integer, parameter :: nwork = pm_work_size

   real(dp) :: logRho, logT, Rho, T, xsum, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT
   real(dp), dimension(:), pointer :: &
      rate_raw =>null(), rate_raw_dT =>null(), rate_raw_dRho =>null(), &
      rate_screened =>null(), rate_screened_dT =>null(), rate_screened_dRho =>null()
   integer :: iounit, decsol_choice, solver_choice

   real(dp), dimension(:), pointer :: times =>null(), dxdt_source_term =>null()
   real(dp), dimension(:), pointer :: &
      log10Ts_f1 =>null(), log10Rhos_f1 =>null(), etas_f1 =>null(), log10Ps_f1 =>null()
   real(dp), dimension(:,:), pointer :: &
      log10Ts_f =>null(), log10Rhos_f =>null(), etas_f =>null(), log10Ps_f =>null()

   real(dp), pointer :: x_initial(:) =>null() ! (species)
   real(dp), pointer :: x_previous(:) =>null() ! (species)
   real(dp), pointer :: pm_work(:) =>null()

   integer :: burn_lwork, net_lwork
   integer :: screening_mode

   real(dp) :: told = 0

   contains

   subroutine read_inlist()
      use rates_lib
      integer :: unit, ierr, status

      ierr = 0

      if (COMMAND_ARGUMENT_COUNT() >= 1) then
         call GET_COMMAND_ARGUMENT(1, inlist_fname, STATUS=status)
         if (status /= 0) inlist_fname = ''
      end if

      open(newunit=unit,file=inlist_fname,status='old',action='read')
      read(unit,nml=bbq)
      close(unit)

      screening_mode = screening_option(screening, ierr)
      if (ierr /= 0) then
         write(*,*) 'read_inlist failed setting screening mode'
         call mesa_error(__FILE__,__LINE__)
      end if     

   end subroutine read_inlist

   subroutine setup_eos()
   ! allocate and load the eos tables
      use eos_def
      use eos_lib
      integer :: ierr
      logical, parameter :: use_cache = .true.
      type(EoS_General_Info),pointer :: rq

      call eos_init('', use_cache, ierr)
      if (ierr /= 0) then
         write(*,*) 'eos_init failed in Setup_eos'
         call mesa_error(__FILE__,__LINE__)
      end if         

      eos_handle = alloc_eos_handle_using_inlist(inlist_fname, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed trying to allocate eos handle'
         call mesa_error(__FILE__,__LINE__)
      end if    

      call get_eos_ptr(eos_handle,rq, ierr )
      if (ierr /= 0) then
         write(*,*) 'failed trying to get eos handle'
         call mesa_error(__FILE__,__LINE__)
      end if  

   end subroutine setup_eos

   subroutine load_libs
      use const_lib
      use const_def, only: mesa_dir
      use chem_lib
      use rates_lib, only: rates_init, rates_warning_init
      integer :: ierr
      character (len=64) :: my_mesa_dir

      my_mesa_dir = ''         
      call const_init(my_mesa_dir,ierr)     
      if (ierr /= 0) then
         write(*,*) 'const_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if        

      call math_init()

      ierr = 0
      call chem_init('isotopes.data', ierr)
      if (ierr /= 0) then
         write(*,*) 'chem_init failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      call rates_init('reactions.list', '', 'rate_tables', .false., .false., &
                     '', '', '',ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

      call rates_warning_init(.true., 10d0)

   end subroutine load_libs

   subroutine net_setup()
      integer :: ierr
      
      call read_inlist()

      call load_libs()

      call setup_eos()

      allocate(net_iso(num_chem_isos), chem_id(num_chem_isos))
      
      call net_init(ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      net_handle = alloc_net_handle(ierr)
      if (ierr /= 0) then
         write(*,*) 'alloc_net_handle failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call net_start_def(net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_start_def failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call read_net_file(net_name, net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'read_net_file failed ', trim(net_name)
         call mesa_error(__FILE__,__LINE__)
      end if

      call net_finish_def(net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_finish_def failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call net_ptr(net_handle, g, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_ptr failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      species = g% num_isos
      num_reactions = g% num_reactions

      allocate(reaction_id(num_reactions))

      call net_setup_tables(net_handle, cache_suffix, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_setup_tables failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call get_chem_id_table(net_handle, species, chem_id, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_chem_id_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call get_net_iso_table(net_handle, net_iso, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_net_iso_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      call get_reaction_id_table(net_handle, num_reactions, reaction_id, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_reaction_id_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call net_set_fe56ec_fake_factor(net_handle, 1d-7,3d9, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_set_fe56ec_fake_factor failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      allocate( &
            xin_log(species), d_eps_nuc_dx(species), y_log(species), &
            dxdt(species), d_dxdt_dRho(species), d_dxdt_dT(species), d_dxdt_dx(species, species))

      netinfo => netinfo_target

      allocate( &
         rate_factors(num_reactions), &
         x_initial(species), x_previous(species), ending_x(species), times(num_times), &
         log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
         etas_f1(4*num_times), log10Ps_f1(4*num_times), &
         stat=ierr)
      if (ierr /= 0) then
         write(*,*) 'allocate failed for Do_One_Zone_Burn'
         call mesa_error(__FILE__,__LINE__)
      end if

      rate_factors(:) = 1

      call get_net_reaction_table_ptr(net_handle, net_reaction_ptr, ierr)
      if (ierr /= 0) then
         write(*,*) 'bad net?  get_net_reaction_table_ptr failed'
         return
      end if

      log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
      log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
      etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)
      log10Ps_f(1:4,1:num_times) => log10Ps_f1(1:4*num_times)

      decsol_choice = decsol_option(small_mtx_decsol, ierr)
      if (ierr /= 0) then
         write(*,*) 'ERROR: unknown small_mtx_decsol ' // trim(small_mtx_decsol)
         return
      end if

      solver_choice = solver_option(which_solver, ierr)
      if (ierr /= 0) then
         write(*,*) 'ERROR: unknown value for which_solver ' // trim(which_solver)
         return
      end if
      
      burn_lwork = net_1_zone_burn_work_size(net_handle,ierr)
      if (ierr /= 0) return
      net_lwork = net_work_size(net_handle,ierr)
      if (ierr /= 0) return
      allocate(net_work_array(net_lwork), burn_work_array(burn_lwork))


   end subroutine net_setup


   subroutine do_burn(burn_logT, burn_logRho, log_burn_tend, xin,&
                     avg_eps_nuc, eps_neu_total, xout, eps_nuc_categories,&
                     ierr)
      use eos_lib
      use eos_def
      use utils_lib

      real(dp),intent(in) :: burn_logT, burn_logRho, log_burn_tend
      real(dp),intent(in) :: xin(species)

      real(dp),intent(out) :: avg_eps_nuc, eps_neu_total
      real(dp), intent(out) :: xout(species)
      real(dp), intent(out), target :: eps_nuc_categories(num_categories)

      integer, intent(out) :: ierr
      
      real(dp) :: burn_eta, starting_logT, burn_tend

      logical,parameter :: burn_dbg=.false.

      integer :: i,r

      real(dp) :: res(num_eos_basic_results) ! (num_eos_basic_results)         
      real(dp) :: d_dlnd(num_eos_basic_results) ! (num_eos_basic_results) 
      real(dp) :: d_dlnT(num_eos_basic_results) ! (num_eos_basic_results)
      real(dp) :: d_dxa(num_eos_basic_results,species) ! (num_eos_d_dxa_results,species)

      logT = burn_logT
      T = exp10(logT)
      logRho = burn_logRho
      Rho = exp10(logRho)
      ierr = 0

      burn_tend = exp10(log_burn_tend)
      
      call eosDT_get( &
         eos_handle, species, chem_id, net_iso, xin, &
         Rho, logRho, T, logT, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)

      burn_eta = res(i_eta)

      if(is_bad(burn_eta)) return

      times(1) = burn_tend
      log10Ts_f(1,1) = logT
      log10Ts_f(2:4,1) = 0
      log10Rhos_f(1,1) = logRho
      log10Rhos_f(2:4,1) = 0
      etas_f(1:1,1) = burn_eta
      etas_f(2:4,1) = 0

      starting_logT = logT
               
      dxdt_source_term => null()

      ending_x = 0d0

      call net_1_zone_burn( &
         net_handle, eos_handle, species, num_reactions, 0d0, burn_tend, xin, &
         num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
         rate_factors, weak_rate_factor, &
         std_reaction_Qs, std_reaction_neuQs, &
         screening_mode,  & 
         stptry, max_steps, eps, odescal, &
         .true., .false., burn_dbg, burn_finish_substep, &
         burn_lwork, burn_work_array, & 
         net_lwork, net_work_array, & 
         ending_x, eps_nuc_categories, &
         avg_eps_nuc, eps_neu_total, &
         nfcn, njac, nstep, naccpt, nrejct, ierr)

      contains
   
      subroutine burn_finish_substep(step, time, y, ierr)
         integer,intent(in) :: step
         real(dp), intent(in) :: time, y(:)
         integer, intent(out) :: ierr 
         ierr = 0            

      end subroutine burn_finish_substep
         
   end subroutine do_burn


   subroutine write_isos(filename)
      integer :: fisos,j
      character(len=*),intent(in) :: filename

      open(newunit=fisos,file=filename,status='replace',action='write')
      do j=1,species
         write(fisos,'(A)') chem_isos% name(g% chem_id(j))
      end do
      close(fisos)

   end subroutine write_isos



end module bbq_lib


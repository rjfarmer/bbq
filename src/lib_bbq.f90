module bbq_lib
   use chem_lib
   use chem_def
   use net_def
   use net_lib
   use rates_lib
   use rates_def
   use const_lib
   use const_def
   use utils_lib
   use math_lib
   use mtx_def
   use eos_def
   use eos_lib
   use num_lib, only: solver_option
   use mtx_lib, only: decsol_option

   use ctrls

   implicit none

   private
   public :: inputs_t,outputs_t
   public :: net_setup, write_isos, write_iso_names, write_output, do_burn

   real(dp), parameter :: ABS_ERR=1d-2 ! Mesa uses this

   type inputs_t
      integer :: id
      real(dp) :: time, logT, logRho
      real(dp), allocatable :: xa(:)
   end type inputs_t

   type outputs_t
      real(dp) :: eps_nuc, eps_neu
      real(dp) :: eps_nuc_categories(num_categories)
      real(dp), allocatable :: xa(:)
      integer :: nfcn    ! number of function evaluations
      integer :: njac    ! number of jacobian evaluations
      integer :: nstep   ! number of computed steps
      integer :: naccpt  ! number of accepted steps
      integer :: nrejct  ! number of rejected steps
   end type outputs_t

   character (len=64),parameter :: cache_suffix ='0'


   contains

   subroutine read_inlist(bbq_in)
      integer :: unit, ierr, status
      character(len=strlen) :: inlist
      type(bbq_t) :: bbq_in

      ierr = 0

      inlist = 'inlist'
      if (COMMAND_ARGUMENT_COUNT() >= 1) then
         call GET_COMMAND_ARGUMENT(1, inlist, STATUS=status)
         if (status /= 0) inlist = ''
      end if

      call load_bbq_inputs(inlist, bbq_in, ierr)
      if(ierr/=0) return

      call load_nuclear_inputs(bbq_in% inlist, bbq_in% nuclear, ierr)
      if(ierr/=0) return

   end subroutine read_inlist

   subroutine setup_eos(bbq_in)
   ! allocate and load the eos tables
      type(bbq_t) :: bbq_in
      integer :: ierr
      logical, parameter :: use_cache = .true.

      call eos_init('', use_cache, ierr)
      if (ierr /= 0) then
         write(*,*) 'eos_init failed in Setup_eos'
         call mesa_error(__FILE__,__LINE__)
      end if         

      bbq_in% state% eos_handle = alloc_eos_handle_using_inlist(bbq_in% inlist, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed trying to allocate eos handle'
         call mesa_error(__FILE__,__LINE__)
      end if    


   end subroutine setup_eos

   subroutine load_libs(bbq_in)
      type(bbq_t),target :: bbq_in
      type(nuclear_t), pointer :: nuc => null()
      integer :: ierr
      character (len=64) :: my_mesa_dir

      nuc => bbq_in% nuclear

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

      call rates_init(nuc% net_reaction_filename, nuc% jina_reaclib_filename,&
                     nuc% rate_tables_dir, &
                     nuc% use_suzuki_weak_rates, nuc% use_special_weak_rates, &
                     nuc% special_weak_states_file, nuc% special_weak_transitions_file,&
                     '',&
                     ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

      call rates_warning_init(nuc% warn_rates_for_high_temp ,nuc% max_safe_logT_for_rates)

      nuc% screening_opt = screening_option(nuc% screening_mode, ierr)
      if (ierr /= 0) then
         write(*,*) 'read_inlist failed setting screening mode'
         call mesa_error(__FILE__,__LINE__)
      end if   

   end subroutine load_libs

   subroutine net_setup(bbq_in)
      type(bbq_t),target :: bbq_in
      integer :: ierr, j
      type(state_t), pointer :: s => null()
      type(nuclear_t), pointer :: nuc => null()
      
      call read_inlist(bbq_in)

      call load_libs(bbq_in)

      call setup_eos(bbq_in)

      s => bbq_in% state
      nuc => bbq_in% nuclear

      allocate(s% net_iso(num_chem_isos), &
               s% chem_id(num_chem_isos))
      
      call net_init(ierr)
      if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
      
      s%  net_handle = alloc_net_handle(ierr)
      if (ierr /= 0) then
         write(*,*) 'alloc_net_handle failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call net_start_def(s% net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_start_def failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call read_net_file(bbq_in% net_name,s%  net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'read_net_file failed ', trim(bbq_in% net_name)
         call mesa_error(__FILE__,__LINE__)
      end if

      call net_finish_def(s% net_handle, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_finish_def failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call net_ptr(s% net_handle, s% g, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_ptr failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      s% species = s% g% num_isos
      s% num_reactions = s% g% num_reactions

      allocate(s% reaction_id(s% num_reactions))
      
      call get_chem_id_table(s% net_handle, s% species, s% chem_id, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_chem_id_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      call get_net_iso_table(s% net_handle, s% net_iso, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_net_iso_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      call get_reaction_id_table(s% net_handle,s% num_reactions,&
                                 s% reaction_id, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_reaction_id_table failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
      s% netinfo => s% netinfo_target

      allocate( &
         s% rate_factors(s% num_reactions),&
         stat=ierr)
      if (ierr /= 0) then
         write(*,*) 'allocate failed for Do_One_Zone_Burn'
         call mesa_error(__FILE__,__LINE__)
      end if

      call set_rate_factors(s, nuc, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in set_rate_factors'
         call mesa_error(__FILE__,__LINE__)
      end if  

      call set_global_nuc_opts(s, nuc)


      call get_net_reaction_table_ptr(s% net_handle, &
                                       s% net_reaction_ptr, ierr)
      if (ierr /= 0) then
         write(*,*) 'bad net?  get_net_reaction_table_ptr failed'
         return
      end if      

      call set_rate_factors(s, nuc, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in set_rate_factors'
         call mesa_error(__FILE__,__LINE__)
      end if  

      call set_global_nuc_opts(s, nuc)

      call net_set_fe56ec_fake_factor(s% net_handle,&
                                      nuc% fe56ec_fake_factor,&
                                      nuc% min_T_for_fe56ec_fake_factor,&
                                      ierr)
      if (ierr /= 0) then
         write(*,*) 'net_set_fe56ec_fake_factor failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      call net_setup_tables(s% net_handle, nuc% rate_cache_suffix, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_setup_tables failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      s% neut_id = -1
      s% prot_id = -1
      do j=1,s% species
         if(trim(chem_isos% name(s% g% chem_id(j)))=='neut') s% neut_id = j
         if(trim(chem_isos% name(s% g% chem_id(j)))=='prot') s% prot_id = j
      end do

   end subroutine net_setup


   subroutine set_rate_factors(state, nuc, ierr)
      type(state_t) :: state
      type(nuclear_t) :: nuc
      integer, intent(out) :: ierr
      integer :: j, i, ir
      logical :: error
            
      state% rate_factors(:) = 1
      if (nuc% num_special_rate_factors <= 0) return
            
      do i=1,nuc% num_special_rate_factors
         if (len_trim(nuc% reaction_for_special_factor(i)) == 0) cycle
         ir = rates_reaction_id(nuc% reaction_for_special_factor(i))
         j = 0
         if (ir > 0) j = state% net_reaction_ptr(ir)
         if (j <= 0) then
            write(*,*) 'Failed to find reaction_for_special_factor ' // &
            trim(nuc% reaction_for_special_factor(i)), &
            j, nuc% special_rate_factor(i)
            error = .true.
            cycle
         end if
         state% rate_factors(j) = nuc% special_rate_factor(i)
         write(*,*) 'set special rate factor for ' // &
               trim(nuc% reaction_for_special_factor(i)), &
               j, nuc% special_rate_factor(i)
      end do

      if(error) call mesa_error(__FILE__,__LINE__)


      call read_rates_from_files(nuc% reaction_for_special_factor, nuc% filename_of_special_rate, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in read_rates_from_files'
         return
      end if
      
   end subroutine set_rate_factors


   subroutine set_global_nuc_opts(state, nuc)
      use rates_def
      type(state_t) :: state
      type(nuclear_t) :: nuc
      integer :: ierr


      if (abs(nuc% T9_weaklib_full_off - T9_weaklib_full_off) > 1d-6) then
         write(*,*) 'set T9_weaklib_full_off', nuc% T9_weaklib_full_off
         T9_weaklib_full_off = nuc% T9_weaklib_full_off
      end if
      
      if (abs(nuc% T9_weaklib_full_on - T9_weaklib_full_on) > 1d-6) then
         write(*,*) 'set T9_weaklib_full_on', nuc% T9_weaklib_full_on
         T9_weaklib_full_on = nuc% T9_weaklib_full_on
      end if
      
      if (nuc% weaklib_blend_hi_Z /= weaklib_blend_hi_Z) then
         write(*,*) 'set weaklib_blend_hi_Z', nuc% weaklib_blend_hi_Z
         weaklib_blend_hi_Z = nuc% weaklib_blend_hi_Z
      end if
      
      if (abs(nuc% T9_weaklib_full_off_hi_Z - T9_weaklib_full_off_hi_Z) > 1d-6) then
         write(*,*) 'set T9_weaklib_full_off_hi_Z', nuc% T9_weaklib_full_off_hi_Z
         T9_weaklib_full_off_hi_Z = nuc% T9_weaklib_full_off_hi_Z
      end if
      
      if (abs(nuc% T9_weaklib_full_on_hi_Z - T9_weaklib_full_on_hi_Z) > 1d-6) then
         write(*,*) 'set T9_weaklib_full_on_hi_Z', nuc% T9_weaklib_full_on_hi_Z
         T9_weaklib_full_on_hi_Z = nuc% T9_weaklib_full_on_hi_Z
      end if

      ! set up coulomb corrections for the special weak rates
      which_mui_coulomb = get_mui_value(nuc% ion_coulomb_corrections)
      which_vs_coulomb = get_vs_value(nuc% electron_coulomb_corrections)
      

      call net_set_logTcut(state% net_handle, nuc% net_logTcut_lo, nuc% net_logTcut_lim, ierr)
      if (ierr /= 0) return

   end subroutine set_global_nuc_opts




   subroutine do_burn(in, out, bbq_in,  ierr)
      type(inputs_t),intent(in) :: in
      type(outputs_t),intent(inout) :: out
      type(bbq_t),target :: bbq_in
      type(nuclear_t), pointer :: nuc=>null()
      integer, intent(out) :: ierr
      
      real(dp) :: eta, logRho, logT, Rho, T, xsum, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, avg_eps_nuc, &
               eps_neu

      logical,parameter :: burn_dbg=.false.

      integer :: i,r, species
      integer,parameter :: num_times=1

      real(dp) :: res(num_eos_basic_results) ! (num_eos_basic_results)         
      real(dp) :: d_dlnd(num_eos_basic_results) ! (num_eos_basic_results) 
      real(dp) :: d_dlnT(num_eos_basic_results) ! (num_eos_basic_results)
      real(dp),allocatable :: d_dxa(:,:) ! (num_eos_d_dxa_results,species)

      real(dp), dimension(:), pointer :: times =>null(), dxdt_source_term =>null()
      real(dp), dimension(:), pointer :: &
          log10Ts_f1 =>null(), log10Rhos_f1 =>null(), etas_f1 =>null(), log10Ps_f1 =>null()
      real(dp), dimension(:,:), pointer :: &
          log10Ts_f =>null(), log10Rhos_f =>null(), etas_f =>null(), log10Ps_f =>null()
      real(dp), dimension(:), pointer :: ending_x=>null()

      ! Validate inputs
      if(in% time < 0) then
         write(*,*) "Got a negative timestep ",in% time
         stop 1
      end if

      if(in% logT > 11) then
         write(*,*) "Temperature is too high for MESA's physics ",in% logT
         stop 1
      end if

      nuc => bbq_in% nuclear

      logT = in% logT
      T = exp10(logT)
      logRho = in% logRho
      Rho = exp10(logRho)
      ierr = 0
      species  = bbq_in% state% species

      
      allocate(d_dxa(num_eos_basic_results,species))
      
      call eosDT_get( &
         bbq_in% state% eos_handle, species, bbq_in% state% chem_id, bbq_in% state% net_iso, in% xa, &
         Rho, logRho, T, logT, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)

      eta = res(i_eta)

      if(is_bad(eta)) return


      allocate( &
         times(num_times), &
         log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
         etas_f1(4*num_times), log10Ps_f1(4*num_times), &
         ending_x(species),&
         stat=ierr)
      if (ierr /= 0) then
         write(*,*) 'allocate failed for Do_One_Zone_Burn'
         call mesa_error(__FILE__,__LINE__)
      end if

      log10Ts_f(1:4,1:num_times) => log10Ts_f1(1:4*num_times)
      log10Rhos_f(1:4,1:num_times) => log10Rhos_f1(1:4*num_times)
      etas_f(1:4,1:num_times) => etas_f1(1:4*num_times)
      log10Ps_f(1:4,1:num_times) => log10Ps_f1(1:4*num_times)

      times(1) = in% time
      log10Ts_f(1,1) = logT
      log10Ts_f(2:4,1) = 0
      log10Rhos_f(1,1) = logRho
      log10Rhos_f(2:4,1) = 0
      etas_f(1:1,1) = eta
      etas_f(2:4,1) = 0
               
      dxdt_source_term => null()

      call net_1_zone_burn( &
         bbq_in% state% net_handle, bbq_in% state% eos_handle, bbq_in% state% species, &
         bbq_in% state% num_reactions, 0d0, times(1), in% xa, &
         num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
         bbq_in% state% rate_factors, nuc% weak_rate_factor, &
         std_reaction_Qs, std_reaction_neuQs, &
         nuc% screening_opt,  & 
         bbq_in% stptry, bbq_in% max_steps, bbq_in% eps, bbq_in% odescal, &
         .true., .false., burn_dbg, burn_finish_substep, &
         ending_x, out% eps_nuc_categories, &
         avg_eps_nuc, eps_neu, &
         out% nfcn, out% njac, out% nstep, &
         out% naccpt, out% nrejct, ierr)

      out% eps_nuc = avg_eps_nuc * in% time
      out% xa = ending_x
      out% eps_neu = eps_neu

      call check_output(in, out, bbq_in, ierr)

      deallocate(times,log10Ts_f1, log10Rhos_f1,etas_f1, log10Ps_f1, &
               ending_x)

      contains
   
      subroutine burn_finish_substep(step, time, y, ierr)
         integer,intent(in) :: step
         real(dp), intent(in) :: time, y(:)
         integer, intent(out) :: ierr 
         ierr = 0            

      end subroutine burn_finish_substep
         
   end subroutine do_burn

   subroutine check_output(in, out, bbq_in, ierr)
      type(inputs_t),intent(in) :: in
      type(outputs_t),intent(inout) :: out
      type(bbq_t),target :: bbq_in
      integer, intent(inout) :: ierr

      if(out% nstep >= bbq_in% max_steps) ierr = -1

      call do_clean1(bbq_in% state% species, out%xa, ierr)

   end subroutine check_output

   subroutine write_isos(bbq_in, filename)
      type(bbq_t) :: bbq_in
      integer :: fisos,j
      character(len=*),intent(in) :: filename

      open(newunit=fisos,file=filename,status='replace',action='write')
      do j=1, bbq_in% state% species
         write(fisos,'(A)') chem_isos% name(bbq_in% state% g% chem_id(j))
      end do
      close(fisos)

   end subroutine write_isos

   subroutine write_iso_names(bbq_in, unit)
      type(bbq_t) :: bbq_in
      integer :: unit, species,j

      do j=1, bbq_in% state% species
         write(unit,'(A)',advance='no') trim(chem_isos% name(bbq_in% state% g% chem_id(j)))//' '
      end do
      write(unit,*) 

   end subroutine write_iso_names


   subroutine write_output(in, out, fout)
      type(inputs_t) :: in
      type(outputs_t) :: out
      integer :: fout,j

      write(fout,'(I0,1X,2(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') in% id, out% eps_nuc, out% eps_neu

      do j=1, size(out% xa)
         write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') out% xa(j)
      end do
      write(fout,*) 

   end subroutine write_output


   subroutine do_clean1(species, xa, ierr)
      use utils_lib
      integer, intent(in) :: species
      real(dp), intent(inout) :: xa(:) ! (species)
      real(dp),parameter :: max_sum_abs = 10d0
      real(dp),parameter :: xsum_tol =1d-2
      integer, intent(out) :: ierr
      integer :: j
      real(dp) :: xsum
      real(dp), parameter :: tiny_x = 1d-99
      character (len=256) :: message
      if (max_sum_abs > 1) then ! check for crazy values
         xsum = sum(abs(xa(1: species)))
         if (is_bad(xsum) .or. xsum > max_sum_abs) then
            ierr = -1
            return
         end if
      end if
      ierr = 0
      do j = 1, species
         if (xa(j) < tiny_x) xa(j) = tiny_x
         if (xa(j) > 1) xa(j) = 1
      end do         
      xsum = sum(xa(1: species))         
      if (abs(xsum-1) > xsum_tol) then
         ierr = -1
         return
      end if
      xa(1: species) = xa(1: species)/xsum
      do j = 1, species
         if (xa(j) < tiny_x) xa(j) = tiny_x
      end do    

   end subroutine do_clean1


end module bbq_lib


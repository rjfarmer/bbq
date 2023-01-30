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


   type inputs_t
      real(dp) :: time, logT, logRho
      real(dp), allocatable :: xa(:)
   end type inputs_t

   type outputs_t
      real(dp) :: eps_nuc, eps_neu
      real(dp) :: eps_nuc_categories(num_categories)
      real(dp), allocatable :: xa(:)
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
      type(bbq_t) :: bbq_in
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

      bbq_in% state% screening_opt = screening_option(bbq_in% screening_mode, ierr)
      if (ierr /= 0) then
         write(*,*) 'read_inlist failed setting screening mode'
         call mesa_error(__FILE__,__LINE__)
      end if   

   end subroutine load_libs

   subroutine net_setup(bbq_in)
      type(bbq_t),target :: bbq_in
      integer :: ierr, j
      type(rate_t), pointer :: s => null()
      
      call read_inlist(bbq_in)

      call load_libs(bbq_in)

      call setup_eos(bbq_in)

      s => bbq_in% state

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

      call net_setup_tables(s% net_handle, cache_suffix, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_setup_tables failed'
         call mesa_error(__FILE__,__LINE__)
      end if
      
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
      
      call net_set_fe56ec_fake_factor(s% net_handle, 1d-7,3d9, ierr)
      if (ierr /= 0) then
         write(*,*) 'net_set_fe56ec_fake_factor failed'
         call mesa_error(__FILE__,__LINE__)
      end if

      s% netinfo => s% netinfo_target

      allocate( &
         s% rate_factors(s% num_reactions),&
         s% ending_x(s% species), &
         stat=ierr)
      if (ierr /= 0) then
         write(*,*) 'allocate failed for Do_One_Zone_Burn'
         call mesa_error(__FILE__,__LINE__)
      end if

      s% rate_factors(:) = 1

      call get_net_reaction_table_ptr(s% net_handle, &
                                       s% net_reaction_ptr, ierr)
      if (ierr /= 0) then
         write(*,*) 'bad net?  get_net_reaction_table_ptr failed'
         return
      end if

      s% neut_id = -1
      s% prot_id = -1
      do j=1,s% species
         if(trim(chem_isos% name(s% g% chem_id(j)))=='neut') s% neut_id = j
         if(trim(chem_isos% name(s% g% chem_id(j)))=='prot') s% prot_id = j
      end do

   end subroutine net_setup


   subroutine do_burn(in, out, bbq_in,  ierr)
      type(inputs_t),intent(in) :: in
      type(outputs_t),intent(out) :: out
      type(bbq_t) :: bbq_in
      integer, intent(out) :: ierr
      
      real(dp) :: eta, logRho, logT, Rho, T, xsum, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, avg_eps_nuc

      logical,parameter :: burn_dbg=.false.

      integer :: i,r
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


      logT = in% logT
      T = exp10(logT)
      logRho = in% logRho
      Rho = exp10(logRho)
      ierr = 0

      allocate(d_dxa(num_eos_basic_results,bbq_in% state% species))
      
      call eosDT_get( &
         bbq_in% state% eos_handle, bbq_in% state% species, bbq_in% state% chem_id, bbq_in% state% net_iso, in% xa, &
         Rho, logRho, T, logT, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)

      eta = res(i_eta)

      if(is_bad(eta)) return


      allocate( &
         times(num_times), &
         log10Ts_f1(4*num_times), log10Rhos_f1(4*num_times), &
         etas_f1(4*num_times), log10Ps_f1(4*num_times), &
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

      bbq_in% state% ending_x = 0d0

      call net_1_zone_burn( &
         bbq_in% state% net_handle, bbq_in% state% eos_handle, bbq_in% state% species, &
         bbq_in% state% num_reactions, 0d0, times(1), in% xa, &
         num_times, times, log10Ts_f1, log10Rhos_f1, etas_f1, dxdt_source_term, &
         bbq_in% state% rate_factors, bbq_in% weak_rate_factor, &
         std_reaction_Qs, std_reaction_neuQs, &
         bbq_in% state% screening_opt,  & 
         bbq_in% stptry, bbq_in% max_steps, bbq_in% eps, bbq_in% odescal, &
         .true., .false., burn_dbg, burn_finish_substep, &
         bbq_in% state% ending_x, out% eps_nuc_categories, &
         avg_eps_nuc, out% eps_neu, &
         bbq_in% state% nfcn, bbq_in% state% njac, bbq_in% state% nstep, &
         bbq_in% state% naccpt, bbq_in% state% nrejct, ierr)

      out% eps_nuc = avg_eps_nuc * in% time
      out% xa = bbq_in% state% ending_x

      contains
   
      subroutine burn_finish_substep(step, time, y, ierr)
         integer,intent(in) :: step
         real(dp), intent(in) :: time, y(:)
         integer, intent(out) :: ierr 
         ierr = 0            

      end subroutine burn_finish_substep
         
   end subroutine do_burn


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


   subroutine write_output(out, fout)
      type(outputs_t) :: out
      integer :: fout,j


      write(fout,'(2(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') out% eps_nuc, out% eps_neu

      do j=1, size(out% xa)
         write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') out% xa(j)
      end do
      write(fout,*) 

   end subroutine write_output


end module bbq_lib


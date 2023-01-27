module random_lib
   use bbq_lib
   use math_lib
   implicit none

   real(dp) :: log_time_min = -12d0, log_time_max = 0d0
   real(dp) :: log_temp_min = 8d0, log_temp_max = 10.0d0
   real(dp) :: log_rho_min = 0d0, log_rho_max= 10d0
   real(dp) :: log_xa_min = -30d0, log_xa_max = 0d0
   real(dp) :: neut_prot_limit_frac = 0.05
   character(len=strlen) :: output_starting_filename='',output_ending_filename=''

   integer :: num_samples=-1, seed = -1

   namelist /random/ log_time_min, log_time_max, log_temp_min, log_temp_max, &
                       log_rho_min,  log_rho_max, log_xa_min, log_xa_max, &
                       neut_prot_limit_frac, &
                       output_starting_filename, output_ending_filename,&
                       num_samples, seed


   real(dp),allocatable :: xin(:)
   real(dp), pointer :: vec(:)
   character(len=4096) :: line
   integer :: neut_id, prot_id
   real(dp) :: time,logt_in,logrho_in,sum,log_time
   integer :: i,j,iostat, k, fin, fout, finput

   real(dp) :: avg_eps_nuc, eps_neu_total
   real(dp), allocatable :: xout(:)
   real(dp), target :: eps_nuc_categories(num_categories)

   logical :: output_in=.true.

   contains

   subroutine run_random(inlist)
      character(len=*),intent(in) :: inlist
      real(dp) :: r
      integer :: total, ierr, n
      integer, allocatable :: seed_arr(:)

      call read_random_inlist(inlist)

      if(seed>0) then
         call RANDOM_SEED(size=n)
         allocate(seed_arr(n))
         seed_arr = seed
         call RANDOM_SEED(put=seed_arr)

         call random_init(.true., .false.)
      else
         call random_init(.false., .false.)
      end if
      call random_number(r)

      call random_setup()
         
      total = 0
      k=0
      do 
         if(num_samples>0 .and. total>num_samples) exit
         total = total+1

         call random_number(r)
         log_time = flat_r(log_time_min,log_time_max,r)

         call random_number(r)
         logt_in = flat_r(log_temp_min, log_temp_max,r)

         call random_number(r)
         logrho_in = flat_r(log_rho_min, log_rho_max,r)

         sum = 0d0
         xin=0d0
         do j=1,species
            call random_number(r)
            if(j==neut_id .or. j==prot_id) then
               ! Limit the amount of free neutrons or protons
               xin(j) = exp10(flat_r(log_xa_min,log10(neut_prot_limit_frac),r))
            else
               xin(j) = exp10(flat_r(log_xa_min,log_xa_max,r))
            end if
            sum = sum + xin(j)
         end do

         xin = xin/sum

         call do_random_burn(ierr)
         if(ierr/=0) return

      end do


   end subroutine run_random


   subroutine read_random_inlist(inlist)
      character(len=*), intent(in) :: inlist
      integer :: unit, ierr, status

      ierr = 0

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=random)
      close(unit)

   end subroutine read_random_inlist

   real(dp) function flat_r(minx,maxx,ran)
      real(dp),intent(in) :: minx,maxx,ran

      flat_r = (ran * (maxx-minx)) + minx

   end function flat_r


   subroutine random_setup()
      integer :: fisos

      allocate(xin(species),vec(species+3),xout(species))

      open(newunit=fin,file=output_starting_filename,status='replace',action='write')
      open(newunit=fout,file=output_ending_filename,status='replace',action='write')


      if(write_iso_list) then
         call write_isos(iso_list_filename)
      end if

      k = 1
         
      neut_id = -1
      prot_id = -1
      do j=1,species
         if(trim(chem_isos% name(g% chem_id(j)))=='neut') neut_id = j
         if(trim(chem_isos% name(g% chem_id(j)))=='prot') prot_id = j
      end do

      !Write header
      write(fout,'(A)',advance='no') '# eps_nuc eps_neu '
      do j=1,species
          write(fout,'(A)',advance='no') trim(chem_isos% name(g% chem_id(j)))//' '
      end do
      write(fout,*)

      !Write header
      write(fin,'(A)',advance='no') '# log_age logt logrho '
      do j=1,species
         write(fin,'(A)',advance='no') trim(chem_isos% name(g% chem_id(j)))//' '
      end do
      write(fin,*)


   end subroutine random_setup


   subroutine do_random_burn(ierr)
      integer :: ierr

      ierr=0
      xout = 0d0
      call do_burn(logt_in, logrho_in, log_time, xin,&
                  avg_eps_nuc, eps_neu_total, xout, eps_nuc_categories, ierr )
      if(ierr/=0) return


      write(fin,'(3(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') log_time, logT, logRho
      write(fout,'(2(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') avg_eps_nuc*10**log_time, eps_neu_total

      do j=1,species
         if(output_in) write(fin,'(1pe26.16)', ROUND='COMPATIBLE',ADVANCE='no') xin(j)
         write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') xout(j)
      end do
      write(fin,*)
      write(fout,*)


      if(mod(k,flush_freq)==0) then
         close(fin)
         close(fout)

         open(newunit=fin,file=output_starting_filename,status='old', position="append", action="write")
         open(newunit=fout,file=output_ending_filename,status='old', position="append", action="write")
         k = 1
      else
         k = k+1
      end if

   end subroutine do_random_burn


end module random_lib
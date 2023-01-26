module sampler
   use bbq_lib
   implicit none

   real(dp) :: log_time_min = -12d0, log_time_max = 0d0
   real(dp) :: log_temp_min = 8d0, log_temp_max = 10.0d0
   real(dp) :: log_rho_min = 0d0, log_rho_max= 10d0
   real(dp) :: log_xa_min = -30d0, log_xa_max = 0d0
   real(dp) :: neut_prot_limit_frac = 0.05
   integer  :: flush_freq = 50
   character(len=strlen) :: input_filename='',output_starting_filename='',output_ending_filename='',&
                            iso_list_filename = ''

   logical :: random_sampling = .false., write_iso_list = .true.
   integer :: num_samples=-1

   namelist /sampling/ log_time_min, log_time_max, log_temp_min, log_temp_max, &
                       log_rho_min,  log_rho_max, log_xa_min, log_xa_max, &
                       neut_prot_limit_frac, flush_freq, &
                       input_filename, output_starting_filename, output_ending_filename,&
                       random_sampling, num_samples, &
                       iso_list_filename


   real(dp),allocatable :: xin(:)
   real(dp), pointer :: vec(:)
   character(len=4096) :: line
   integer :: neut_id, prot_id
   real(dp) :: r, time,logt_in,logrho_in,sum,log_time
   integer :: i,j,n,iostat, k, fin, fout, finput

   real(dp) :: avg_eps_nuc, eps_neu_total
   real(dp), allocatable :: xout(:)
   real(dp), target :: eps_nuc_categories(num_categories)


   contains



   subroutine run_sampler_from_file(inlist)
      character(len=*), intent(in) :: inlist
      integer :: ierr

      call read_sampler_inlist(inlist)

      call sampler_setup()
         
      k = 0
      do 
         xin=0d0

         ! Read data
         read(finput,'(A)',IOSTAT=iostat) line
         if(iostat/=0) stop 
         vec = 0d0
         call str_to_vector(line, vec, n, ierr)
         
         if(n/=species+3) stop 'Bad number of elemenets in row'

         log_time = vec(1)
         logt_in = vec(2)
         logrho_in = vec(3)
         xin(:) = vec(4:species+3)

         call do_sampler_burn(ierr)
         if(ierr/=0) return
      end do

   end subroutine run_sampler_from_file

   subroutine run_sampler_random(inlist)
      character(len=*),intent(in) :: inlist
      real(dp) :: r
      integer :: total, ierr

      call read_sampler_inlist(inlist)

      call random_init(.false., .false.)
      call random_number(r)

      call sampler_setup()
         
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

         call do_sampler_burn(ierr)
         if(ierr/=0) return

      end do


   end subroutine run_sampler_random


   subroutine read_sampler_inlist(inlist)
      character(len=*), intent(in) :: inlist
      integer :: unit, ierr, status

      ierr = 0

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=sampling)
      close(unit)

   end subroutine read_sampler_inlist

   real(dp) function flat_r(minx,maxx,ran)
      real(dp),intent(in) :: minx,maxx,ran

      flat_r = (ran * (maxx-minx)) + minx

   end function flat_r


   subroutine sampler_setup()
      integer :: fisos

      allocate(xin(species),vec(species+3),xout(species))

      open(newunit=fin,file=output_starting_filename,status='replace',action='write')
      open(newunit=fout,file=output_ending_filename,status='replace',action='write')


      if(write_iso_list) then
         open(newunit=fisos,file=iso_list_filename,status='replace',action='write')
         write(fisos,'(A)') chem_isos% name(g% chem_id(j))
         close(fisos)
      end if

      ! Read existing data file
      open(newunit=finput,file=input_filename,status='old',action='read')

      k = 1
         
      neut_id = -1
      prot_id = -1
      do j=1,species
         if(trim(chem_isos% name(g% chem_id(j)))=='neut') neut_id = j
         if(trim(chem_isos% name(g% chem_id(j)))=='prot') prot_id = j
      end do

   end subroutine sampler_setup


   subroutine do_sampler_burn(ierr)
      integer :: ierr

      ierr=0
      call do_burn(logt_in, logrho_in, log_time, xin,&
                  avg_eps_nuc, eps_neu_total, xout, eps_nuc_categories, ierr )
      if(ierr/=0) return


      if(ierr==0) then
         write(fin,'(2(1pe26.16))', ROUND='COMPATIBLE',ADVANCE='no') log_time, logT, logRho
         write(fout,'(22(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') avg_eps_nuc*10**log_time, eps_neu_total
         do j=1,species
            write(fin,'(2(1pe26.16))', ROUND='COMPATIBLE',ADVANCE='no') xin(j)
            write(fout,'(22(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') xout(j)
         end do
         write(fin,*)
         write(fout,*)
      end if

      if(mod(k,flush_freq)==0) then
         close(fin)
         close(fout)
         open(newunit=fin,file=output_starting_filename,status='old', position="append", action="write")
         open(newunit=fout,file=output_ending_filename,status='old', position="append", action="write")
         k = 1
      else
         k = k+1
      end if

   end subroutine do_sampler_burn


end module sampler
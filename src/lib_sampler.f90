module sampler_lib
   use bbq_lib
   use math_lib
   implicit none

   character(len=strlen) :: input_filename='',output_filename=''

   logical :: random_sampling = .false.
   integer :: num_samples=-1

   namelist /sampling/ input_filename, output_filename

   real(dp),allocatable :: xin(:)
   real(dp), pointer :: vec(:)
   character(len=4096) :: line
   integer :: neut_id, prot_id
   real(dp) :: time,logt_in,logrho_in,sum,log_time
   integer :: i,j,n,iostat, k, fin, fout, finput

   real(dp) :: avg_eps_nuc, eps_neu_total
   real(dp), allocatable :: xout(:)
   real(dp), target :: eps_nuc_categories(num_categories)

   logical :: output_in=.true.

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
         if(iostat/=0) return
         vec = 0d0

         call str_to_vector(line, vec, n, ierr)
         if(n/=species+3) call mesa_error(__FILE__,__LINE__, 'Bad number of elemenets in row')

         log_time = vec(1)
         logt_in = vec(2)
         logrho_in = vec(3)
         xin(:) = vec(4:species+3)

         call do_sampler_burn(ierr)
         if(ierr/=0) return
      end do

   end subroutine run_sampler_from_file

   subroutine read_sampler_inlist(inlist)
      character(len=*), intent(in) :: inlist
      integer :: unit, ierr, status

      ierr = 0

      open(newunit=unit,file=inlist,status='old',action='read')
      read(unit,nml=sampling)
      close(unit)

   end subroutine read_sampler_inlist


   subroutine sampler_setup()
      integer :: fisos

      allocate(xin(species),vec(species+3),xout(species))

      open(newunit=fout,file=output_filename,status='replace',action='write')

      if(write_iso_list) then
         call write_isos(iso_list_filename)
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
      xout = 0d0
      call do_burn(logt_in, logrho_in, log_time, xin,&
                  avg_eps_nuc, eps_neu_total, xout, eps_nuc_categories, ierr )
      if(ierr/=0) return


      write(fout,'(2(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') avg_eps_nuc*10**log_time, eps_neu_total

      do j=1,species
         write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') xout(j)
      end do
      write(fout,*)

      if(mod(k,flush_freq)==0) then
         close(fout)
         open(newunit=fout,file=output_filename,status='old', position="append", action="write")
         k = 1
      else
         k = k+1
      end if

   end subroutine do_sampler_burn


end module sampler_lib
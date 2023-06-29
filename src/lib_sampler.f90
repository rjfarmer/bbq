module sampler_lib
   use bbq_lib
   use math_lib
   use iso_fortran_env
   use ctrls

   implicit none

   private
   public :: run_sampler_from_file


   contains

   subroutine run_sampler_from_file(bbq_in)
      type(bbq_t) :: bbq_in
      type(sample_t) :: sample_in
      integer :: ierr, finput,iostat
      type(inputs_t) :: in
      type(outputs_t) :: out 

      character(len=512) :: tmp
      character(len=:),allocatable :: line
      real(dp),pointer :: vec(:) =>null()
      integer :: species, n, fout,line_len, sz, total

      call load_sampling_inputs(bbq_in% inlist, sample_in, ierr)

      call sampler_setup(bbq_in, sample_in)
         
      ! Read existing data file
      open(newunit=finput,file=sample_in% input_filename,status='old',action='read')

      open(newunit=fout,file=sample_in% output_filename,status='old', position="append", action="write")


      species = bbq_in% state% species 
      allocate(in% xa(species),vec(3+species))


      ! Probe the length of the data lines
      line_len =0 
      do
         read(finput,'(A)',iostat=iostat,size=sz,advance='no') tmp
         line_len = line_len+sz
         if(IS_IOSTAT_END(iostat) .or.IS_IOSTAT_EOR(iostat)) exit
      end do
      rewind(finput)
      ! Add buffer padding in case some lines change length
      allocate(character(len=line_len*2) :: line)

      total = 0
      do 
         total = total+1
         in% id = total

         in% xa = 0d0
         ! Read data
         read(finput,'(A)',iostat=iostat) line
         if (iostat == iostat_end) exit

         call str_to_vector(line,vec,n,ierr)

         in% time = vec(1)
         in% logT = vec(2)
         in% logRho = vec(3)

         if(sample_in% uniform_composition) then
            in% xa = 1.d0/size(in% xa)
         else
            in% xa = vec(4:)
         end if

         call do_sampler_burn(bbq_in, sample_in, in, out, fout, ierr)
         if(ierr/=0) return
      end do

   end subroutine run_sampler_from_file

   subroutine sampler_setup(bbq_in, sample_in)
      type(bbq_t) :: bbq_in
      type(sample_t) :: sample_in
      integer :: fout

      open(newunit=fout,file=sample_in% output_filename,status='replace',action='write')

      write(fout,'(A)',advance='no') 'id eps_nuc eps_neu '
      call write_iso_names(bbq_in, fout)
      close(fout)

   end subroutine sampler_setup


   subroutine do_sampler_burn(bbq_in, sample_in, in, out, fout, ierr)
      type(bbq_t) :: bbq_in
      type(sample_t) :: sample_in
      type(inputs_t) :: in
      type(outputs_t) :: out 
      integer :: ierr, fout

      ierr=0
      call do_burn(in, out, bbq_in, ierr)
      if(ierr/=0) return

      call write_output(in, out, fout)

   end subroutine do_sampler_burn


end module sampler_lib

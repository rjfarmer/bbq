module hydrostatic_lib
   use iso_fortran_env
   use bbq_lib
   use ctrls

   implicit none

   private
   public :: run_hydrostatic

   contains

   subroutine run_hydrostatic(bbq_in)
      type(bbq_t) :: bbq_in
      type(hydrostatic_t) :: hydrostatic_in
      type(inputs_t),allocatable :: in(:)
      type(outputs_t) :: out

      integer :: i,j, ierr, fout, fcomp, num_lines
      character(len=8) :: fmt
      character(len=256) :: filename
      real(dp) :: total_time

      call load_hydrostatic_inputs(bbq_in% inlist, hydrostatic_in, ierr)

      call hydrostatic_setup(bbq_in, hydrostatic_in, in, out)

      open(newunit=fout,file=hydrostatic_in% output_filename,status='old', position="append", action="write")

      num_lines = size(in)
      total_time = 0

      ! Add initial line
      call output_hydrostatic(fout, total_time, in(1), out)

      do j=1,num_lines
         in(j)% xa = out% xa
        ! write(*,*) in(j)% xa
         call do_hydrostatic_burn(in(j), out, bbq_in, total_time, fout, ierr )
         if(ierr/=0) return
      end do

   end subroutine run_hydrostatic


   subroutine hydrostatic_setup(bbq_in, hydrostatic_in, in, out)
      type(bbq_t) :: bbq_in
      type(hydrostatic_t) :: hydrostatic_in
      type(inputs_t), allocatable :: in(:)
      type(outputs_t) :: out

      integer :: i,j,fcomp,stat, finput, fout
      integer :: num_lines, species, minTime, maxTime


      if(hydrostatic_in% times_from_file) then
         open(newunit=finput,file=hydrostatic_in% input_filename,action='read')

         num_lines = 0
         do
            read(finput,*,iostat=stat)
            if (stat == iostat_end) exit
            num_lines=num_lines+1
         end do
         rewind(finput)

         allocate(in(num_lines))

         do i=1,num_lines
            read(finput,*) in(i)% time

         end do

         close(finput)

      else

         if(hydrostatic_in% num_times<0) then
            write(*,*) "Bad number of times set ",hydrostatic_in% num_times
         end if

         allocate(in(hydrostatic_in% num_times))

         do i=1, hydrostatic_in% num_times
            in(i)% time = hydrostatic_in% min_time + (i-1)*(hydrostatic_in% max_time - hydrostatic_in% min_time)/(hydrostatic_in% num_times-1)
         end do


         if(hydrostatic_in% log_time) then
            do i=1, hydrostatic_in% num_times
               in(i)% time = exp10(in(i)% time)
            end do
         end if

      end if

      in(:)% logT = hydrostatic_in% logT
      in(:)% logRho = hydrostatic_in% logRho


      species = bbq_in% state% species
      allocate(out% xa(species))
      ! Read in composition
      open(newunit=fcomp,file=hydrostatic_in% input_composition_filename,action='read')

      do j=1,species
         read(fcomp,*) out% xa(j)
      end do
      close(fcomp)

      open(newunit=fout,file=hydrostatic_in% output_filename,status='replace',action='write')

      write(fout,'(A)',advance='no') 'age dt eps_nuc eps_neu '

      !Write header
      call write_iso_names(bbq_in, fout)
      close(fout)

   end subroutine hydrostatic_setup

   subroutine do_hydrostatic_burn(in, out, bbq_in, total_time, fout, ierr)
      type(inputs_t) :: in
      type(outputs_t) :: out
      type(bbq_t) :: bbq_in
      real(dp) :: total_time
      integer :: fout
      integer :: ierr,j

      ierr=0
      call do_burn(in, out, bbq_in, ierr )
      if(ierr/=0) return

      total_time = total_time + in% time
      call output_hydrostatic(fout, total_time,in,out)


   end subroutine do_hydrostatic_burn

   subroutine output_hydrostatic(fout,total_time, in, out)
      type(inputs_t) :: in
      type(outputs_t) :: out
      real(dp) :: total_time
      integer :: fout, j

      write(fout,'(4(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') total_time, in% time, out% eps_nuc, out% eps_neu 

      do j=1, size(out% xa)
            write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') out% xa(j)
      end do
      write(fout,*)

   end subroutine output_hydrostatic


end module hydrostatic_lib

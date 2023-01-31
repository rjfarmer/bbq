module profile_lib
   use iso_fortran_env
   use bbq_lib
   use ctrls

   implicit none

   private
   public :: run_profile

   contains

   subroutine run_profile(bbq_in)
      type(bbq_t) :: bbq_in
      type(profile_t) :: profile_in
      type(inputs_t),allocatable :: in(:)
      type(outputs_t) :: out

      integer :: i,j, ierr, fout, fcomp, num_lines
      character(len=8) :: fmt
      character(len=256) :: filename
      real(dp) :: total_time

      call load_profile_inputs(bbq_in% inlist, profile_in, ierr)

      call profile_setup(bbq_in, profile_in, in)

      open(newunit=fout,file=profile_in% output_filename,status='old', position="append", action="write")

      num_lines = size(in)
      total_time = 0

      ! Initaly setup the composition
      out% xa = in(1)% xa
      ! Add initial line
      call output_profile(fout, total_time, in(1), out)

      do i=1,profile_in% num_loops
         do j=1,num_lines
            in(j)% xa = out% xa
            call do_profile_burn(in(j), out, bbq_in, total_time, fout, ierr )
            if(ierr/=0) return
            write(*,*) "Loop",i,"of",profile_in% num_loops, "zone",j,"of",size(in)   
         end do
         if(.not. profile_in% reflective_boundaries) exit

         do j=num_lines-1,2,-1
            in(j)% xa = out% xa
            call do_profile_burn(in(j), out, bbq_in, total_time, fout, ierr )
            if(ierr/=0) return
            write(*,*) "Loop",i,"of",profile_in% num_loops, "zone",j,"of",size(in)   
         end do

         ! Force a sync
         close(fout)
         open(newunit=fout,file=profile_in% output_filename,status='old', position="append", action="write")

         if(profile_in% write_comp_every_loop) then
            write(fmt,'(I0)') i
            filename = 'comp_'//trim(fmt)//'.txt'
            open(newunit=fcomp,file=filename,status='replace',action="write")

            do j=1,size(in(1)% xa)
                  write(fcomp,'(1pe26.16,1X)', ROUND='COMPATIBLE') out% xa(j)
            end do

            close(fcomp)
         end if

      end do

   end subroutine run_profile


   subroutine profile_setup(bbq_in, profile_in, in)
      type(bbq_t) :: bbq_in
      type(profile_t) :: profile_in
      type(inputs_t), allocatable :: in(:)

      integer :: i,j,fcomp,stat, finput, fout
      integer :: num_lines, species

      open(newunit=finput,file=profile_in% input_filename,action='read')

      num_lines = 0
      do
         read(finput,*,iostat=stat)
         if (stat == iostat_end) exit
         num_lines=num_lines+1
      end do
      rewind(finput)

      allocate(in(num_lines))

      do i=1,num_lines
         read(finput,*) in(i)% time, in(i)% logt, in(i)% logrho
      end do

      close(finput)

      species = bbq_in% state% species
      allocate(in(1)% xa(species))
      ! Read in composition
      open(newunit=fcomp,file=profile_in% input_composition_filename,action='read')

      do j=1,species
         read(fcomp,*) in(1)% xa(j)
      end do
      close(fcomp)

      open(newunit=fout,file=profile_in% output_filename,status='replace',action='write')

      write(fout,'(A)',advance='no') '# age dt logt logrho '

      !Write header
      call write_iso_names(bbq_in, fout)
      close(fout)

   end subroutine profile_setup

   subroutine do_profile_burn(in, out, bbq_in, total_time, fout, ierr)
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
      call output_profile(fout, total_time,in,out)


   end subroutine do_profile_burn

   subroutine output_profile(fout,total_time, in, out)
      type(inputs_t) :: in
      type(outputs_t) :: out
      real(dp) :: total_time
      integer :: fout, j

      write(fout,'(4(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') total_time, in% time, in% logT, in% logRho

      do j=1, size(in% xa)
            write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') out% xa(j)
      end do
      write(fout,*)

   end subroutine output_profile


end module profile_lib

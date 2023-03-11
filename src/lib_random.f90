module random_lib
   use bbq_lib
   use math_lib
   use num_lib
   use ctrls

   implicit none

   private
   public :: run_random

   contains

   subroutine run_random(bbq_in)
      type(bbq_t) :: bbq_in
      type(random_t) :: random_in
      type(inputs_t) :: in
      type(outputs_t) :: out 


      real(dp) :: r,sum
      integer :: total, ierr, n, species, j, fin, i
      integer :: seed_val, fout

      call load_random_inputs(bbq_in% inlist, random_in, ierr)

      if(random_in% seed>0) then
         seed_val = random_in% seed
      else
         call get_seed_for_random(seed_val)
      end if

      ! Burn in random number generator
      do i=1,10
         r = get_dp_uniform_01(seed_val)
      end do

      call random_setup(bbq_in, random_in)
         
      total = 0
      species = bbq_in% state% species
      allocate(in% xa(species))

      open(newunit=fin,file=random_in% output_starting_filename,status='old', position="append", action="write")
      open(newunit=fout,file=random_in% output_ending_filename,status='old', position="append", action="write")


      do 
         if(random_in% num_samples>0 .and. total> random_in% num_samples) exit
         total = total+1

         r = get_dp_uniform_01(seed_val)
         in% time = exp10(flat_r(random_in% log_time_min,random_in% log_time_max,r))

         r = get_dp_uniform_01(seed_val)
         in% logT = flat_r(random_in% log_temp_min, random_in% log_temp_max,r)

         r = get_dp_uniform_01(seed_val)
         in% logRho = flat_r(random_in% log_rho_min, random_in% log_rho_max,r)

         sum = 0d0
         in% xa = 0d0

         do j=1,species
            r = get_dp_uniform_01(seed_val)
            if(j==bbq_in% state% neut_id .or. bbq_in% state% prot_id == j) then
               ! Limit the amount of free neutrons or protons
               in% xa(j) = exp10(flat_r(random_in% log_xa_min,log10(random_in% neut_prot_limit_frac),r))
            else
               in% xa(j) = exp10(flat_r(random_in% log_xa_min,random_in% log_xa_max,r))
            end if
            sum = sum + in% xa(j)
         end do

         in% xa = in% xa/sum

         call do_random_burn(bbq_in, random_in, in, out, fin, fout, ierr)
         if(ierr/=0) return

         if(random_in% sync_freq>0) then
            if(mod(total,random_in% sync_freq)==0) then
               close(fin)
               close(fout)

               open(newunit=fin,file=random_in% output_starting_filename,status='old', position="append", action="write")
               open(newunit=fout,file=random_in% output_ending_filename,status='old', position="append", action="write")
         
            end if
         end if

      end do

   end subroutine run_random

   real(dp) function flat_r(minx, maxx, ran)
      real(dp),intent(in) :: minx, maxx, ran

      flat_r = (ran * (maxx-minx)) + minx

   end function flat_r


   subroutine random_setup(bbq_in, random_in)
      type(bbq_t) :: bbq_in
      type(random_t) :: random_in
      integer :: fin, fout
      integer :: fisos, species

      open(newunit=fin,file=random_in% output_starting_filename,status='replace',action='write')
      open(newunit=fout,file=random_in% output_ending_filename,status='replace',action='write')

      species = bbq_in% state% species

      !Write header
      write(fout,'(A)',advance='no') '# eps_nuc eps_neu '
      call write_iso_names(bbq_in, fout)
      close(fout)

      !Write header
      write(fin,'(A)',advance='no') '# dt logt logrho '
      call write_iso_names(bbq_in, fin)
      close(fin)


   end subroutine random_setup


   subroutine do_random_burn(bbq_in, random_in, in, out, fin, fout, ierr)
      type(bbq_t) :: bbq_in
      type(random_t) :: random_in
      type(inputs_t) :: in
      type(outputs_t) :: out 
      integer :: ierr
      integer :: fin, fout, j


      ierr=0
      call do_burn(in, out, bbq_in, ierr )
      if(ierr/=0) return

      call write_output(out, fout)

      write(fin,'(3(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') in% time, in% logT, in% logRho

      do j=1, size(in% xa)
         write(fin,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') in% xa(j)
      end do


   end subroutine do_random_burn

end module random_lib

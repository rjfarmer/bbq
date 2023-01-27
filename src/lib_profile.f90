module profile_lib
    use iso_fortran_env
    use bbq_lib

    implicit none

    character(len=strlen) :: input_filename='',output_filename='',input_composition_filename=''
    logical :: reflective_boundaries=.true.
    integer :: num_loops = 1


    namelist /profile/ input_filename,output_filename, reflective_boundaries, num_loops,&
                       input_composition_filename


    integer :: finput, fout

    real(dp),allocatable :: xin(:),xout(:),logts(:),logrhos(:),logtimes(:)
    real(dp), pointer :: vec(:)
    integer :: num_lines, loop_iter
    real(dp) :: total_t=0


    contains

    subroutine run_profile(inlist)
        character(len=*) :: inlist
        integer :: i,j, ierr

        call read_profile_inlist(inlist)

        call profile_setup()


        if(reflective_boundaries) then
            loop_iter = 1
            call do_profile_burn(ierr)
            do i=1,num_loops    
                do j=2,num_lines
                    write(*,*) "Loop",i,"of",num_loops, "zone",j,"of",num_lines                    
                    loop_iter = j
                    if(reflective_boundaries) then
                        if(mod(i,2)==0) then
                            loop_iter = num_lines - j +1
                        end if
                    end if
                    call do_profile_burn(ierr)
                    if(ierr/=0) return
                end do
                ! Force a sync
                close(fout)
                open(newunit=fout,file=output_filename,status='old', position="append", action="write")

            end do


        else
            do j=1,num_lines
                loop_iter = j
                call do_profile_burn(ierr)
                if(ierr/=0) return
            end do
        end if

    end subroutine run_profile


    subroutine profile_setup()
        integer :: i,j,fcomp,stat

        allocate(xin(species),vec(3),xout(species))


        if(write_iso_list) then
            call write_isos(iso_list_filename)
         end if

        open(newunit=finput,file=input_filename,action='read')

        num_lines = 0
        do
            read(finput,*,iostat=stat)
            if (stat == iostat_end) exit
            num_lines=num_lines+1
        end do
        rewind(finput)

        allocate(logts(num_lines),logrhos(num_lines),logtimes(num_lines))

        do i=1,num_lines
            read(finput,*) logtimes(i),logts(i),logrhos(i)
         end do

        close(finput)

        ! Read in composition
        open(newunit=fcomp,file=input_composition_filename,action='read')

        do j=1,species
            read(fcomp,*) xin(j)
        end do
        close(fcomp)

        open(newunit=fout,file=output_filename,status='replace',action='write')

   
    end subroutine profile_setup

    subroutine do_profile_burn(ierr)
        integer :: ierr,j
        real(dp) :: avg_eps_nuc, eps_neu_total
        real(dp), target :: eps_nuc_categories(num_categories)

        ierr=0
        xout = 0d0
        call do_burn(logts(loop_iter), logrhos(loop_iter), logtimes(loop_iter), xin,&
                    avg_eps_nuc, eps_neu_total, xout, eps_nuc_categories, ierr )
        if(ierr/=0) return
  
        write(fout,'(4(1pe26.16,1X))', ROUND='COMPATIBLE',ADVANCE='no') total_t,10**logtimes(loop_iter), logts(loop_iter), logrhos(loop_iter)

        do j=1,species
            write(fout,'(1pe26.16,1X)', ROUND='COMPATIBLE',ADVANCE='no') xout(j)
        end do
        write(fout,*)

        total_t = total_t + 10**logtimes(loop_iter)
        xin = xout

    end subroutine do_profile_burn



    subroutine read_profile_inlist(inlist)
        character(len=*), intent(in) :: inlist
        integer :: unit, ierr, status
  
        ierr = 0
  
        open(newunit=unit,file=inlist,status='old',action='read')
        read(unit,nml=profile)
        close(unit)
  
     end subroutine read_profile_inlist

end module profile_lib
module zpares_wrapper
    use zpares
    use iso_c_binding
    type(zpares_prm) :: prm_global

    contains 

    subroutine initialize_zpares_prm(L,N,M,Lmax) bind(c, name = 'initialize_zpares_prm')
        implicit none
        integer,intent(in)::L,N,M,Lmax
        call zpares_init(prm_global)
        prm_global%L = L
        prm_global%N = N
        prm_global%M = M
        prm_global%Lmax = Lmax
        
        return
    end subroutine

    subroutine get_ZPARES_TASK(tasks) bind(c,name='get_ZPARES_TASK')
        implicit none
        integer,intent(out)::tasks(7)
        tasks(1) = ZPARES_TASK_FINISH
        tasks(2) = ZPARES_TASK_FACTO
        tasks(3) = ZPARES_TASK_SOLVE
        tasks(4) = ZPARES_TASK_SOLVE_H
        tasks(5) = ZPARES_TASK_MULT_A
        tasks(6) = ZPARES_TASK_NONE
        tasks(7) = ZPARES_TASK_MULT_B
        write(*,*) tasks
        return
    end

    subroutine wrapper_zpares_zrciheev &
        (ncv,Lmax,mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res,  &
        itask,xs,ws,nc,info) bind(c,name = 'wrapper_zpares_zrciheev')
        implicit none
        integer,intent(in)::mat_size,ncv
        real(8), intent(in) :: emin, emax
        real(8), intent(inout) :: res(1:ncv), eigval(1:ncv)
        complex(8), intent(out)::  X(1:mat_size,1:ncv)
        integer,intent(out)::num_ev(1),info(1),itask(1)
        complex(8),intent(inout)::z(1)
        integer,intent(in)::Lmax
        complex(8),intent(inout)::rwork(mat_size, Lmax),cwork(mat_size, Lmax)
        integer,intent(inout)::xs(1),ws(1),nc(1)
        
        call zpares_zrciheev(prm_global,mat_size, z(1), rwork, cwork, emin, emax, num_ev(1), eigval, X, res, info(1))
        itask(1) = prm_global%itask
        xs(1) = prm_global%xs
        ws(1) = prm_global%ws
        nc(1) = prm_global%nc

    end

    subroutine wrapper_zpares_zrcigegv &
        (ncv,Lmax,mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res,  &
        itask,xs,ws,nc,info) bind(c,name = 'wrapper_zpares_zrcigegv')
        implicit none
        integer,intent(in)::mat_size,ncv
        real(8), intent(in) :: right
        complex(8),intent(in)::left
        real(8), intent(inout) :: res(1:ncv)
        complex(8), intent(out)::  X(1:mat_size,1:ncv),eigval(1:ncv)
        integer,intent(out)::num_ev(1),info(1),itask(1)
        complex(8),intent(inout)::z(1)
        integer,intent(in)::Lmax
        complex(8),intent(inout)::rwork(mat_size, Lmax),cwork(mat_size, Lmax)
        integer,intent(inout)::xs(1),ws(1),nc(1)
        
        call zpares_zrcigegv(prm_global,mat_size, z(1), rwork, cwork, left, right, num_ev(1), eigval, X, res, info(1))
        itask(1) = prm_global%itask
        xs(1) = prm_global%xs
        ws(1) = prm_global%ws
        nc(1) = prm_global%nc

    end

    subroutine wrapper_zpares_zrcigeev &
        (ncv,Lmax,mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res,  &
        itask,xs,ws,nc,info) bind(c,name = 'wrapper_zpares_zrcigeev')
        implicit none
        integer,intent(in)::mat_size,ncv
        real(8), intent(in) :: right
        complex(8),intent(in)::left
        real(8), intent(inout) :: res(1:ncv)
        complex(8), intent(out)::  X(1:mat_size,1:ncv),eigval(1:ncv)
        integer,intent(out)::num_ev(1),info(1),itask(1)
        complex(8),intent(inout)::z(1)
        integer,intent(in)::Lmax
        complex(8),intent(inout)::rwork(mat_size, Lmax),cwork(mat_size, Lmax)
        integer,intent(inout)::xs(1),ws(1),nc(1)
        
        call zpares_zrcigeev(prm_global,mat_size, z(1), rwork, cwork, left, right, num_ev(1), eigval, X, res, info(1))
        itask(1) = prm_global%itask
        xs(1) = prm_global%xs
        ws(1) = prm_global%ws
        nc(1) = prm_global%nc

    end

    integer function zpares_get_ncv_wrapper() bind(c, name = 'zpares_get_ncv_wrapper')
        implicit none

        zpares_get_ncv_wrapper = zpares_get_ncv(prm_global)

        return 
    end function
end module zpares_wrapper

#include "def_macros.h"

subroutine MACRO_INSERT_PRFX(zpares_,rcigegv) &
(prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
  use zpares_aux
  implicit none
#ifdef REALMAT
  logical, parameter :: real_matrix = .true.
#else
  logical, parameter :: real_matrix = .false.
#endif
#include "set_rule.f90"
  optional set_rule
  type(zpares_prm), intent(inout), target :: prm
  integer, intent(in) :: nrow_local
  integer, intent(inout) :: info
  integer, intent(inout) :: num_ev
  REAL_TYPE, intent(in) :: right
  REAL_TYPE, intent(out) :: res(*)
  COMPLEX_TYPE, intent(in) :: left
  COMPLEX_TYPE, intent(inout) :: z, cwork(nrow_local,*), eigval(*)
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  
  integer :: i, j, k, N, M, low_comm, low_comm_rank, low_comm_size, high_comm, high_comm_rank, high_comm_size, comm_self
  integer :: Lmax, Lstep, Lold, rhs_size, counter, imax, quad_type, extract, LmaxM, LM, jx, jw, ierr, infola
  integer, allocatable :: real_unsym_map(:)
  integer, pointer :: itask, L, ws, xs, nc, iter, state, quad_idx, x_offset, x_ncol, mode, num_ev_est, num_basis
  logical :: calc_res, standard, Hermitian, B_pos_def, real_unsym
  REAL_TYPE :: asp_ratio, delta, res_max, tol
  COMPLEX_TYPE :: zeta, weight, tmpc
  MATRIX_TYPE :: m_num_ev_est, dummy_U(1,1), dummy_VT(1,1)
  integer, allocatable :: indices(:)
  REAL_TYPE, allocatable :: AXnorm(:), BXnorm(:), sigma(:), tmprvec(:)
  MATRIX_TYPE, allocatable :: tmp_m1d(:), tmp_mat(:,:), tmp_mat2(:,:), H0(:,:)
  MATRIX_TYPE, pointer :: Mu(:,:)
  logical, allocatable :: flags(:)

  if ( prm%state == ZPARES_STATE_INIT1 .and. prm%iter == 0 ) then
     info = ZPARES_INFO_UNDETERMINED
     call initialize(prm)
     if ( .not. check_inputs(prm) ) then
        info = ZPARES_INFO_INVALID_PARAM
        goto 999 ! Exception
     end if
     if ( prm%user_source ) prm%mode = ZPARES_MODE_SS        
  end if
  
  ! set invariable variables
  N = prm%N
  M = prm%M
  Lmax = prm%Lmax
  Lstep = prm%Lstep
  low_comm = prm%low_comm
  high_comm = prm%high_comm
  imax = prm%imax
  quad_type = prm%quad_type
  calc_res = prm%calc_res
  extract = prm%extract
  delta = prm%delta
  asp_ratio = prm%asp_ratio
  tol = prm%tol
  Hermitian = prm%Hermitian
  standard = prm%standard
  B_pos_def = prm%B_pos_def
  LmaxM = Lmax*M

  ! set variable variables
  itask => prm%itask
  L => prm%L
  iter => prm%iter
  ws => prm%ws
  xs => prm%xs
  nc => prm%nc
  state => prm%state
  quad_idx => prm%quad_idx
  x_offset => prm%x_offset
  x_ncol => prm%x_ncol
  
  mode => prm%mode
  num_ev_est => prm%num_ev_est
  num_basis => prm%num_basis

  call get_rank_and_size(low_comm, low_comm_rank, low_comm_size)
  call get_rank_and_size(high_comm, high_comm_rank, high_comm_size)
  
  select case(state)
  case (ZPARES_STATE_INIT1)
     
     if ( iter == 0 .and. .not. prm%user_source ) then     
        rhs_size = Lmax
     else
        rhs_size = L
     end if
     select case(itask)
     case(ZPARES_TASK_NONE)          

        if ( iter == 0 .and. .not. prm%user_source ) then
           call start_timer(prm%timer_rand)
           call create_hutch_samples(rwork, nrow_local, Lmax, low_comm_rank)
           call stop_timer(prm%timer_rand)
        else
           if ( num_ev > rhs_size ) then
              allocate(tmp_mat(num_ev, rhs_size))
              call create_rand_matrix(tmp_mat, num_ev, rhs_size, 0)
              call GEMM_ALLREDUCE('N', 'N', nrow_local, rhs_size, num_ev, ONE_M, X, nrow_local, &
                   tmp_mat, num_ev, ZERO_M, rwork, nrow_local, ierr)
              deallocate(tmp_mat)
           else if ( num_ev == rhs_size ) then
              rwork(:,1:rhs_size) = X(:,1:rhs_size)
           else
              info = ZPARES_INFO_NOSUF_USER_SRC
              goto 999
           end if
        end if
        if ( .not. standard ) then
           call para_range(rhs_size, high_comm_size, high_comm_rank, ws, nc)
           X(:,1:rhs_size) = rwork(:,1:rhs_size)
           xs = ws
           rwork(:,1:rhs_size) = ZERO_M ! for allreduce
           if ( nc /= 0 ) then             
              itask = ZPARES_TASK_MULT_B
              call start_timer(prm%timer_mult_B)
           else
              state = ZPARES_STATE_INIT2
           end if
        else
           state = ZPARES_STATE_INIT2
        end if

     case(ZPARES_TASK_MULT_B)
        itask = ZPARES_TASK_NONE
        call stop_timer(prm%timer_mult_B)
        state = ZPARES_STATE_INIT2
     end select
     
  case(ZPARES_STATE_INIT2)

     if ( .not. standard ) then
        if ( iter == 0 ) then
           rhs_size = Lmax
        else
           rhs_size = L
        end if
        call ALLREDUCE_SUM_2D(rwork, 1, nrow_local, rhs_size, ierr, high_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
     end if
     X(:,1:LmaxM) = ZERO_M
     call zpares_rci_sub_get_Mu(prm, Mu)
     Mu(:,:) = ZERO_M
     ws = 1
     nc = L
     state = ZPARES_STATE_LINSOLVE       
     
  case(ZPARES_STATE_LINSOLVE) 
     
     call zpares_rci_sub_linsolve(prm, nrow_local, z, rwork, cwork, left, right, X, info, set_rule)
     if ( info /= ZPARES_INFO_UNDETERMINED ) goto 999
     
  case(ZPARES_STATE_ALLREDUCE)
     
     call start_timer(prm%timer_sum)
     do i = 1, M
        call ALLREDUCE_SUM_2D(X, (i-1)*Lmax+ws, nrow_local, nc, ierr, high_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           goto 999
        end if
     end do
     if ( extract == ZPARES_EXTRACT_EM .or. mode /= ZPARES_MODE_SS ) then
        call zpares_rci_sub_get_Mu(prm, Mu)
        do i = 1, 2*M
           call ALLREDUCE_SUM_2D(Mu, (i-1)*Lmax+ws, Lmax, nc, ierr, high_comm)
           if ( ierr /= 0 ) then
              info = ZPARES_INFO_MPI_ERROR
              goto 999
           end if
        end do
     end if
     call stop_timer(prm%timer_sum)

     select case(mode)
     case(ZPARES_MODE_STC)
        state = ZPARES_STATE_STC
     case(ZPARES_MODE_LADD)
        if ( L /= Lmax ) then
           state = ZPARES_STATE_LADD
        else
           state = ZPARES_STATE_SS_SVD
        end if
     case(ZPARES_MODE_SS)
        state = ZPARES_STATE_SS_SVD
     end select
     
  case(ZPARES_STATE_STC)

     call create_hutch_samples(rwork, nrow_local, L, low_comm_rank)
     allocate(tmp_m1d(L))
     call DOT_ALLREDUCE(rwork, X, nrow_local, L, low_comm, tmp_m1d, ierr)
     if ( ierr /= 0 ) then
        info = ZPARES_INFO_MPI_ERROR
        goto 999
     end if
     m_num_ev_est = ZERO_M
     do i = 1, L
        m_num_ev_est = m_num_ev_est + tmp_m1d(i)
     end do
     deallocate(tmp_m1d)
     m_num_ev_est = m_num_ev_est / L
     num_ev_est = nint(abs(m_num_ev_est))
     ! num_ev_est = nint(real(m_num_ev_est,kind(ZERO_R)))
     Lold = L
     L = min(ceiling(2*ONE_R*num_ev_est / M), Lmax)
     if ( Lold /= Lmax .and. L > Lold ) then
        ws = Lold + 1
        nc = L - Lold
        mode = ZPARES_MODE_LADD
        state = ZPARES_STATE_LINSOLVE
     else
        L = Lold
        mode = ZPARES_MODE_SS
        state = ZPARES_STATE_SS_SVD
     end if
     
  case(ZPARES_STATE_LADD)
     
     call zpares_rci_sub_get_Mu(prm, Mu)
     
     if ( L /= Lmax ) then
        j = L + 1
        do i = 2, 2*M
           do k = (i-1)*Lmax+1, (i-1)*Lmax+L
              Mu(:,j) = Mu(:,k)
              j = j + 1
           end do
        end do
     end if
     LM = L*M
     allocate(H0(LM,LM))
     
     call block_Hankel(Lmax, L, M, 0, Mu, H0)
     allocate(sigma(L*M))
     call start_timer(prm%timer_serial_SVD)
     call serial_SVD('N', LM, LM, H0, LM, delta, sigma, dummy_U, 1, dummy_VT, 1, infola, num_basis)
     ! if ( infola /= 0 ) then
     !    info = ZPARES_INFO_LAPACK_ERROR
     !    goto 999
     ! end if
     call stop_timer(prm%timer_serial_SVD)
     deallocate(H0)
     if ( sigma(1) < delta .or. num_basis == LM ) then
        Lold = L
        L = min(L + Lstep, Lmax)
        ws = Lold + 1
        nc = L - Lold
        state = ZPARES_STATE_LINSOLVE
     else          
        mode = ZPARES_MODE_SS
        state = ZPARES_STATE_SS_SVD
     end if
     deallocate(sigma)
     
  case(ZPARES_STATE_SS_SVD)
     
     cwork(:,1:L) = X(:,1:L)
     call zpares_rci_sub_svd(prm, nrow_local, rwork, X, info)
     if ( info /= ZPARES_INFO_UNDETERMINED ) goto 999
     if ( ( prm%sig_val(1) < delta .or. num_basis == L*M .or. prm%force_inner ) .and. iter < imax ) then
        X(:,1:L) = cwork(:,1:L)
        num_ev = L
        iter = iter + 1
        prm%iter_info(iter) = 1
        itask = ZPARES_TASK_NONE
        state = ZPARES_STATE_INIT1
        ! Do iterative refinedent
     else
        state = ZPARES_STATE_INIT_REDUCED_EIG
     end if
     
  case(ZPARES_STATE_INIT_REDUCED_EIG, ZPARES_STATE_FORM_REDUCED_EIG, ZPARES_STATE_SOLVE_REDUCED_EIG)

     select case(extract)
     case(ZPARES_EXTRACT_RR)
        call zpares_rci_sub_Rayleigh_Ritz &
             (prm, nrow_local, rwork, eigval, X, info)
        if ( info /= ZPARES_INFO_UNDETERMINED ) goto 999
     case(ZPARES_EXTRACT_EM)
        call zpares_rci_sub_Hankel_method &
             (prm, nrow_local, rwork, left, right, eigval, X, info, set_rule)
        if ( info /= ZPARES_INFO_UNDETERMINED ) goto 999
     end select

  case(ZPARES_STATE_FIN_REDUCED_EIG)
     
     itask = ZPARES_TASK_NONE
     if ( calc_res ) then
        state = ZPARES_STATE_CALC_RES0
     else
        state = ZPARES_STATE_FINISH
     end if     
     
  case(ZPARES_STATE_CALC_RES0, ZPARES_STATE_CALC_RES1, ZPARES_STATE_CALC_RES2)
     
     call zpares_rci_sub_calc_res(prm, nrow_local, rwork, cwork, eigval, X, res, info)
     if ( info /= ZPARES_INFO_UNDETERMINED ) goto 999

  case(ZPARES_STATE_CHECK_OUTER)
     
     if ( iter < imax .and. num_basis >= L ) then
        allocate(flags(num_basis))
        call inside_ellipse(left, right, asp_ratio, num_basis, eigval, flags, counter)
        res_max = ZERO_R
        do i = 1, num_basis
           if (flags(i)) then
              res_max = max(res_max, res(i))
           end if
        end do
        flags(1:num_basis) = flags(1:num_basis) .and. (prm%indi_spu(1:num_basis) > prm%spu_thres)
        counter = count(flags)
        if ( res_max > tol ) then
           if ( counter > L ) then
              call packing(num_basis, flags, eigval, X, nrow_local, res, prm%indi_spu)
              num_ev = counter
           else
              num_ev = num_basis
           end if
           iter = iter + 1
           prm%iter_info(iter) = 2
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_INIT1
           ! Do iterative refinement
        else
           state = ZPARES_STATE_FINISH
        end if
        deallocate(flags)     
     else
        state = ZPARES_STATE_FINISH
     end if

  case(ZPARES_STATE_FINISH)

     num_ev = num_basis     
     if ( prm%trim_out .and. prm%quad_type == ZPARES_QUAD_ELL_TRAP ) then
        allocate(flags(num_ev))
        call inside_ellipse(left, right, asp_ratio, num_basis, eigval, flags, j)
        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
        num_ev = j
        deallocate(flags)
     end if
     if ( prm%trim_spu ) then
        allocate(flags(num_ev))
        flags(1:num_ev) = prm%indi_spu(1:num_ev) > prm%spu_thres
        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
        num_ev = count(flags)
        deallocate(flags)
     end if
     if ( prm%trim_res ) then
        allocate(flags(num_ev))
        flags(1:num_ev) = res(1:num_ev) > prm%tol
        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
        num_ev = count(flags)
        deallocate(flags)
     end if
     itask = ZPARES_TASK_FINISH
     info = ZPARES_INFO_SUCCESS
     call finalize(prm)

  end select

  return

999 continue
  ! Exception handling
  num_ev = 0
  prm%itask = ZPARES_TASK_FINISH
  call finalize(prm)

  contains
    
    subroutine initialize(prm)
      implicit none
      type(zpares_prm) :: prm
      
      allocate(prm%sig_val(prm%Lmax*prm%M))
      allocate(prm%indi_spu(prm%Lmax*prm%M))
      if ( prm%imax > 0 ) then
         allocate(prm%iter_info(prm%imax))
         prm%iter_info(:) = 0
      end if
      allocate(prm%ptrs%MACRO_ADD_PRFX(projAB)(prm%Lmax*prm%M, 2*prm%Lmax*prm%M))
      allocate(prm%ptrs%MACRO_ADD_PRFX(Mu)(prm%Lmax, 2*prm%Lmax*prm%M))
    end subroutine initialize
    
    subroutine finalize(prm)
      implicit none
      type(zpares_prm) ::  prm
  
      deallocate(prm%ptrs%MACRO_ADD_PRFX(projAB))
      deallocate(prm%ptrs%MACRO_ADD_PRFX(Mu))
    end subroutine finalize

end subroutine


#ifdef REALMAT
subroutine MACRO_INSERT_PRFX(zpares_,rcisygv) &
#else
subroutine MACRO_INSERT_PRFX(zpares_,rcihegv) &
#endif
(prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
  implicit none
#include "set_rule.f90"
  optional set_rule
  type(zpares_prm), intent(inout) :: prm
  integer, intent(in) :: nrow_local
  integer, intent(out) :: info
  integer, intent(inout) :: num_ev
  REAL_TYPE, intent(in) :: left, right
  REAL_TYPE, intent(out) :: res(*), eigval(*)
  COMPLEX_TYPE, intent(inout) :: z, cwork(nrow_local,*)
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)

  COMPLEX_TYPE, allocatable :: eigval_c(:)

  integer :: i, LmaxM

  prm%Hermitian = .true.
  prm%B_pos_def = .true.
  LmaxM = prm%Lmax*prm%M
  allocate(eigval_c(LmaxM))
  eigval_c(1:LmaxM) = cmplx(eigval(1:LmaxM), ZERO_R, kind(ZERO_R))

  call MACRO_INSERT_PRFX(zpares_,rcigegv) &
  (prm, nrow_local, z, rwork, cwork &
  , cmplx(left, ZERO_R, kind(ZERO_R)), right, num_ev, eigval_c, X, res, info, set_rule)

  eigval(1:LmaxM) = real(eigval_c(1:LmaxM), kind(ZERO_R))
  deallocate(eigval_c)
end subroutine


#ifdef REALMAT
subroutine MACRO_INSERT_PRFX(zpares_,rcisyev) &
#else
subroutine MACRO_INSERT_PRFX(zpares_,rciheev) &
#endif
(prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
  implicit none
#include "set_rule.f90"
  optional set_rule
  type(zpares_prm), intent(inout) :: prm
  integer, intent(in) :: nrow_local
  integer, intent(out) :: info
  integer, intent(inout) :: num_ev
  REAL_TYPE, intent(in) :: left, right
  REAL_TYPE, intent(out) :: res(*), eigval(*)
  COMPLEX_TYPE, intent(inout) :: z, cwork(nrow_local,*)
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  
  prm%standard = .true.
#ifdef REALMAT
  call MACRO_INSERT_PRFX(zpares_,rcisygv) &
#else
  call MACRO_INSERT_PRFX(zpares_,rcihegv) &
#endif
  (prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)

end subroutine


subroutine MACRO_INSERT_PRFX(zpares_,rcigeev) &
(prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
  implicit none
#include "set_rule.f90"
  optional set_rule
  type(zpares_prm), intent(inout), target :: prm
  integer, intent(in) :: nrow_local
  integer, intent(out) :: info
  integer, intent(inout) :: num_ev
  REAL_TYPE, intent(in) :: right
  REAL_TYPE, intent(out) :: res(*)
  COMPLEX_TYPE, intent(in) :: left
  COMPLEX_TYPE, intent(inout) :: z, cwork(nrow_local,*), eigval(*)
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)

  prm%standard = .true.
  prm%B_pos_def = .true.
  call MACRO_INSERT_PRFX(zpares_,rcigegv) &
  (prm, nrow_local, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
end subroutine

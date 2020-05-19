#include "def_macros.h"

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_get_projAB)(prm, projAB)
  implicit none
  type(zpares_prm), intent(inout) :: prm
  MATRIX_TYPE, pointer :: projAB(:,:)
  
  projAB => prm%ptrs%MACRO_ADD_PRFX(projAB)
end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_get_Mu)(prm, Mu)
  implicit none
  type(zpares_prm), intent(inout) :: prm
  MATRIX_TYPE, pointer :: Mu(:,:)
  
  Mu => prm%ptrs%MACRO_ADD_PRFX(Mu)
end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_linsolve) &
     (prm, nrow_local, z, rwork, cwork, left, right, X, info, set_rule)
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
  COMPLEX_TYPE, intent(out) :: z
  COMPLEX_TYPE, intent(in) :: left
  REAL_TYPE, intent(in) :: right
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  COMPLEX_TYPE, intent(inout) :: cwork(nrow_local,*)
  integer, intent(inout) :: info
  integer :: i, j, k, L, N, M, Lmax, LM, jw, jx, quad_type, num_linsolve_global, ierr
  integer, pointer :: itask, state, quad_idx, ws, nc, iter
  MATRIX_TYPE, allocatable :: tmp_mat(:,:), H0(:,:)
  REAL_TYPE :: asp_ratio, v1or2
  REAL_TYPE, allocatable :: sigma(:)
  COMPLEX_TYPE, allocatable :: tmp_mat_c(:,:), c_rwork(:,:)
  MATRIX_TYPE, pointer :: Mu(:,:)
  COMPLEX_TYPE :: weight, zeta, tmpc, dummy
  integer :: low_comm, low_comm_rank, low_comm_size, high_comm, high_comm_rank, high_comm_size
  logical :: sym_contour

  N = prm%N
  M = prm%M
  Lmax = prm%Lmax
  low_comm = prm%low_comm
  high_comm = prm%high_comm
  call get_rank_and_size(low_comm, low_comm_rank, low_comm_size)
  call get_rank_and_size(high_comm, high_comm_rank, high_comm_size)
  itask => prm%itask
  state => prm%state
  quad_idx => prm%quad_idx
  ws => prm%ws
  nc => prm%nc
  iter => prm%iter
  quad_type = prm%quad_type
  asp_ratio = prm%asp_ratio
  sym_contour = zpares_rci_sub_get_sym_contour(prm, left)
  if ( sym_contour .and. ( real_matrix .or. prm%Hermitian ) ) then
     num_linsolve_global = N / 2
  else
     num_linsolve_global = N
  end if

  select case(itask)
  case(ZPARES_TASK_NONE)
     quad_idx = 1 + high_comm_rank
     if ( quad_idx <= num_linsolve_global ) then
        select case(quad_type)
        case(ZPARES_QUAD_ELL_TRAP)
           call quad_ell_trap(left, right, N, asp_ratio, quad_idx, zeta, weight, z)
        case(ZPARES_QUAD_USER)
           if ( present(set_rule) ) then
             call set_rule(ZPARES_QU_MODE_RULE, prm%N, quad_idx, left, right, z, weight, zeta, dummy)
           else
              !! error TODO: error handling
           end if
        end select
        if ( high_comm_size >= num_linsolve_global .and. iter >= 1 ) then
           itask = ZPARES_TASK_SOLVE
           call start_timer(prm%timer_solve)
        else
           itask = ZPARES_TASK_FACTO
           call start_timer(prm%timer_fact)
        end if
     else
        state = ZPARES_STATE_ALLREDUCE
     end if
  case(ZPARES_TASK_FACTO)
     call stop_timer(prm%timer_fact)
     cwork(:,ws:ws+nc-1) = rwork(:,ws:ws+nc-1)
     itask = ZPARES_TASK_SOLVE
     call start_timer(prm%timer_solve)
  case(ZPARES_TASK_FACTO_H)
     call stop_timer(prm%timer_fact_H)
     cwork(:,ws:ws+nc-1) = rwork(:,ws:ws+nc-1)
     itask = ZPARES_TASK_SOLVE_H
     call start_timer(prm%timer_solve_H)
  case(ZPARES_TASK_SOLVE, ZPARES_TASK_SOLVE_H)
     if ( itask == ZPARES_TASK_SOLVE ) then
        call stop_timer(prm%timer_solve)
     else
        call stop_timer(prm%timer_solve_H)
     end if
     select case(quad_type)
     case(ZPARES_QUAD_ELL_TRAP)                      
        call quad_ell_trap(left, right, N, asp_ratio, quad_idx, zeta, weight, z)
     case(ZPARES_QUAD_USER)
        if ( present(set_rule) ) then
           call set_rule(ZPARES_QU_MODE_RULE, prm%N, quad_idx, left, right, z, weight, zeta, dummy)
        else
           !! error TODO: error handling
        end if
     end select
     if ( itask == ZPARES_TASK_SOLVE_H ) then
        weight = conjg(weight)
        zeta = conjg(zeta)
     end if
     do i = 1, M
        call start_timer(prm%timer_sum)
        if ( real_matrix ) then
           if ( sym_contour ) then
              v1or2 = ONE_R + ONE_R
           else
              v1or2 = ONE_R              
           end if
           tmpc = v1or2*weight*zeta**(i-1)
           do j = 1, nc
              jx = (i-1)*Lmax+ws+j-1
              jw = ws+j-1
              X(:,jx) = X(:,jx) + real(tmpc*cwork(:,jw), kind(ZERO_R))
           end do
        else
           tmpc = weight*zeta**(i-1)
           do j = 1, nc
              jx = (i-1)*Lmax+ws+j-1
              jw = ws+j-1
              X(:,jx) = X(:,jx) + tmpc*cwork(:,jw)
           end do
        end if
        call stop_timer(prm%timer_sum)
     end do
     if ( prm%extract == ZPARES_EXTRACT_EM .or. prm%mode /= ZPARES_MODE_SS ) then
        allocate(tmp_mat_c(Lmax,nc), c_rwork(nrow_local,Lmax))              
        c_rwork(:,1:Lmax) = rwork(:,1:Lmax)
        call start_timer(prm%timer_gemm_reduce)
        call GEMM_ALLREDUCE('C', 'N', Lmax, nc, nrow_local, ONE_C, c_rwork, nrow_local &
             , cwork(:,ws:ws+nc-1), nrow_local, ZERO_C, tmp_mat_c, Lmax, ierr, low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        call stop_timer(prm%timer_gemm_reduce)
        deallocate(c_rwork)
        call zpares_rci_sub_get_Mu(prm, Mu)
        call start_timer(prm%timer_sum)
        do i= 1, 2*M
           if ( real_matrix ) then
              if ( sym_contour ) then
                 v1or2 = ONE_R + ONE_R
              else
                 v1or2 = ONE_R
              end if
              tmpc = v1or2*weight*zeta**(i-1)
              do j = 1, nc
                 jx = (i-1)*Lmax+ws+j-1
                 jw = ws+j-1
                 Mu(:,jx) = Mu(:,jx) + real(tmpc*tmp_mat_c(:,jw), kind(ZERO_R))
              end do
           else
              tmpc = weight*zeta**(i-1)
              do j = 1, nc
                 jx = (i-1)*Lmax+ws+j-1
                 jw = ws+j-1
                 Mu(:,jx) = Mu(:,jx) + tmpc*tmp_mat_c(:,jw)
              end do
           end if
        end do
        call stop_timer(prm%timer_sum)
        deallocate(tmp_mat_c)
     end if
     if ( itask == ZPARES_TASK_SOLVE_H .or. real_matrix .or. ( .not. ( sym_contour .and. prm%Hermitian ) ) ) then
        if ( quad_idx + high_comm_size > num_linsolve_global ) then
           ! end of solve
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_ALLREDUCE
        else
           quad_idx = quad_idx + high_comm_size
           select case(quad_type)
           case(ZPARES_QUAD_ELL_TRAP)
              call quad_ell_trap(left, right, N, asp_ratio, quad_idx, zeta, weight, z)
           case(ZPARES_QUAD_USER)
              if ( present(set_rule) ) then
                 call set_rule(ZPARES_QU_MODE_RULE, prm%N, quad_idx, left, right, z, weight, zeta, dummy)
              else
                 !! error TODO: error handling
              end if
           end select
           ! zeta and weight will not be used
           if ( high_comm_size >= num_linsolve_global .and. iter >= 1 ) then
              itask = ZPARES_TASK_SOLVE
              call start_timer(prm%timer_solve)
           else
              itask = ZPARES_TASK_FACTO
              call start_timer(prm%timer_fact)
           end if
        end if
     else
        if ( high_comm_size >= num_linsolve_global .and. iter >= 1 ) then
           itask = ZPARES_TASK_SOLVE_H
           call start_timer(prm%timer_solve_H)
        else
           itask = ZPARES_TASK_FACTO_H
           call start_timer(prm%timer_fact_H)
        end if
     end if
  end select  

end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_svd)(prm, nrow_local, rwork, X, info)
  use zpares_aux
  implicit none
  type(zpares_prm), intent(inout), target :: prm
  integer, intent(in) :: nrow_local
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  integer, intent(inout) :: info  
  integer :: i, j, k, L, M, Lmax, LM, LmaxM, ierr, infola
  integer, pointer :: num_basis
  MATRIX_TYPE, allocatable :: tmp_mat(:,:), H0(:,:), H1(:,:)
  MATRIX_TYPE :: dummy_U(1,1), dummy_VT(1,1)
  REAL_TYPE :: delta
  REAL_TYPE, allocatable :: sigma(:)
  MATRIX_TYPE, pointer :: Mu(:,:), projAB(:,:)
  
  L = prm%L
  M = prm%M
  Lmax = prm%Lmax
  LM = prm%L*prm%M
  LmaxM = prm%Lmax*prm%M
  num_basis => prm%num_basis
  delta = prm%delta

  ! pack columns of S
  if ( L /= Lmax ) then
     j = L + 1
     do i = 2, M
        do k = (i-1)*Lmax+1, (i-1)*Lmax+L
           X(:,j) = X(:,k)
           j = j + 1
        end do
     end do
  end if
  
  allocate(sigma(LM))
  select case(prm%extract)
  case(ZPARES_EXTRACT_RR)
     
     allocate(tmp_mat(LM,LM))
     call start_timer(prm%timer_orth)
     do i = 1, prm%n_orth
        call orth_SVD(nrow_local, LM, X, rwork, Lmax, tmp_mat, LM, i==1, prm%low_comm, infola, ierr) 
        ! if ( infola /= 0 ) then
        !    info = ZPARES_INFO_LAPACK_ERROR
        !    return
        ! end if
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
     end do
     call stop_timer(prm%timer_orth)
     if ( prm%n_orth > 0 ) then
        call start_timer(prm%timer_serial_SVD)
        call serial_SVD('R', LM, LM, tmp_mat, LM, delta, sigma, dummy_U, 1, dummy_VT, 1, num_basis, infola)
        ! if ( infola /= 0 ) then
        !    info = ZPARES_INFO_LAPACK_ERROR
        !    return
        ! end if
        prm%sig_val(1:LM) = sigma(1:LM)
        call stop_timer(prm%timer_serial_SVD)
        call basis_rotation('C', nrow_local, num_basis, LM, tmp_mat, LM, rwork, nrow_local*Lmax/(LM), X)
     else
        num_basis = LM
     end if
     deallocate(tmp_mat)
     
  case(ZPARES_EXTRACT_EM)
     
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
     
     call zpares_rci_sub_get_projAB(prm, projAB)        
     allocate(H0(LM, LM), H1(LM, LM))
     call block_Hankel(Lmax, L, M, 0, Mu, H0)
     ! TODO SVD is not needed if it is after LADD
     call start_timer(prm%timer_serial_SVD)
     call serial_SVD('B', LM, LM, H0, LM, delta, sigma, projAB(:,1:LmaxM), LmaxM &
          , projAB(:,LmaxM+1:2*LmaxM), LmaxM, num_basis, infola)
     ! if ( infola /= 0 ) then
     !    info = ZPARES_INFO_LAPACK_ERROR
     !    return
     ! end if
     prm%sig_val(1:LM) = sigma(1:LM)
     call stop_timer(prm%timer_serial_SVD)
     do i = 1, LM
        projAB(1:num_basis,LmaxM+i) = projAB(1:num_basis,LmaxM+i) / sqrt(sigma(1:num_basis))
     end do
     call block_Hankel(Lmax, L, M, 1, Mu, H1)
     call MACRO_ADD_PRFX(GEMM)('N', 'C', LM, num_basis, LM &
          , ONE_M, H1, LM, projAB(1,LmaxM+1), LmaxM, ZERO_M, H0, LM)
     call MACRO_ADD_PRFX(GEMM)('C', 'N', num_basis, num_basis, LM &
          , ONE_M, projAB(1,1), LmaxM, H0, LM, ZERO_M, H1, LM)
     do i = 1, num_basis
        projAB(1:num_basis,i) = H1(1:num_basis,i) / sqrt(sigma(1:num_basis))
     end do
     deallocate(H0, H1)
     
  end select
       
  deallocate(sigma)
end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_Rayleigh_Ritz) &
     (prm, nrow_local, rwork, eigval, X, info)
  use zpares_aux
  implicit none
#ifdef REALMAT
  logical, parameter :: real_matrix = .true.
#else
  logical, parameter :: real_matrix = .false.
#endif
  type(zpares_prm), intent(inout), target :: prm
  integer, intent(in) :: nrow_local
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  COMPLEX_TYPE, intent(inout) :: eigval(*)
  integer, intent(inout) :: info
  integer, pointer :: itask, state, x_offset, x_ncol, num_basis, ws, xs, nc, iter
  integer :: i ,j, Lmax, LmaxM, high_comm_rank, high_comm_size, LDQAQ, LDQBQ, ierr, infola
  ! integer, allocatable :: real_unsym_map(:), indices(:)
  ! REAL_TYPE :: res_max, tol
  ! REAL_TYPE, allocatable :: AXnorm(:), BXnorm(:), tmprvec(:)
  MATRIX_TYPE, allocatable :: tmp_mat(:,:), QAQ(:,:), QBQ(:,:)
  MATRIX_TYPE, pointer :: projAB(:,:)
  logical :: real_unsym

  real_unsym = real_matrix .and. ( .not. prm%Hermitian)

  itask => prm%itask
  state => prm%state
  x_offset => prm%x_offset
  x_ncol => prm%x_ncol
  num_basis => prm%num_basis
  ws => prm%ws
  xs => prm%xs
  nc => prm%nc
  iter => prm%iter
  Lmax = prm%Lmax
  LmaxM = prm%Lmax*prm%M
  call zpares_rci_sub_get_projAB(prm, projAB)
  call get_rank_and_size(prm%high_comm, high_comm_rank, high_comm_size)

  select case(state)
  case(ZPARES_STATE_INIT_REDUCED_EIG)
          
     x_offset = 1
     x_ncol = 0
     projAB(1:num_basis, 1:num_basis) = ZERO_M
     projAB(1:num_basis, LmaxM+1:LmaxM+num_basis) = ZERO_M
     state = ZPARES_STATE_FORM_REDUCED_EIG

  case(ZPARES_STATE_FORM_REDUCED_EIG)

     select case(itask)
     case(ZPARES_TASK_NONE)
        x_offset = x_offset + x_ncol
        if ( x_offset > num_basis ) then
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_SOLVE_REDUCED_EIG
           return
        end if
        x_ncol = min(num_basis - x_offset + 1, Lmax)
        call para_range(x_ncol, high_comm_size, high_comm_rank, ws, nc)
        if ( nc /= 0 ) then
           xs = x_offset + ws - 1
           itask = ZPARES_TASK_MULT_A
           call start_timer(prm%timer_mult_A)
        else
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_FORM_REDUCED_EIG
        end if
     case(ZPARES_TASK_MULT_A)
        call stop_timer(prm%timer_mult_A)
        allocate(tmp_mat(num_basis,nc))
        call start_timer(prm%timer_gemm_reduce)
        call GEMM_ALLREDUCE('C', 'N', num_basis, nc, nrow_local &
             , ONE_M, X(:,1:num_basis), nrow_local, rwork(:,ws:ws+nc-1), nrow_local &
             , ZERO_M, tmp_mat, num_basis, ierr, prm%low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        call stop_timer(prm%timer_gemm_reduce)
        projAB(1:num_basis,xs:xs+nc-1)  = tmp_mat(:,:)
        deallocate(tmp_mat)
        if ( prm%standard ) then
           rwork(:,ws:ws+nc-1) = X(:,xs:xs+nc-1)
        end if
        itask = ZPARES_TASK_MULT_B
        call start_timer(prm%timer_mult_B)
     case(ZPARES_TASK_MULT_B)
        call stop_timer(prm%timer_mult_B)
        allocate(tmp_mat(num_basis,nc))
        call start_timer(prm%timer_gemm_reduce)
        call GEMM_ALLREDUCE('C', 'N', num_basis, nc, nrow_local &
             , ONE_M, X(:,1:num_basis), nrow_local, rwork(:,ws:ws+nc-1), nrow_local &
             , ZERO_M, tmp_mat, num_basis, ierr, prm%low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        projAB(1:num_basis,LmaxM+xs:LmaxM+xs+nc-1) = tmp_mat(:,:)
        call stop_timer(prm%timer_gemm_reduce)
        deallocate(tmp_mat)
        itask = ZPARES_TASK_NONE
     end select
     
  case(ZPARES_STATE_SOLVE_REDUCED_EIG)
     
     call ALLREDUCE_SUM_2D(projAB, 1, LmaxM, 2*LmaxM, ierr, prm%high_comm)
     if ( ierr /= 0 ) then
        info = ZPARES_INFO_MPI_ERROR
        return
     end if
     ! now projAB(1:num_basis,1:num_basis) = Q'*A*Q
     ! now projAB(1:num_basis,Lmax*M+1:Lmax*M+num_basis) = Q'*B*Q
     allocate(QAQ(1:num_basis,1:num_basis))
     QAQ = projAB(1:num_basis, 1:num_basis)
     allocate(QBQ(1:num_basis,1:num_basis))
     QBQ = projAB(1:num_basis, LmaxM+1:LmaxM+num_basis)
     call start_timer(prm%timer_reduced_eig)
     LDQAQ = num_basis
     LDQBQ = num_basis
     if ( prm%B_pos_def .and. prm%Hermitian ) then
     ! if ( prm%B_pos_def .and. prm%Hermitian .and. ( prm%real_matrix .and. sym_contour .or. ( .not. prm%real_matrix ) ) ) then
        call HEGV_reduced_eig(nrow_local, num_basis, QAQ, LDQAQ, QBQ, LDQBQ, eigval, infola)
     else
        call GEGV_reduced_eig(nrow_local, num_basis, QAQ, LDQAQ, QBQ, LDQBQ, eigval, infola)
     end if
     ! if ( infola /= 0 ) then
     !    info = ZPARES_INFO_LAPACK_ERROR
     !    return
     ! end if
     ! now QAQ is the (right) eigenvectors of the reduced eigenvalue problem
     call stop_timer(prm%timer_reduced_eig)
     deallocate(QBQ)
     call start_timer(prm%timer_sub_rot)
     call basis_rotation('N', nrow_local, num_basis, num_basis, QAQ, LDQAQ, rwork, nrow_local*Lmax/num_basis, X)
     call stop_timer(prm%timer_sub_rot)
     if ( prm%n_orth > 0 ) then
        if ( real_unsym ) then
           do i = 1, num_basis
              if ( aimag(eigval(i)) == ZERO_R ) then
                 prm%indi_spu(i) = sum(QAQ(:,i)**2) / sum(QAQ(:,i)**2/prm%sig_val(1:num_basis))
              else if ( aimag(eigval(i)) > ZERO_R ) then
                 prm%indi_spu(i) = sum(QAQ(:,i)**2 + QAQ(:,i+1)**2) / sum((QAQ(:,i)**2+QAQ(:,i+1)**2)/prm%sig_val(1:num_basis))
                 prm%indi_spu(i+1) = prm%indi_spu(i)
              end if
           end do
        else
           do i = 1, num_basis
              prm%indi_spu(i) = sum(abs(QAQ(:,i))**2) / sum(abs(QAQ(:,i))**2/prm%sig_val(1:num_basis))
           end do
        end if
     else
        prm%indi_spu(1:num_basis) = 1d0
     end if
     prm%indi_spu(1:num_basis) = prm%indi_spu(1:num_basis) / maxval(prm%indi_spu(1:num_basis))
     deallocate(QAQ)
     state = ZPARES_STATE_FIN_REDUCED_EIG
     
  end select
end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_Hankel_method) &
     (prm, nrow_local, rwork, left, right, eigval, X, info, set_rule)
  use zpares_aux
  implicit none
#include "set_rule.f90"
  optional set_rule
#ifdef REALMAT
  logical, parameter :: real_matrix = .true.
#else
  logical, parameter :: real_matrix = .false.
#endif
  type(zpares_prm), intent(inout), target :: prm  
  integer, intent(in) :: nrow_local  
  integer, intent(inout) :: info
  COMPLEX_TYPE, intent(in) :: left
  REAL_TYPE, intent(in) :: right
  COMPLEX_TYPE, intent(out) :: eigval(*)
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)  
  integer :: i, L, M, Lmax, num_basis, dummy1, dummy2, LM, LmaxM, infola
  MATRIX_TYPE, pointer :: projAB(:,:)
  MATRIX_TYPE, allocatable :: tmp_mat(:,:)
  REAL_TYPE :: radius
  COMPLEX_TYPE :: center, dummy3, dummy4, dummy5
  logical :: real_unsym  

  select case(prm%state)
  case(ZPARES_STATE_INIT_REDUCED_EIG)

     ! do nothing
     prm%state = ZPARES_STATE_SOLVE_REDUCED_EIG

  case(ZPARES_STATE_SOLVE_REDUCED_EIG)

     num_basis = prm%num_basis
     L = prm%L
     M = prm%M
     Lmax = prm%Lmax
     LM = prm%L*prm%M
     LmaxM = prm%Lmax*prm%M
     real_unsym = real_matrix .and. ( .not. prm%Hermitian)
     call zpares_rci_sub_get_projAB(prm, projAB)
     call start_timer(prm%timer_reduced_eig)
     call GEEV_reduced_eig(nrow_local, num_basis, projAB(:,1:num_basis), LmaxM, eigval, infola)
     ! if ( infola /= 0 ) then
     !    info = ZPARES_INFO_LAPACK_ERROR
     !    return
     ! end if
     call stop_timer(prm%timer_reduced_eig)

     allocate(tmp_mat(LM, num_basis))
     tmp_mat(1:num_basis, 1:num_basis) = projAB(1:num_basis, 1:num_basis)     
     if ( real_unsym ) then
        do i = 1, num_basis           
           if ( aimag(eigval(i)) == ZERO_R ) then
              prm%indi_spu(i) = sum(tmp_mat(:,i)**2) / sum(tmp_mat(:,i)**2/prm%sig_val(1:num_basis))
           else if ( aimag(eigval(i)) > ZERO_R ) then
              prm%indi_spu(i) = sum(tmp_mat(:,i)**2 + tmp_mat(:,i+1)**2) &
                   / sum((tmp_mat(:,i)**2+tmp_mat(:,i+1)**2)/prm%sig_val(1:num_basis))
              prm%indi_spu(i+1) = prm%indi_spu(i)
           end if
        end do
     else
        do i = 1, num_basis
           prm%indi_spu(i) = sum(abs(tmp_mat(:,i))**2) / sum(abs(tmp_mat(:,i))**2/prm%sig_val(1:num_basis))
        end do
     end if
     prm%indi_spu(1:num_basis) = prm%indi_spu(1:num_basis) / maxval(prm%indi_spu(1:num_basis))

     call MACRO_ADD_PRFX(GEMM)('C','N', LM, num_basis, num_basis &
          , ONE_M, projAB(1,LmaxM+1), LmaxM, projAB(1,1), LmaxM, ZERO_M, tmp_mat, LM)
     call start_timer(prm%timer_sub_rot)
     call basis_rotation('N', nrow_local, num_basis, LM, tmp_mat, LM, rwork, nrow_local*Lmax/(LM), X)
     call stop_timer(prm%timer_sub_rot)
     select case(prm%quad_type)
     case(ZPARES_QUAD_ELL_TRAP)
        call calc_center_radius(left, right, center, radius)
        ! if ( real_matrix ) then
        !    do i = 1, num_basis
        !       if ( aimag(eigval(i)) >= 0d0 ) then
        !          eigval(i) = center + radius*eigval(i))
        !       else
        !          eigval(i) = conjg(center + radius*eigval(i))
        !       end if
        !    end do
        ! else
           eigval(1:num_basis) = center + radius*eigval(1:num_basis)     
        ! end if
     case(ZPARES_QUAD_USER)
        if ( present(set_rule) ) then
           do i = 1, num_basis
              call set_rule(ZPARES_QU_MODE_BACK, dummy1, dummy2, left, right, dummy3, dummy4, dummy5, eigval(i))
           end do
        else
           !! error TODO: error handling
        end if
     end select
     prm%state = ZPARES_STATE_FIN_REDUCED_EIG
     deallocate(tmp_mat)

  end select

end subroutine

subroutine MACRO_INSERT_PRFX(zpares_,rci_sub_calc_res) &
     (prm, nrow_local, rwork, cwork, eigval, X, res, info)
  use zpares_aux
#ifdef REALMAT
  logical, parameter :: real_matrix = .true.
#else
  logical, parameter :: real_matrix = .false.
#endif
  type(zpares_prm), intent(inout), target :: prm
  integer, intent(in) :: nrow_local
  MATRIX_TYPE, intent(inout) :: rwork(nrow_local,*), X(nrow_local,*)
  REAL_TYPE, intent(out) :: res(*)
  COMPLEX_TYPE, intent(inout) :: cwork(nrow_local,*), eigval(*)
  integer, intent(inout) :: info
  integer, pointer :: itask, state, x_offset, x_ncol, num_basis, ws, xs, nc, iter  
  integer :: i ,j, Lmax, counter, high_comm_rank, high_comm_size, ierr
  integer, allocatable :: real_unsym_map(:), indices(:)
  REAL_TYPE :: res_max, tol
  REAL_TYPE, allocatable :: AXnorm(:), BXnorm(:), tmprvec(:)
  MATRIX_TYPE, allocatable :: tmp_mat(:,:)
  MATRIX_TYPE, pointer :: projAB(:,:)
  logical :: real_unsym
  logical, allocatable :: flags(:)

  itask => prm%itask
  state => prm%state
  x_offset => prm%x_offset
  x_ncol => prm%x_ncol
  num_basis => prm%num_basis
  ws => prm%ws
  xs => prm%xs
  nc => prm%nc
  iter => prm%iter
  Lmax = prm%Lmax
  real_unsym = real_matrix .and. ( .not. prm%Hermitian )
  call zpares_rci_sub_get_projAB(prm, projAB)
  call get_rank_and_size(prm%high_comm, high_comm_rank, high_comm_size)
  
  select case(state)
  case(ZPARES_STATE_CALC_RES0)
     
     x_offset = 1
     x_ncol = 0
     state = ZPARES_STATE_CALC_RES1
     res(1:num_basis) = ZERO_R
     projAB(1:num_basis,1:2) = ZERO_M     
     
  case(ZPARES_STATE_CALC_RES1)

     select case(itask)
     case(ZPARES_TASK_NONE)
        x_offset = x_offset + x_ncol
        if ( x_offset > num_basis ) then
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_CALC_RES2
           return
        end if
        x_ncol = min(num_basis - x_offset + 1, Lmax)
        if ( real_unsym ) then
           if ( aimag(eigval(x_offset+x_ncol-1)) > ZERO_R ) x_ncol = x_ncol - 1
           counter = 0 ! count number of unique values : real(lambda) + abs(imag(lambda))*1i
           do i = x_offset, x_offset+x_ncol-1
              if ( aimag(eigval(i)) <= ZERO_R ) counter = counter + 1
           end do
           allocate(real_unsym_map(counter))
           j = 1
           do i = 1, counter
              real_unsym_map(i) = j
              if ( aimag(eigval(x_offset+j-1)) == ZERO_R ) then                    
                 j = j + 1
              else
                 j = j + 2
              end if
           end do
        else
           counter = x_ncol
        end if
        call para_range(counter, high_comm_size, high_comm_rank, ws, nc)
        if ( real_unsym ) then
           if ( nc /= 0 ) then
              j = ws+nc-1
              ws = real_unsym_map(ws)
              nc = real_unsym_map(j) - ws + 1
              if ( aimag(eigval(x_offset + real_unsym_map(j) - 1)) /= ZERO_R ) then                    
                 nc = nc + 1
              end if
           end if
           deallocate(real_unsym_map)
        end if
        if ( nc /= 0 ) then
           xs = x_offset + ws - 1
           itask = ZPARES_TASK_MULT_A
           call start_timer(prm%timer_mult_A)
        else
           itask = ZPARES_TASK_NONE
           state = ZPARES_STATE_CALC_RES1
        end if
     case(ZPARES_TASK_MULT_A)
        call stop_timer(prm%timer_mult_A)
        ! rwork = AX
        allocate(AXnorm(nc))
        call start_timer(prm%timer_res_norm)
        call norm2_blk(rwork(:,ws:ws+nc-1), AXnorm, nrow_local, nc, ierr, prm%low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        call stop_timer(prm%timer_res_norm)
        projAB(xs:xs+nc-1,1) = AXnorm(1:nc)
        deallocate(AXnorm)
        cwork(:,ws:ws+nc-1) = rwork(:,ws:ws+nc-1)
        if ( prm%standard ) then
           rwork(:,ws:ws+nc-1) = X(:,xs:xs+nc-1)
        end if
        itask = ZPARES_TASK_MULT_B
        call start_timer(prm%timer_mult_B)
     case(ZPARES_TASK_MULT_B)   
        call stop_timer(prm%timer_mult_B)
        ! rwork = BX
        allocate(BXnorm(nc))
        call start_timer(prm%timer_res_norm)
        call norm2_blk(rwork(:,ws:ws+nc-1), BXnorm, nrow_local, nc, ierr, prm%low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        call stop_timer(prm%timer_res_norm)
        projAB(xs:xs+nc-1,2) = BXnorm(1:nc)
        if ( real_unsym ) then
           i = 1
           allocate(tmprvec(nrow_local))
           do while ( i <= nc )
              j = ws+i-1
              if ( aimag(eigval(x_offset+j-1)) == ZERO_R ) then
                 rwork(:,j) = cwork(:,j) - eigval(x_offset+j-1)*rwork(:,j)
                 i = i + 1
              else
                 tmprvec(:) = rwork(:,j)
                 rwork(:,j) = real(cwork(:,j),kind(ZERO_R)) - real(eigval(x_offset+j-1),kind(ZERO_R))*tmprvec(:) &
                      + aimag(eigval(x_offset+j-1))*rwork(:,j+1) ! real part of residual
                 rwork(:,j+1) = real(cwork(:,j+1)) - aimag(eigval(x_offset+j-1))*tmprvec(:) &
                      - real(eigval(x_offset+j-1),kind(ZERO_R))*rwork(:,j+1) ! imag part of residual                
                 i = i + 2                 
              end if              
           end do
           deallocate(tmprvec)
        else
           do i = 1, nc
              j = ws+i-1
              rwork(:,j) = cwork(:,j) - eigval(x_offset+j-1)*rwork(:,j)
           end do
        end if
        ! now rwork is the residual
        call start_timer(prm%timer_res_norm)
        call norm2_blk(rwork(:,ws:ws+nc-1), res(xs:xs+nc-1), nrow_local, nc, ierr, prm%low_comm)
        if ( ierr /= 0 ) then
           info = ZPARES_INFO_MPI_ERROR
           return
        end if
        call stop_timer(prm%timer_res_norm)
        deallocate(BXnorm)
        itask = ZPARES_TASK_NONE
     end select
     
  case(ZPARES_STATE_CALC_RES2)
     
     call ALLREDUCE_SUM_1D(res, num_basis, ierr, prm%high_comm)
     if ( ierr /= 0 ) then
        info = ZPARES_INFO_MPI_ERROR
        return
     end if
     allocate(AXnorm(num_basis), BXnorm(num_basis))
     AXnorm(1:num_basis) = projAB(1:num_basis,1)
     BXnorm(1:num_basis) = projAB(1:num_basis,2)
     call ALLREDUCE_SUM_1D(AXnorm, num_basis, ierr, prm%high_comm)
     if ( ierr /= 0 ) then
        info = ZPARES_INFO_MPI_ERROR
        return
     end if
     call ALLREDUCE_SUM_1D(BXnorm, num_basis, ierr, prm%high_comm)
     if ( ierr /= 0 ) then
        info = ZPARES_INFO_MPI_ERROR
        return
     end if
     if ( real_unsym ) then
        ! if A,B are real unsymmertic matrices and symmetric contour path with respect to the real-axis
        ! eigenvectors are stored DGEGV's way
        i = 1
        do while ( i <= num_basis )
           if ( aimag(eigval(i)) == ZERO_R ) then
              res(i) = res(i) / (AXnorm(i) + abs(eigval(i))*BXnorm(i))
              i = i + 1
           else
              res(i) = sqrt(res(i)**2 + res(i+1)**2)
              AXnorm(i) = sqrt(AXnorm(i)**2 + AXnorm(i+1)**2)
              BXnorm(i) = sqrt(BXnorm(i)**2 + BXnorm(i+1)**2)
              res(i) = res(i) / ( AXnorm(i) + abs(eigval(i)) * BXnorm(i) )
              res(i+1) = res(i)
              i = i + 2                 
           end if
        end do
     else
        res(1:num_basis) = res(1:num_basis) / (AXnorm(1:num_basis) + abs(eigval(1:num_basis))*BXnorm(1:num_basis))
     end if
     deallocate(AXnorm, BXnorm)
     state = ZPARES_STATE_CHECK_OUTER

  end select
end subroutine

#ifndef REALMAT
logical function MACRO_INSERT_PRFX(zpares_,rci_sub_get_sym_contour)(prm, left) result(ret)
  implicit none
  type(zpares_prm) :: prm
  COMPLEX_TYPE, intent(in) :: left
  if ( prm%quad_type == ZPARES_QUAD_USER ) then
     ret = prm%sym_contour
  else
     ret = (prm%quad_type == ZPARES_QUAD_ELL_TRAP) .and. (aimag(left) == ZERO_R)
  end if
end function
#endif



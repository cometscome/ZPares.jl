#include "def_macros.h"

#ifndef REALMAT
#ifdef SINGLE
  subroutine c_MUMPS_wrap(mumps_par)
    implicit none
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
    call CMUMPS(mumps_par)
  end subroutine c_MUMPS_wrap
#else
  subroutine z_MUMPS_wrap(mumps_par)
    implicit none
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
    call ZMUMPS(mumps_par)
  end subroutine z_MUMPS_wrap
#endif
#endif


!> @brief Sparse driver routine for general unsymmetric matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] rowptrB : Row pointers of matrix B
!! @param[in] colindB : Column indices of matrix B
!! @param[in] valB : Values of matrix B
!! @param[in] left left edge of ellipse
!! @param[in] right real part of right edge of ellipse
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpsgegv) &
  (prm, mat_size, rowptrA, colindA, valA, rowptrB, colindB, valB, left, right, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
#include "set_rule.f90"
    optional set_rule
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    REAL_TYPE, intent(in) :: right
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(out) :: eigval(*)
    COMPLEX_TYPE, intent(in) :: left
    MATRIX_TYPE, intent(in) :: valA(*), valB(*)
    MATRIX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M
    COMPLEX_TYPE :: z
    MATRIX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0

    integer :: i
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 0
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(rwork(mat_size, Lmax), cwork(mat_size*Lmax))

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zB_A(mat_size, rowptrA, colindA, rowptrB, colindB)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.

    call idx_COO_zB_A(mat_size, rowptrA, colindA, rowptrB, colindB, mumps_par%IRN, mumps_par%JCN)    
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcigegv) &
       (prm, mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zB_A(mat_size, z, rowptrA, colindA, valA, rowptrB, colindB, valB, mumps_par%A)
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_MULT_A)

          call matvec_csr(rowptrA(mat_size+1)-1,mat_size, rowptrA, colindA, valA, prm%nc, &
               X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       case(ZPARES_TASK_MULT_B)
          
          call matvec_csr(rowptrB(mat_size+1)-1,mat_size, rowptrB, colindB, valB, prm%nc, &
               X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
    
  end subroutine


#ifdef REALMAT
!> @brief Sparse driver routine for real symmmetric matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] rowptrB : Row pointers of matrix B
!! @param[in] colindB : Column indices of matrix B
!! @param[in] valB : Values of matrix B
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpssygv) &
  (prm, mat_size, rowptrA, colindA, valA, rowptrB, colindB, valB, emin, emax, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#include "set_rule.f90"
    optional set_rule
#ifdef MPI
    include 'mpif.h'
#endif
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*)
    REAL_TYPE, intent(in) :: valA(*), valB(*)
    REAL_TYPE, intent(out) :: X(mat_size,*)
    REAL_TYPE, intent(out) :: eigval(*)
    
    integer :: L, Lmax, M
    COMPLEX_TYPE :: z
    integer, allocatable :: diagptrA(:), diagptrB(:)
    REAL_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 2
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(diagptrA(mat_size), diagptrB(mat_size) &
         , rwork(mat_size, Lmax), cwork(mat_size*Lmax))
    
    call get_diagptr(mat_size, rowptrA, colindA, diagptrA)
    call get_diagptr(mat_size, rowptrB, colindB, diagptrB)

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zB_A(mat_size, rowptrA, colindA, rowptrB, colindB)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.

    call idx_COO_zB_A(mat_size, rowptrA, colindA, rowptrB, colindB, mumps_par%IRN, mumps_par%JCN)
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcisygv) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zB_A(mat_size, z, rowptrA, colindA, valA, rowptrB, colindB, valB, mumps_par%A)
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_MULT_A)

          call matvec_csr_hermite(mat_size, rowptrA, colindA, valA, diagptrA, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       case(ZPARES_TASK_MULT_B)
          
          call matvec_csr_hermite(mat_size, rowptrB, colindB, valB, diagptrB, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(diagptrA, diagptrB, rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
  end subroutine


#else


!> @brief Sparse driver routine for complex Hermitian matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] rowptrB : Row pointers of matrix B
!! @param[in] colindB : Column indices of matrix B
!! @param[in] valB : Values of matrix B
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpshegv) &
  (prm, mat_size, rowptrA, colindA, valA, rowptrB, colindB, valB, emin, emax, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#include "set_rule.f90"
    optional set_rule
#ifdef MPI
    include 'mpif.h'
#endif
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(in) :: valA(*), valB(*)
    COMPLEX_TYPE, intent(out) :: X(mat_size,*)
    REAL_TYPE, intent(out) :: eigval(*)
    
    integer :: L, Lmax, M
    COMPLEX_TYPE :: z
    integer, allocatable :: diagptrA(:), diagptrB(:)
    COMPLEX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0

    integer :: i
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 0
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(diagptrA(mat_size), diagptrB(mat_size) &
         , rwork(mat_size, Lmax), cwork(mat_size*Lmax))

    call get_diagptr(mat_size, rowptrA, colindA, diagptrA)
    call get_diagptr(mat_size, rowptrB, colindB, diagptrB)

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zB_A_h2f(mat_size, rowptrA, colindA, rowptrB, colindB)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.

    call idx_COO_zB_A_h2f(mat_size, rowptrA, colindA, rowptrB, colindB, mumps_par%IRN, mumps_par%JCN)          
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )       
       call MACRO_INSERT_PRFX(zpares_,rcihegv) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zB_A_h2f(mat_size, z, rowptrA, colindA, valA, rowptrB, colindB, valB, mumps_par%A)          
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE_H)

          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%RHS(:) = conjg(mumps_par%RHS(:))
          mumps_par%ICNTL(9) = 0
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)
          mumps_par%RHS(:) = conjg(mumps_par%RHS(:))
          mumps_par%ICNTL(9) = 1

       case(ZPARES_TASK_MULT_A)

          call matvec_csr_hermite(mat_size, rowptrA, colindA, valA, diagptrA, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       case(ZPARES_TASK_MULT_B)
          
          call matvec_csr_hermite(mat_size, rowptrB, colindB, valB, diagptrB, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))
          
       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(diagptrA, diagptrB, rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
  end subroutine
#endif


!> @brief Sparse driver routine for general unsymmetric matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] left left edge of ellipse
!! @param[in] right real part of right edge of ellipse
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpsgeev) &
  (prm, mat_size, rowptrA, colindA, valA, left, right, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#include "set_rule.f90"
    optional set_rule
#ifdef MPI
    include 'mpif.h'
#endif
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*)
    REAL_TYPE, intent(in) :: right
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(in) :: left
    COMPLEX_TYPE, intent(out) :: eigval(*)
    MATRIX_TYPE, intent(in) :: valA(*)
    MATRIX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M
    integer, allocatable :: diagptrA(:)
    COMPLEX_TYPE :: z
    MATRIX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0

    integer :: i
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 0
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(diagptrA(mat_size), rwork(mat_size, Lmax), cwork(mat_size*Lmax))
    
    call get_diagptr(mat_size, rowptrA, colindA, diagptrA)

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zI_A(mat_size, rowptrA, diagptrA)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.

    call idx_COO_zI_A(mat_size, rowptrA, colindA, diagptrA, mumps_par%IRN, mumps_par%JCN)
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )       
       call MACRO_INSERT_PRFX(zpares_,rcigeev) &
       (prm, mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zI_A(mat_size, z, rowptrA, colindA, valA, diagptrA, mumps_par%A)
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_MULT_A)

          call matvec_csr(rowptrA(mat_size+1)-1,mat_size, rowptrA, colindA, valA, prm%nc, &
               X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(diagptrA, rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
    
  end subroutine


#ifdef REALMAT
!> @brief Sparse driver routine for real symmetric matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpssyev) &
  (prm, mat_size, rowptrA, colindA, valA, emin, emax, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#include "set_rule.f90"
    optional set_rule
#ifdef MPI
    include 'mpif.h'
#endif
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*)
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*)
    REAL_TYPE, intent(in) :: valA(*)
    REAL_TYPE, intent(out) :: X(mat_size,*)
    REAL_TYPE, intent(out) :: eigval(*)
    
    integer :: L, Lmax, M
    COMPLEX_TYPE :: z
    integer, allocatable :: diagptrA(:)
    REAL_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 2
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(diagptrA(mat_size), rwork(mat_size, Lmax), cwork(mat_size*Lmax))
    
    call get_diagptr(mat_size, rowptrA, colindA, diagptrA)

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zI_A(mat_size, rowptrA, diagptrA)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.
    
    call idx_COO_zI_A(mat_size, rowptrA, colindA, diagptrA, mumps_par%IRN, mumps_par%JCN)

    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcisyev) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zI_A(mat_size, z, rowptrA, colindA, valA, diagptrA, mumps_par%A)
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_MULT_A)

          call matvec_csr_hermite(mat_size, rowptrA, colindA, valA, diagptrA, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(diagptrA, rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
  end subroutine


#else


!> @brief Sparse driver routine for complex Hermitian matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] rowptrA : Row pointers of matrix A
!! @param[in] colindA : Column indices of matrix A
!! @param[in] valA : Values of matrix A
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,mpsheev) &
  (prm, mat_size, rowptrA, colindA, valA, emin, emax, num_ev, eigval, X, res, info, set_rule)
    use zpares
    implicit none
#include "set_rule.f90"
    optional set_rule
#ifdef MPI
    include 'mpif.h'
#endif
#ifdef SINGLE
    include 'cmumps_struc.h'
    type(cmumps_struc) :: mumps_par
#else
    include 'zmumps_struc.h'
    type(zmumps_struc) :: mumps_par
#endif
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    integer, intent(in) :: rowptrA(*), colindA(*)
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(in) :: valA(*)
    COMPLEX_TYPE, intent(out) :: X(mat_size,*)
    REAL_TYPE, intent(out) :: eigval(*)
    
    integer :: L, Lmax, M
    COMPLEX_TYPE :: z
    integer, allocatable :: diagptrA(:)
    COMPLEX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable, target :: cwork(:)
    double precision :: t0, t = 0d0

    integer :: i
  
    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    !-- Control of parallelism.
#ifdef MPI
    mumps_par%COMM = MPI_COMM_SELF  ! Communicator for MUMPS.
#endif
    mumps_par%PAR  = 1 ! 1: Host is involved in factorization/solve phases.    
    !-- Matrix type.
    mumps_par%SYM  = 0
    !-- Initialization of MUMPS.
    mumps_par%JOB = -1
    call MUMPS_wrap(mumps_par)
    
    allocate(diagptrA(mat_size), rwork(mat_size, Lmax), cwork(mat_size*Lmax))

    call get_diagptr(mat_size, rowptrA, colindA, diagptrA)

    call set_mumps_icntl(mumps_par%ICNTL)
    mumps_par%N = mat_size
    mumps_par%NZ = count_nnz_zI_A_h2f(mat_size, rowptrA, diagptrA)
    mumps_par%LRHS = mumps_par%N

    allocate ( mumps_par%A(mumps_par%NZ), &  ! Nonzero elements of A.
         mumps_par%IRN(mumps_par%NZ), &  ! Row index of A.
         mumps_par%JCN(mumps_par%NZ) )   ! Column index of A.

    call idx_COO_zI_A_h2f(mat_size, rowptrA, colindA, diagptrA, mumps_par%IRN, mumps_par%JCN)          
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )       
       call MACRO_INSERT_PRFX(zpares_,rciheev) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          call val_COO_zI_A_h2f(mat_size, z, rowptrA, colindA, valA, diagptrA, mumps_par%A)          
          mumps_par%JOB = 4
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE)
          
          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)

       case(ZPARES_TASK_SOLVE_H)

          mumps_par%RHS => cwork(mat_size*(prm%ws-1)+1:mat_size*(prm%ws+prm%nc-1))
          mumps_par%NRHS = prm%nc
          mumps_par%RHS(:) = conjg(mumps_par%RHS(:))
          mumps_par%ICNTL(9) = 0
          mumps_par%JOB = 3
          call MUMPS_wrap(mumps_par)
          mumps_par%RHS(:) = conjg(mumps_par%RHS(:))
          mumps_par%ICNTL(9) = 1

       case(ZPARES_TASK_MULT_A)

          call matvec_csr_hermite(mat_size, rowptrA, colindA, valA, diagptrA, prm%nc &
               , X(:,prm%xs:prm%xs+prm%nc-1), rwork(:,prm%ws:prm%ws+prm%nc-1))

       end select
    end do

    mumps_par%JOB = -2
    call MUMPS_wrap(mumps_par)
    deallocate(diagptrA, rwork, cwork, mumps_par%A, mumps_par%IRN, mumps_par%JCN)
  end subroutine
#endif

  subroutine MACRO_ADD_PRFX(_val_COO_zB_A) &
  (mat_size, z, rowptrA, colindA, valA, rowptrB, colindB, valB, valout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    MATRIX_TYPE, intent(in) :: valA(*), valB(*)
    COMPLEX_TYPE, intent(in) :: z
    COMPLEX_TYPE, intent(out) :: valout(*)
    
    integer :: i, j, jA, jB, posA, posB, pos_zB_A, posAend, posBend
    
    pos_zB_A = 1
    do i = 1, mat_size
       posA = rowptrA(i)
       posB = rowptrB(i)
       posAend = rowptrA(i+1) - 1
       posBend = rowptrB(i+1) - 1
       do
          if ( posA <= posAend ) then
             jA = colindA(posA)
          else
             jA = mat_size + 1
          end if          
          if ( posB <= posBend ) then
             jB = colindB(posB)
          else
             jB = mat_size + 1
          end if

          if ( jA == jB ) then
             valout(pos_zB_A) = z*valB(posB) - valA(posA)
             posA = posA + 1
             posB = posB + 1
          else if ( jB < jA ) then
             valout(pos_zB_A) = z*valB(posB)
             posB = posB + 1
          else if ( jA < jB ) then
             valout(pos_zB_A) = -valA(posA)
             posA = posA + 1
          end if
          pos_zB_A = pos_zB_A + 1
                    
          if ( posA > posAend .and. posB > posBend ) then
             exit
          end if          
       end do
    end do
  end subroutine


#ifndef REALMAT
  subroutine MACRO_ADD_PRFX(_val_COO_zB_A_h2f) &
  (mat_size, z, rowptrA, colindA, valA, rowptrB, colindB, valB, valout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    COMPLEX_TYPE, intent(in) :: valA(*), valB(*)
    COMPLEX_TYPE, intent(in) :: z
    COMPLEX_TYPE, intent(out) :: valout(*)
    
    integer :: i, j, jA, jB, posA, posB, pos_zB_A, posAend, posBend
    
    pos_zB_A = 1
    do i = 1, mat_size
       posA = rowptrA(i)
       posB = rowptrB(i)
       posAend = rowptrA(i+1) - 1
       posBend = rowptrB(i+1) - 1
       do
          if ( posA <= posAend ) then
             jA = colindA(posA)
          else
             jA = mat_size + 1
          end if          
          if ( posB <= posBend ) then
             jB = colindB(posB)
          else
             jB = mat_size + 1
          end if

          if ( jA == jB ) then
             valout(pos_zB_A) = z*valB(posB) - valA(posA)
             if ( jA /= i ) then
                pos_zB_A = pos_zB_A + 1
                valout(pos_zB_A) = z*conjg(valB(posB)) - conjg(valA(posA))
             end if
             posA = posA + 1
             posB = posB + 1
          else if ( jB < jA ) then
             valout(pos_zB_A) = z*valB(posB)
             if ( jB /= i ) then
                pos_zB_A = pos_zB_A + 1
                valout(pos_zB_A) = z*conjg(valB(posB))
             end if
             posB = posB + 1
          else if ( jA < jB ) then
             valout(pos_zB_A) = -valA(posA)
             if ( jA /= i ) then
                pos_zB_A = pos_zB_A + 1
                valout(pos_zB_A) = -conjg(valA(posA))
             end if
             posA = posA + 1
          end if
          pos_zB_A = pos_zB_A + 1

          if ( posA > posAend .and. posB > posBend ) then
             exit
          end if          
       end do
    end do
  end subroutine
#endif


  subroutine MACRO_ADD_PRFX(_val_COO_zI_A) &
  (mat_size, z, rowptrA, colindA, valA, diagptrA, valout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), diagptrA(*)
    MATRIX_TYPE, intent(in) :: valA(*)
    COMPLEX_TYPE, intent(in) :: z
    COMPLEX_TYPE, intent(out) :: valout(*)
    
    integer :: i, j, pos
    
    pos = 1
    do i = 1, mat_size
       do j = rowptrA(i), diagptrA(i) - 1
          valout(pos) = -valA(j)
          pos = pos + 1
       end do
       do j = diagptrA(i) + 1, rowptrA(i+1) - 1
          valout(pos) = -valA(j)
          pos = pos + 1
       end do
       if ( diagptrA(i) /= rowptrA(i+1) ) then ! there is diagonal entry
          valout(pos) = z - valA(diagptrA(i))
       else
          valout(pos) = z
       end if
       pos = pos + 1          
    end do
  end subroutine


#ifndef REALMAT
  subroutine MACRO_ADD_PRFX(_val_COO_zI_A_h2f) &
  (mat_size, z, rowptrA, colindA, valA, diagptrA, valout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), diagptrA(*)
    MATRIX_TYPE, intent(in) :: valA(*)
    COMPLEX_TYPE, intent(in) :: z
    COMPLEX_TYPE, intent(out) :: valout(*)
    
    integer :: i, j, pos
    
    pos = 1
    do i = 1, mat_size
       do j = rowptrA(i), diagptrA(i) - 1
          valout(pos) = -valA(j)
          pos = pos + 1
          valout(pos) = -conjg(valA(j))
          pos = pos + 1
       end do
       do j = diagptrA(i) + 1, rowptrA(i+1) - 1
          valout(pos) = -valA(j)
          pos = pos + 1
          valout(pos) = -conjg(valA(j))
          pos = pos + 1
       end do
       if ( diagptrA(i) /= rowptrA(i+1) ) then ! there is diagonal entry
          valout(pos) = z - valA(diagptrA(i))
       else
          valout(pos) = z
       end if
       pos = pos + 1          
    end do
  end subroutine
#endif


  subroutine MACRO_ADD_PRFX(_matvec_csr) &
  (nnz, mat_size, rowptr, colind, val, nvec, X, Y)
    implicit none
    integer, intent(in) :: nnz, mat_size, rowptr(*), colind(*), nvec
    MATRIX_TYPE, intent(in) :: val(*), X(mat_size,*)
    MATRIX_TYPE, intent(out) :: Y(mat_size,*)
    
    integer :: i, j, k

    do k = 1, nvec
       Y(:,k) = ZERO_M
    end do

    !$omp parallel do default(shared) private(i,j,k)
    do i = 1, mat_size
       do j = rowptr(i), rowptr(i+1) - 1
          do k = 1, nvec
             Y(i,k) = Y(i,k) + val(j)*X(colind(j),k)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine


  subroutine MACRO_ADD_PRFX(_matvec_csr_hermite) &
  (mat_size, rowptr, colind, val, diagptr, nvec, X, Y)
    implicit none
    integer, intent(in) :: mat_size, rowptr(*), colind(*), diagptr(*), nvec
    MATRIX_TYPE, intent(in) :: val(*), X(mat_size,*)
    MATRIX_TYPE, intent(out) :: Y(mat_size,*)
    
    integer :: i, j, k

    do k = 1, nvec
       Y(:,k) = ZERO_M
    end do
    
    !$omp parallel do default(shared) private(i,j,k)
    do i = 1, mat_size
       do j = rowptr(i), rowptr(i+1) - 1
          do k = 1, nvec
             Y(i,k) = Y(i,k) + val(j)*X(colind(j),k)
          end do
       end do
    end do
    !$omp end parallel do

    !$omp parallel do default(shared) private(i,j,k)
    do k = 1, nvec
       do i = 1, mat_size
          do j = rowptr(i), diagptr(i) - 1
#ifdef REALMAT
             Y(colind(j),k) = Y(colind(j),k) + val(j)*X(i,k)
#else
             Y(colind(j),k) = Y(colind(j),k) + conjg(val(j))*X(i,k)
#endif
          end do
       
          do j = diagptr(i) + 1, rowptr(i+1) - 1
#ifdef REALMAT
             Y(colind(j),k) = Y(colind(j),k) + val(j)*X(i,k)
#else
             Y(colind(j),k) = Y(colind(j),k) + conjg(val(j))*X(i,k)
#endif
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine

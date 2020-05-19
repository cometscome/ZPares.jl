#include "def_macros.h"


!> @brief Dense driver routine for general unsymmetric matrix
!!
!! Solve generalized eigenvalue problem \f$A x = \lambda B x\f$
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDA Leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB Leading dimension of B
!! @param[in] left left edge of ellipse
!! @param[in] right real part of right edge of ellipse
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnsgegv) &
     (prm, mat_size, A, LDA, B, LDB, left, right, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA, LDB
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    REAL_TYPE, intent(in) :: right
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(in) :: left
    COMPLEX_TYPE, intent(out) :: eigval(*)
    MATRIX_TYPE, intent(in) :: A(LDA,*), B(LDB,*)
    MATRIX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, infola
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z
    MATRIX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable :: cwork(:,:), C(:,:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax), &
         C(mat_size, mat_size), ipiv(mat_size))

    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcigegv) &
       (prm, mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          C = z*B(1:mat_size,1:mat_size) - A(1:mat_size,1:mat_size)          
          call MACRO_ADD_C_PRFX(GETRF) &
          (mat_size, mat_size, C, mat_size, ipiv, infola)              

       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(GETRS) &
          ('N', mat_size, prm%nc, C, mat_size, ipiv &
               , cwork(1,prm%ws), mat_size, infola)
          
       case(ZPARES_TASK_MULT_A)

          call MACRO_ADD_PRFX(GEMM) &
          ('N', 'N', mat_size, prm%nc, mat_size &
               , ONE_M, A, LDA, X(1,prm%xs), mat_size &
               , ZERO_M, rwork(1,prm%ws), mat_size)

       case(ZPARES_TASK_MULT_B)

          call MACRO_ADD_PRFX(GEMM) &
          ('N', 'N', mat_size, prm%nc, mat_size &
               , ONE_M, B, LDB, X(1,prm%xs), mat_size &
               , ZERO_M, rwork(1,prm%ws), mat_size)
          
       end select
    end do

    deallocate(rwork, cwork, C, ipiv)
  end subroutine


#ifdef REALMAT
!> @brief Dense driver routine for real symmetric matrix 
!! @param[in,out] prm optional parameters
!! @param[in] UPLO = 'U': Upper triangle part of A and B are stored; = 'L': Lower triangle part of A and B are stored.
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDA Leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB Leading dimension of B
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnssygv) &
  (prm, UPLO, mat_size, A, LDA, B, LDB, emin, emax, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA, LDB
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    character, intent(in) :: UPLO
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*), eigval(*)
    REAL_TYPE, intent(in) :: A(LDA,*), B(LDB,*)
    REAL_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, i, j, infola, lwork
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z, optlwork
    REAL_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable :: cwork(:,:), C(:,:), workla(:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M    
    
    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax) &
         , C(mat_size, mat_size), ipiv(mat_size))
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcisygv) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          if ( UPLO == 'L' .or. UPLO == 'l' ) then
             do i = 1, mat_size
                do j = 1, i
                   C(i,j) = z*B(i,j) - A(i,j)
                end do
             end do
          else
             do i = 1, mat_size
                do j = 1, i
                   C(j,i) = z*B(j,i) - A(j,i)
                end do
             end do
          end if          
          call MACRO_ADD_C_PRFX(SYTRF) &
          (UPLO, mat_size, C, mat_size, ipiv, optlwork, -1, infola)
          lwork = int(optlwork)
          allocate(workla(lwork))
          call MACRO_ADD_C_PRFX(SYTRF) &
          (UPLO, mat_size, C, mat_size, ipiv, workla, lwork, infola)
          deallocate(workla)

       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(SYTRS) &
          (UPLO, mat_size, prm%nc, C, mat_size, ipiv &
               , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_MULT_A)
          
          call MACRO_ADD_PRFX(SYMM) &
          ('L', UPLO, mat_size, prm%nc &
                  , ONE_R, A, LDA, X(1,prm%xs), mat_size &
                  , ZERO_R, rwork(1,prm%ws), mat_size)

       case(ZPARES_TASK_MULT_B)
          
          call MACRO_ADD_PRFX(SYMM) &
          ('L', UPLO, mat_size, prm%nc &
          , ONE_R, B, LDB, X(1,prm%xs), mat_size &
          , ZERO_R, rwork(1,prm%ws), mat_size)
          
       end select
    end do

    deallocate(rwork, cwork, C, ipiv)
  end subroutine


#else


!> @brief Dense driver routine for complex Hermitian matrix
!! @param[in,out] prm optional parameters
!! @param[in] UPLO = 'U': Upper triangle part of A and B are stored; = 'L': Lower triangle part of A and B are stored.
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDA Leading dimension of A
!! @param[in] B matrix B
!! @param[in] LDB Leading dimension of B
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnshegv) &
  (prm, UPLO, mat_size, A, LDA, B, LDB, emin, emax, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA, LDB
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    character, intent(in) :: UPLO
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*), eigval(*)
    COMPLEX_TYPE, intent(in) :: A(LDA,*), B(LDB,*)
    COMPLEX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, i, j, infola
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z
    COMPLEX_TYPE, allocatable :: rwork(:,:), cwork(:,:), C(:,:), tmpvec(:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M
    
    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax) &
         , C(mat_size, mat_size), ipiv(mat_size), tmpvec(mat_size))
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcihegv) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
                    

          if ( UPLO == 'L' .or. UPLO == 'l' ) then

             do i = 1, mat_size
                do j = 1, i-1
                   C(i,j) = z*B(i,j) - A(i,j)
                   C(j,i) = z*conjg(B(i,j)) - conjg(A(i,j))
                end do
                C(i,i) = z*B(i,i) - A(i,i)
             end do
             
          elseif ( UPLO == 'U' .or. UPLO == 'u' ) then

             do i = 1, mat_size
                do j = 1, i-1
                   C(j,i) = z*B(j,i) - A(j,i)
                   C(i,j) = z*conjg(B(j,i)) - conjg(A(j,i))
                end do
                C(i,i) = z*B(i,i) - A(i,i)
             end do
             
          end if
          call MACRO_ADD_C_PRFX(GETRF) &
          (mat_size, mat_size, C, mat_size, ipiv, infola)


       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(GETRS) &
          ('N', mat_size, prm%nc, C, mat_size, ipiv &
          , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_SOLVE_H)
          
          call MACRO_ADD_C_PRFX(GETRS) &
          ('C', mat_size, prm%nc, C, mat_size, ipiv &
          , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_MULT_A)
          
          call MACRO_ADD_PRFX(HEMM) &
          ('L', UPLO, mat_size, prm%nc &
          , ONE_C, A, LDA, X(1,prm%xs), mat_size &
          , ZERO_C, rwork(1,prm%ws), mat_size)

       case(ZPARES_TASK_MULT_B)

          call MACRO_ADD_PRFX(HEMM) &
          ('L', UPLO, mat_size, prm%nc &
          , ONE_C, B, LDB, X(1,prm%xs), mat_size &
          , ZERO_C, rwork(1,prm%ws), mat_size)

       end select

    end do
    
    deallocate(rwork, cwork, C, ipiv, tmpvec)
  end subroutine
#endif



!> @brief Dense driver routine for general unsymmetric matrix
!! @param[in,out] prm optional parameters
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDB Leading dimension of A
!! @param[in] left left edge of ellipse
!! @param[in] right real part of right edge of ellipse
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnsgeev) &
  (prm, mat_size, A, LDA, left, right, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    REAL_TYPE, intent(in) :: right
    REAL_TYPE, intent(out) :: res(*)
    COMPLEX_TYPE, intent(in) :: left
    COMPLEX_TYPE, intent(out) :: eigval(*)
    MATRIX_TYPE, intent(in) :: A(LDA,*)
    MATRIX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, infola, i
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z
    MATRIX_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable :: cwork(:,:), C(:,:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax), &
         C(mat_size, mat_size), ipiv(mat_size))

    do while ( prm%itask /= ZPARES_TASK_FINISH )       
       call MACRO_INSERT_PRFX(zpares_,rcigeev) &
       (prm, mat_size, z, rwork, cwork, left, right, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          C = -A(1:mat_size,1:mat_size)
          do i = 1, mat_size
             C(i,i) = C(i,i) + z
          end do
          call MACRO_ADD_C_PRFX(GETRF) &
          (mat_size, mat_size, C, mat_size, ipiv, infola)              

       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(GETRS) &
          ('N', mat_size, prm%nc, C, mat_size, ipiv &
               , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_MULT_A)
          
          call MACRO_ADD_PRFX(GEMM) &
          ('N', 'N', mat_size, prm%nc, mat_size &
               , ONE_M, A, LDA, X(1,prm%xs), mat_size &
               , ZERO_M, rwork(1,prm%ws), mat_size)

       end select
    end do

    deallocate(rwork, cwork, C, ipiv)
  end subroutine


#ifdef REALMAT
!> @brief Dense driver routine for real symmetric matrix 
!! @param[in,out] prm optional parameters
!! @param[in] UPLO = 'U': Upper triangle part of A and B are stored; = 'L': Lower triangle part of A and B are stored.
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDA Leading dimension of A
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnssyev) &
  (prm, UPLO, mat_size, A, LDA, emin, emax, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    character, intent(in) :: UPLO
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*), eigval(*)
    REAL_TYPE, intent(in) :: A(LDA,*)
    REAL_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, i, j, infola, lwork
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z, optlwork
    REAL_TYPE, allocatable :: rwork(:,:)
    COMPLEX_TYPE, allocatable :: cwork(:,:), C(:,:), workla(:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M    
    
    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax) &
         , C(mat_size, mat_size), ipiv(mat_size))
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rcisyev) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          
          if ( UPLO == 'L' .or. UPLO == 'l' ) then
             do i = 1, mat_size
                do j = 1, i
                   C(i,j) = -A(i,j)
                end do
             end do
          else
             do i = 1, mat_size
                do j = 1, i
                   C(j,i) = -A(j,i)
                end do
             end do
          end if
          do i = 1, mat_size
             C(i,i) = C(i,i) + z
          end do
          call MACRO_ADD_C_PRFX(SYTRF) &
          (UPLO, mat_size, C, mat_size, ipiv, optlwork, -1, infola)
          lwork = int(optlwork)
          allocate(workla(lwork))
          call MACRO_ADD_C_PRFX(SYTRF) &
          (UPLO, mat_size, C, mat_size, ipiv, workla, lwork, infola)
          deallocate(workla)

       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(SYTRS) &
          (UPLO, mat_size, prm%nc, C, mat_size, ipiv &
          , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_MULT_A)
          
          call MACRO_ADD_PRFX(SYMM) &
          ('L', UPLO, mat_size, prm%nc &
          , ONE_R, A, LDA, X(1,prm%xs), mat_size &
          , ZERO_R, rwork(1,prm%ws), mat_size)

       end select
    end do

    deallocate(rwork, cwork, C, ipiv)
  end subroutine


#else


!> @brief Dense driver routine for complex Hermitian matrix
!! @param[in,out] prm optional parameters
!! @param[in] UPLO = 'U': Upper triangle part of A and B are stored; = 'L': Lower triangle part of A and B are stored.
!! @param[in] mat_size matrix size
!! @param[in] A matrix A
!! @param[in] LDA Leading dimension of A
!! @param[in] emin Lower bound of the interval
!! @param[in] emax Upper bound of the interval
!! @param[out] num_ev number of eigenvalues
!! @param[out] eigval eigenvalues
!! @param[out] X eigenvectors
!! @param[out] res residuals
!! @param[out] info infomations
  subroutine MACRO_INSERT_PRFX(zpares_,dnsheev) &
  (prm, UPLO, mat_size, A, LDA, emin, emax, num_ev, eigval, X, res, info, set_rule)
    implicit none
#include "set_rule.f90"
    optional set_rule
    type(zpares_prm), intent(inout) :: prm
    integer, intent(in) :: mat_size, LDA
    integer, intent(inout) :: num_ev
    integer, intent(out) :: info
    character, intent(in) :: UPLO
    REAL_TYPE, intent(in) :: emin, emax
    REAL_TYPE, intent(out) :: res(*), eigval(*)
    COMPLEX_TYPE, intent(in) :: A(LDA,*)
    COMPLEX_TYPE, intent(out) :: X(mat_size,*)
    
    integer :: L, Lmax, M, i, j, infola
    integer, allocatable :: ipiv(:)
    COMPLEX_TYPE :: z
    COMPLEX_TYPE, allocatable :: rwork(:,:), cwork(:,:), C(:,:), tmpvec(:)

    L = prm%L
    Lmax = prm%Lmax
    M = prm%M

    allocate(rwork(mat_size, Lmax), cwork(mat_size, Lmax) &
         , C(mat_size, mat_size), ipiv(mat_size), tmpvec(mat_size))
    
    do while ( prm%itask /= ZPARES_TASK_FINISH )
       call MACRO_INSERT_PRFX(zpares_,rciheev) &
       (prm, mat_size, z, rwork, cwork, emin, emax, num_ev, eigval, X, res, info, set_rule)
       select case(prm%itask)
       case(ZPARES_TASK_FACTO)
          if ( UPLO == 'L' .or. UPLO == 'l' ) then

             do i = 1, mat_size
                do j = 1, i-1
                   C(i,j) = -A(i,j)
                   C(j,i) = -conjg(A(i,j))
                end do
                C(i,i) = z - A(i,i)
             end do
             
          elseif ( UPLO == 'U' .or. UPLO == 'u' ) then

             do i = 1, mat_size
                do j = 1, i-1
                   C(j,i) = -A(j,i)
                   C(i,j) = -conjg(A(j,i))
                end do
                C(i,i) = z - A(i,i)
             end do
             
          end if
          call MACRO_ADD_C_PRFX(GETRF) &
          (mat_size, mat_size, C, mat_size, ipiv, infola)

       case(ZPARES_TASK_SOLVE)

          call MACRO_ADD_C_PRFX(GETRS) &
          ('N', mat_size, prm%nc, C, mat_size, ipiv &
          , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_SOLVE_H)
          
          call MACRO_ADD_C_PRFX(GETRS) &
          ('C', mat_size, prm%nc, C, mat_size, ipiv &
          , cwork(1,prm%ws), mat_size, infola)

       case(ZPARES_TASK_MULT_A)
          
          call MACRO_ADD_PRFX(HEMM) &
          ('L', UPLO, mat_size, prm%nc &
          , ONE_C, A, LDA, X(1,prm%xs), mat_size &
          , ZERO_C, rwork(1,prm%ws), mat_size)

       end select

    end do
    
    deallocate(rwork, cwork, C, ipiv, tmpvec)
  end subroutine
#endif

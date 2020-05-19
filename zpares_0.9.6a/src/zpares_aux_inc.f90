#include "def_macros.h"

  subroutine MACRO_ADD_PRFX(_ALLREDUCE_SUM_1D) &
  (array, n, ierr, mpi_comm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: n
    integer, intent(out) :: ierr
    integer, intent(in), optional :: mpi_comm
    integer :: comm
    MATRIX_TYPE, intent(inout) :: array(*)
#ifdef MPI
    if ( present(mpi_comm) ) then
       comm = mpi_comm       
    else
       comm = MPI_COMM_SELF
    end if
    call MPI_ALLREDUCE(MPI_IN_PLACE, array, n, MPI_TYPE, MPI_SUM, comm, ierr)
#else
    ierr = 0
#endif
  end subroutine

  subroutine MACRO_ADD_PRFX(_ALLREDUCE_SUM_2D) &
  (array, ptr, m, n, ierr, mpi_comm)
    implicit none    
#ifdef MPI
    include 'mpif.h'
#endif    
    integer, intent(in) :: ptr, m, n
    integer, intent(out) :: ierr
    integer, intent(in), optional :: mpi_comm
    MATRIX_TYPE, intent(inout) :: array(m, *)
    integer :: comm
#ifdef MPI
    if ( present(mpi_comm ) ) then
       comm = mpi_comm
    else
       comm = MPI_COMM_SELF
    end if
    call MPI_ALLREDUCE(MPI_IN_PLACE, array(1,ptr), m*n, MPI_TYPE, MPI_SUM, comm, ierr)
#else
    ierr = 0
#endif
  end subroutine

  subroutine MACRO_ADD_PRFX(_norm2_blk) &
  (V, norm, nrow, ncol, ierr, mpi_comm)
    implicit none
    integer, intent(in) :: nrow, ncol, mpi_comm
    integer, intent(out) :: ierr
    MATRIX_TYPE, intent(in) :: V(nrow,ncol)
    REAL_TYPE, intent(out) :: norm(ncol)  
    integer :: i
    REAL_TYPE :: temp(ncol)

    norm = ZERO_R
    do i=1,nrow
#ifdef REALMAT
       temp = V(i,:)
#else
       temp = abs(V(i,:))
#endif
       norm = norm + temp*temp
    end do
    call ALLREDUCE_SUM_1D(norm, ncol, ierr, mpi_comm)
    norm = sqrt(norm)
  end subroutine

  subroutine MACRO_ADD_PRFX(GEMM_ALLREDUCE) &
  (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, ierr, mpi_comm)
    implicit none
    character, intent(in) :: transa, transb
    integer, intent(in) ::  m, n, k, lda, ldb, ldc
    integer, intent(out) :: ierr
    integer, intent(in), optional :: mpi_comm
    MATRIX_TYPE, intent(in) :: alpha, A(lda,*), B(ldb,*), beta
    MATRIX_TYPE, intent(inout) :: C(ldc,*)
    character :: transa2, transb2
    
    transa2 = transa
    transb2 = transb
#ifdef REALMAT
    if ( transa == 'C' .or. transa == 'c' ) transa2 = 'T'
    if ( transb == 'C' .or. transb == 'c' ) transb2 = 'T'
#endif

    call MACRO_ADD_PRFX(GEMM) &
    (transa2, transb2, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    
    call ALLREDUCE_SUM_2D(C, 1, ldc, n, ierr, mpi_comm)
  end subroutine


#ifdef REALMAT
  subroutine MACRO_ADD_PRFX(_quad_ell_trap) &
    (left, right, N, asp_ratio, quad_idx, zeta, weight, z)
    implicit none
    COMPLEX_TYPE, intent(in) :: left
    REAL_TYPE, intent(in) :: right, asp_ratio
    integer, intent(in) :: N, quad_idx
    COMPLEX_TYPE,intent(out) :: zeta, weight, z

    COMPLEX_TYPE :: i, t, center
    REAL_TYPE :: pi, radius

    call calc_center_radius(left, right, center, radius)
    pi = (ONE_R+ONE_R+ONE_R+ONE_R)*atan(ONE_R)
    i = (ZERO_R,ONE_R)
    t = (ONE_R+ONE_R) * pi / N * ((quad_idx - 1) + 1/(ONE_R+ONE_R))
    zeta = cos(t) + i * asp_ratio * sin(t)
    weight = radius * (i * sin(t) + asp_ratio * cos(t)) / N
    z = center + radius * zeta
  end subroutine
#endif

  subroutine MACRO_ADD_PRFX(_create_hutch_samples) &
  (V, nrow, ncol, rank)
    implicit none
    integer, intent(in) :: nrow, ncol, rank
    MATRIX_TYPE, intent(out) :: V(nrow,*)
    integer :: i, j
    call create_rand_matrix(V, nrow, ncol, rank)
    do j = 1, ncol
       do i = 1, nrow
          V(i,j) = cmplx(sign(ONE_R, real(V(i,j), kind(ZERO_R))), ZERO_R, kind(ZERO_R))
       end do
    end do    
  end subroutine

#ifdef REALMAT
  
  subroutine MACRO_ADD_PRFX(_create_rand_matrix) &
  (V, nrow, ncol, rank)
    implicit none
    integer, intent(in) :: nrow, ncol, rank
    REAL_TYPE, intent(out) :: V(nrow,*)
    integer :: iseed(4)
    iseed(1) = modulo(rank-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(rank-4096, 4096) ! must be between 0 and 4095
    iseed(3) = modulo(rank, 4096) ! must be between 0 and 4095
    iseed(4) = 1 ! must be between 0 and 4095 and odd
    call MACRO_ADD_PRFX(LARNV) &
         (2, iseed, nrow*ncol, V)    
  end subroutine

#else

  subroutine MACRO_ADD_PRFX(_create_rand_matrix) &
  (V, nrow, ncol, rank)
    implicit none
    integer, intent(in) :: nrow, ncol, rank
    COMPLEX_TYPE, intent(out) :: V(nrow,*)
    integer :: iseed(4)

    REAL_TYPE, allocatable :: tmpV(:,:)
    allocate(tmpV(nrow,ncol))

    iseed(1) = modulo(rank-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(rank-4096, 4096) ! must be between 0 and 4095
    iseed(3) = modulo(rank, 4096) ! must be between 0 and 4095
    iseed(4) = 1 ! must be between 0 and 4095 and odd
! #ifdef SINGLE
!     call CLARNV &
! #else
!     call ZLARNV &
! #endif
!     (2, iseed, nrow*ncol, V)    
#ifdef SINGLE
    call SLARNV &
#else
    call DLARNV &
#endif
    (2, iseed, nrow*ncol, tmpV)
    V(:,1:ncol) = tmpV(:,1:ncol)
  end subroutine
#endif

  subroutine MACRO_ADD_PRFX(_orth_SVD) &
  (nrow_local, ncol, Q, ts_work, wrk_ncol, sqmat, LDSQ, first_time, mpi_comm, infola, mpierr)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    integer, intent(in) :: nrow_local, ncol, wrk_ncol, mpi_comm, LDSQ
    integer, intent(out) :: infola, mpierr
    logical, intent(in) :: first_time
    MATRIX_TYPE, intent(inout) :: Q(nrow_local,*), sqmat(LDSQ,*)    
    MATRIX_TYPE, intent(out) :: ts_work(nrow_local,*) ! just work arrays
    REAL_TYPE, allocatable :: sigma(:)
    COMPLEX_TYPE, allocatable :: csigma(:)
    MATRIX_TYPE, allocatable :: work(:), tmpmat(:,:), tmpmat2(:,:)
    integer :: i, j, width
    REAL_TYPE :: tmp_r, epsilon
    MATRIX_TYPE :: tmp
    REAL_TYPE, external :: MACRO_ADD_R_PRFX(LAMCH)
    
    infola = 0
    mpierr = 0
    allocate(tmpmat(ncol, ncol), tmpmat2(ncol, ncol), sigma(ncol), csigma(ncol))

    do i = 1, ncol, wrk_ncol
       width = min(ncol - i + 1, wrk_ncol)
       j = i + width - 1
       ts_work(:,1:width) = Q(:,i:j)
       call GEMM_ALLREDUCE('C', 'N', ncol, width, nrow_local, ONE_M, Q, nrow_local &
            , ts_work, nrow_local, ZERO_M, tmpmat(:,i:j), ncol, mpierr, mpi_comm)
       if ( mpierr /= 0 ) then
          return
       end if
    end do

    call HEEV_reduced_eig(nrow_local, ncol, tmpmat, ncol, csigma, infola)
    ! if ( infola /= 0 ) then
    !    return
    ! end if
    sigma(:) = real(csigma(:), kind(ZERO_R))

    ! reverse order
    do i = 1, ncol/2
       tmp_r = sigma(i)
       sigma(i) = sigma(ncol-i+1)
       sigma(ncol-i+1) = tmp_r
       do j = 1, ncol
          tmp = tmpmat(j,i)
          tmpmat(j,i) = tmpmat(j,ncol-i+1)
          tmpmat(j,ncol-i+1) = tmp
       end do       
    end do

    epsilon = sigma(1) * MACRO_ADD_R_PRFX(LAMCH)('E')

    do i = 1, ncol
       if ( sigma(i) < ZERO_R ) then
          sigma(i) = epsilon
       end if
    end do
    
    sigma(1:ncol) = sqrt(abs(sigma(1:ncol)))

    ! Q := Q*tmpmat
    call basis_rotation('N', nrow_local, ncol, ncol, tmpmat, ncol &
         , ts_work, nrow_local*wrk_ncol/ncol, Q)

    do i = 1, ncol
       Q(:,i) = Q(:,i) / sigma(i)
    end do

    if ( first_time ) then
       sqmat(1:ncol,1:ncol) = tmpmat(1:ncol,1:ncol)
    else
       ! sqmat := sqmat*tmpmat
       tmpmat2(1:ncol,1:ncol) = sqmat(1:ncol,1:ncol)
       call GEMM_ALLREDUCE('N', 'N', ncol, ncol, ncol, ONE_M, tmpmat2, ncol, tmpmat, ncol &
            , ZERO_M, sqmat, LDSQ, mpierr)
    end if
    do i = 1, ncol
       sqmat(:,i) = sqmat(:,i) * sigma(i)
    end do

    deallocate(tmpmat, tmpmat2, sigma, csigma)
  end subroutine

  subroutine MACRO_ADD_PRFX(DOT_ALLREDUCE) &
  (X, Y, nrow_local, ncol, mpi_comm, dot, ierr)
    implicit none
    integer, intent(in) :: nrow_local, ncol, mpi_comm
    integer, intent(out) :: ierr
    MATRIX_TYPE, intent(in) :: X(nrow_local, *), Y(nrow_local, *)
    MATRIX_TYPE, intent(out) :: dot(*)
    integer :: i, j
    
    dot(1:ncol) = ZERO_M
    do j = 1, ncol
       do i = 1, nrow_local
#ifdef REALMAT
          dot(j) = dot(j) + x(i,j)*y(i,j)
#else
          dot(j) = dot(j) + conjg(x(i,j))*y(i,j)
#endif
       end do
    end do
    call ALLREDUCE_SUM_1D(dot, ncol, ierr, mpi_comm)
  end subroutine


  subroutine MACRO_ADD_PRFX(_block_Hankel) &
  (Lmax, L, M, shift, Mu, H)
    implicit none
    integer, intent(in) :: Lmax, L, M, shift
    MATRIX_TYPE, intent(in) :: Mu(Lmax,*)
    MATRIX_TYPE, intent(out) :: H(L*M,*)
    integer :: i
    do i = 1, M
       H((i-1)*L+1:i*L,1:L*M) = Mu(1:L, (i+shift-1)*L+1:(i+shift-1)*L+L*M)
    end do
 end subroutine


  subroutine MACRO_ADD_PRFX(_serial_SVD) &
  (NLRB, nrow, ncol, H, LDH, delta, sigma, U, LDU, VT, LDVT, num_rank, info)
    implicit none
    character, intent(in) :: NLRB
    integer, intent(in) :: nrow, ncol, LDH, LDU, LDVT
    REAL_TYPE, intent(in) :: delta
    MATRIX_TYPE, intent(in) :: H(LDH,*)
    REAL_TYPE, intent(out) :: sigma(*)
    MATRIX_TYPE, intent(out) :: U(LDU,*), VT(LDVT,*)
    integer, intent(out) :: num_rank, info
    integer :: sig_size, lwork, infola
    REAL_TYPE :: sigma_max, tmp_r
    MATRIX_TYPE, allocatable :: work(:)
    MATRIX_TYPE :: optlwork
    character :: jobU, jobVT
#ifndef REALMAT
    REAL_TYPE, allocatable :: rwork(:)
    allocate(rwork(5*max(nrow,ncol)))
#endif
    sig_size = min(nrow,ncol)

    if ( NLRB == 'N' ) then
       jobU = 'N'
       jobVT = 'N'
    else if ( NLRB == 'L' ) then
       jobU = 'O'
       jobVT = 'N'
    else if ( NLRB == 'R' ) then
       jobU = 'N'
       jobVT = 'O'
    else if ( NLRB == 'B' ) then
       jobU = 'S'
       jobVT = 'S'
    end if

    call MACRO_ADD_PRFX(GESVD) &
    (jobU, jobVT, nrow, ncol, H, LDH, sigma, U, LDU, VT, LDVT, optlwork, -1 &
#ifdef REALMAT
    , infola)
#else
    , rwork, infola)
#endif
    lwork = int(optlwork)
    allocate(work(lwork))

    call MACRO_ADD_PRFX(GESVD) &
    (jobU, jobVT, nrow, ncol, H, LDH, sigma, U, LDU, VT, LDVT,  work, lwork &
#ifdef REALMAT
    , infola)
#else
    , rwork, infola)
#endif    
    sigma_max = sigma(1)
    do num_rank = 1, sig_size
       if ( sigma(num_rank) < delta*sigma_max ) then
          exit
       end if
    end do
    num_rank = num_rank - 1
    deallocate(work)
#ifndef REALMAT
    deallocate(rwork)
#endif
    info = 0 !!!!! TODO: error handling
  end subroutine


  subroutine MACRO_ADD_PRFX(_LAPACK_QR) &
  (nrow, ncol, Q, R)
    implicit none
    integer, intent(in) :: nrow, ncol
    MATRIX_TYPE, intent(inout) :: Q(nrow,*)
    MATRIX_TYPE, intent(out) :: R(ncol,*)
    MATRIX_TYPE, allocatable :: TAU(:), work(:)
    MATRIX_TYPE :: optlwork
    integer :: info, i, j, lwork

    allocate(TAU(min(nrow, ncol)))
    
    call MACRO_ADD_PRFX(GEQRF) &
    (nrow, ncol, Q, nrow, TAU, optlwork, -1, info)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(GEQRF) &
    (nrow, ncol, Q, nrow, TAU, work, lwork, info)
    deallocate(work)
    do i = 1, ncol
       do j = 1, ncol
          if ( i <= j ) then
             R(i,j) = Q(i,j)
          else
             R(i,j) = ZERO_R
          end if
       end do
    end do
#ifdef REALMAT
    call MACRO_ADD_PRFX(ORGQR) &
#else
    call MACRO_ADD_PRFX(UNGQR) &
#endif
    (nrow, ncol, ncol, Q, nrow, TAU, optlwork, -1, info)
    lwork = int(optlwork)
    allocate(work(lwork))
#ifdef REALMAT
    call MACRO_ADD_PRFX(ORGQR) &
#else
    call MACRO_ADD_PRFX(UNGQR) &
#endif
    (nrow, ncol, ncol, Q, nrow, TAU, work, lwork, info)
    deallocate(TAU , work)
  end subroutine


#ifdef REALMAT
  subroutine MACRO_ADD_PRFX(GEGV_reduced_eig) &
  (nrow_local, nb, A, LDA, B, LDB, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, LDB, nb
    integer, intent(out) :: info
    REAL_TYPE, intent(inout) :: A(LDA,*), B(LDB,*)
    COMPLEX_TYPE, intent(out) ::  eigval(nb)
    integer :: i, lwork, infola
    REAL_TYPE, allocatable :: VR(:,:), alpha_r(:), alpha_i(:), beta(:), work(:)
    REAL_TYPE :: dummy_VL(1,1), optlwork

    allocate(VR(nb, nb), alpha_r(nb), alpha_i(nb), beta(nb))
    call MACRO_ADD_PRFX(GEGV) &
    ('N', 'V', nb, A, LDA, B, LDB, alpha_r, alpha_i, beta, dummy_VL, 1, VR, nb, optlwork, -1, infola)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(GEGV) &
    ('N', 'V', nb, A, LDA, B, LDB, alpha_r, alpha_i, beta, dummy_VL, 1, VR, nb, work, lwork, infola)
    eigval(1:nb) = cmplx(alpha_r(1:nb)/beta(1:nb), alpha_i(1:nb)/beta(1:nb), kind(ZERO_R))
    A(1:nb,1:nb) = VR(1:nb,1:nb)
    deallocate(VR, alpha_r, alpha_i, beta, work)
    info = 0 !!!!! TODO: error handling
  end subroutine

#else

  subroutine MACRO_ADD_PRFX(GEGV_reduced_eig) &
  (nrow_local, nb, A, LDA, B, LDB, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, LDB, nb
    integer, intent(out) :: info
    COMPLEX_TYPE, intent(inout) :: A(LDA,*), B(LDB,*)
    COMPLEX_TYPE, intent(out) ::  eigval(nb)
    integer :: lwork, infola
    REAL_TYPE, allocatable :: rwork(:)
    COMPLEX_TYPE, allocatable :: VR(:,:), beta(:), work(:)
    COMPLEX_TYPE :: dummy_VL(1,1), optlwork

    allocate(VR(nb, nb), beta(nb), rwork(8*nb))
    call MACRO_ADD_PRFX(GEGV) &
    ('N', 'V', nb, A, LDA, B, LDB, eigval, beta, dummy_VL, 1, VR, nb, optlwork, -1, rwork, infola)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(GEGV) &
    ('N', 'V', nb, A, LDA, B, LDB, eigval, beta, dummy_VL, 1, VR, nb, work, lwork, rwork, infola)
    eigval(1:nb) = eigval(1:nb) / beta(1:nb)
    A(1:nb,1:nb) = VR(1:nb,1:nb)
    deallocate(VR, beta, work, rwork)
    info = 0 !!!!! TODO: error handling
  end subroutine

#endif


#ifdef REALMAT

  subroutine MACRO_ADD_PRFX(GEEV_reduced_eig) &
  (nrow_local, nb, A, LDA, eigval, info)
    implicit none
    integer, intent(in) :: nrow_local, LDA, nb
    integer, intent(out) :: info
    REAL_TYPE, intent(inout) :: A(LDA,*)
    COMPLEX_TYPE, intent(out) ::  eigval(*)
    integer :: i, lwork, infola
    REAL_TYPE, allocatable :: VR(:,:), wr(:), wi(:), work(:)
    REAL_TYPE :: dummy_VL(1,1), optlwork

    allocate(VR(nb, nb), wr(nb), wi(nb))    
    call MACRO_ADD_PRFX(GEEV) &
    ('N', 'V', nb, A, LDA, wr, wi, dummy_VL, 1, VR, nb, optlwork, -1, infola)
    lwork  = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(GEEV) &
    ('N', 'V', nb, A, LDA, wr, wi, dummy_VL, 1, VR, nb, work, lwork, infola)
    eigval(1:nb) = cmplx(wr(1:nb), wi(1:nb), kind(ZERO_R))
    A(1:nb,1:nb) = VR(1:nb,1:nb)
    deallocate(VR, wr, wi, work)
    info = 0 !!!!! TODO: error handling
  end subroutine

#else

  subroutine MACRO_ADD_PRFX(GEEV_reduced_eig) &
  (nrow_local, nb, A, LDA, eigval, info)
    implicit none
    integer, intent(in) :: nrow_local, LDA, nb
    integer, intent(out) :: info
    COMPLEX_TYPE, intent(inout) :: A(LDA,*)
    COMPLEX_TYPE, intent(out) ::  eigval(*)
    integer :: lwork, infola
    REAL_TYPE,allocatable :: rwork(:)
    COMPLEX_TYPE,allocatable :: VR(:,:), work(:)
    COMPLEX_TYPE :: dummy_VL(1,1), optlwork

    allocate(VR(nb, nb), rwork(2*nb))
    call MACRO_ADD_PRFX(GEEV) &
    ('N', 'V', nb, A, LDA, eigval, dummy_VL, 1, VR, nb, optlwork, -1, rwork, infola)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(GEEV) &
    ('N', 'V', nb, A, LDA, eigval, dummy_VL, 1, VR, nb, work, lwork, rwork, infola)
    A(1:nb,1:nb) = VR(1:nb,1:nb)
    deallocate(VR, work, rwork)
    info = 0 !!!!! TODO: error handling
  end subroutine

#endif


#ifdef REALMAT

  subroutine MACRO_ADD_PRFX(SYGV_reduced_eig) &
    (nrow_local, nb, A, LDA, B, LDB, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, LDB
    integer, intent(inout) :: nb
    integer, intent(out) :: info
    REAL_TYPE, intent(inout) :: A(LDA,*), B(LDB,*)
    COMPLEX_TYPE, intent(out) ::  eigval(*)
    integer :: lwork, counter, infola
    REAL_TYPE, allocatable :: r_eigval(:), work(:), tmpB(:,:)
    REAL_TYPE :: optlwork

    allocate(r_eigval(nb))    
    infola = -1
    counter = 0
    do while ( infola /= 0 .and. counter <= 5 ) 
       allocate(tmpB(nb,nb))
       tmpB(1:nb,1:nb) = B(1:nb,1:nb)
       call MACRO_ADD_PRFX(SYGV) &
            (1, 'V', 'U', nb, A, LDA, tmpB, nb, r_eigval, optlwork, -1, infola)
       lwork = int(optlwork)
       allocate(work(lwork))
       call MACRO_ADD_PRFX(SYGV) &
            (1, 'V', 'U', nb, A, LDA, tmpB, nb, r_eigval, work, lwork, infola)
       if ( infola > nb ) then
          nb = infola - nb - 1
       end if       
       deallocate(tmpB, work)
       counter = counter + 1
    end do
    eigval(1:nb) = cmplx(r_eigval(1:nb), ZERO_R, kind(ZERO_R))
    deallocate(r_eigval)
    info = 0 !!!!! TODO: error handling
  end subroutine

#else

  subroutine MACRO_ADD_PRFX(HEGV_reduced_eig) &
  (nrow_local, nb, A, LDA, B, LDB, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, LDB
    integer, intent(inout) :: nb
    integer, intent(out) :: info
    COMPLEX_TYPE, intent(inout) :: A(LDA,*), B(LDB,*)
    COMPLEX_TYPE, intent(out) ::  eigval(*)
    integer :: lwork, infola, counter
    REAL_TYPE, allocatable :: rwork(:), r_eigval(:)
    COMPLEX_TYPE,allocatable :: work(:), tmpB(:,:)
    COMPLEX_TYPE :: optlwork

    allocate(r_eigval(nb))
    infola = -1
    counter = 0
    do while ( infola /= 0 .and. counter <= 5 )
       allocate(tmpB(nb,nb), rwork(3*nb-2))
       tmpB(1:nb,1:nb) = B(1:nb,1:nb)
       call MACRO_ADD_PRFX(HEGV) &
            (1, 'V', 'U', nb, A, LDA, tmpB, nb, r_eigval, optlwork, -1, rwork, infola)
       lwork = int(optlwork)
       allocate(work(lwork))
       call MACRO_ADD_PRFX(HEGV) &
            (1, 'V', 'U', nb, A, LDA, tmpB, nb, r_eigval, work, lwork, rwork, infola)
       if ( infola > nb ) then
          nb = infola - nb - 1
       end if       
       deallocate(tmpB, work, rwork)
       counter = counter + 1
    end do
    eigval(1:nb) = cmplx(r_eigval(1:nb), ZERO_R, kind(ZERO_R))
    deallocate(r_eigval)
    info = 0 !!!!! TODO: error handling
  end subroutine
#endif


#ifdef REALMAT
  subroutine MACRO_ADD_PRFX(SYEV_reduced_eig) &
  (nrow_local, nb, A, LDA, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, nb
    integer, intent(out) :: info
    REAL_TYPE, intent(inout) :: A(LDA,*)
    COMPLEX_TYPE, intent(out) ::  eigval(nb)
    integer :: lwork, infola
    REAL_TYPE, allocatable :: r_eigval(:)
    REAL_TYPE, allocatable :: work(:)
    REAL_TYPE :: optlwork

    allocate(r_eigval(nb))
    call MACRO_ADD_PRFX(SYEV) &
    ('V', 'U', nb, A, LDA, r_eigval, optlwork, -1, infola)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(SYEV) &
    ('V', 'U', nb, A, LDA, r_eigval, work, lwork, infola)
    eigval = cmplx(r_eigval, ZERO_R, kind(ZERO_R))
    deallocate(work, r_eigval)    
    info = 0 !!!!! TODO: error handling
  end subroutine

#else

  subroutine MACRO_ADD_PRFX(HEEV_reduced_eig) &
  (nrow_local, nb, A, LDA, eigval, info)
    implicit none
    integer, intent(in) ::  nrow_local, LDA, nb
    integer, intent(out) :: info
    COMPLEX_TYPE, intent(in) :: A(LDA,*)
    COMPLEX_TYPE, intent(out) ::  eigval(nb)
    integer :: lwork, infola
    REAL_TYPE, allocatable :: rwork(:), r_eigval(:)
    COMPLEX_TYPE,allocatable :: work(:)
    COMPLEX_TYPE :: optlwork

    allocate(rwork(3*nb-2), r_eigval(nb))
    call MACRO_ADD_PRFX(HEEV) &
    ('V', 'U', nb, A, LDA, r_eigval, optlwork, -1, rwork, infola)
    lwork = int(optlwork)
    allocate(work(lwork))
    call MACRO_ADD_PRFX(HEEV) &
    ('V', 'U', nb, A, LDA, r_eigval, work, lwork, rwork, infola)
    eigval = cmplx(r_eigval, ZERO_R, kind(ZERO_R))
    deallocate(work, rwork, r_eigval)
    info = 0 !!!!! TODO: error handling
  end subroutine
#endif


  subroutine MACRO_ADD_PRFX(_basis_rotation) &
  (transRV, m, n, k, reduced_vec, LDRV, work, LDWK, orgvec)
    implicit none
    character, intent(in) :: transRV
    integer, intent(in) ::  m, n, k, LDWK, LDRV
    MATRIX_TYPE, intent(in) :: reduced_vec(LDRV,*)
    MATRIX_TYPE, intent(inout) :: orgvec(m,*)
    MATRIX_TYPE, intent(out) :: work(LDWK,*) ! just work arrays
    character :: transRV2

    integer :: i, gemm_m
    
    transRV2 = transRV
#ifdef REALMAT
    if ( transRV == 'C' ) then
       transRV2 = 'T'
    end if
#endif

    do i = 1, m, LDWK
       gemm_m = min(m - i + 1, LDWK)
       work(1:gemm_m,1:k) = orgvec(i:i+gemm_m-1,1:k)
       call MACRO_ADD_PRFX(GEMM) &
       ('N', transRV, gemm_m, n, k, ONE_M, work, LDWK, reduced_vec, LDRV, ZERO_M, orgvec(i,1), m)
       ! orgvec(i:i+gemm_m-1,1:n) = work2(1:gemm_m,1:n)
    end do
  end subroutine

#ifdef REALMAT
  subroutine MACRO_ADD_PRFX(_inside_ellipse) &
  (left, right, asp_ratio, m, eigval, flags, num_true)
    implicit none
    COMPLEX_TYPE, intent(in) :: left
    REAL_TYPE, intent(in) :: right, asp_ratio
    COMPLEX_TYPE, intent(in) :: eigval(*)
    integer, intent(in) :: m
    logical, intent(out) :: flags(*)
    integer, intent(out) :: num_true    
    integer :: i
    COMPLEX_TYPE :: z, center
    REAL_TYPE :: radius
    
    call calc_center_radius(left, right, center, radius)
    num_true = 0
    do i = 1, m
       z = (eigval(i) - center) / radius
       if ( real(z,kind(ZERO_R))**2 + aimag(z)**2/(asp_ratio**2) <= ONE_R ) then
          flags(i) = .true.
          num_true = num_true + 1
       else
          flags(i) = .false.
       end if
    end do
  end subroutine
#endif
  
  subroutine MACRO_ADD_PRFX(_packing) &
  (m, flags, eigval, eigvec, LDEV, res, indi_spu)
    implicit none
    integer, intent(in) :: m, LDEV
    logical, intent(in) :: flags(:)
    COMPLEX_TYPE, intent(inout) :: eigval(*)
    MATRIX_TYPE, intent(inout) :: eigvec(LDEV,*)
    REAL_TYPE, intent(inout), optional :: res(*)
    double precision, intent(inout), optional :: indi_spu(*)    
    integer :: i, ptr
    
    ptr = 1
    do i = 1, m
       if ( flags(i) ) then
          if ( ptr /= i ) then
             eigval(ptr) = eigval(i)
             eigvec(:,ptr) = eigvec(:,i)
             if ( present(res) ) then
                res(ptr) = res(i)
             end if
             if ( present(indi_spu) ) then
                indi_spu(ptr) = indi_spu(i)
             end if
          end if
          ptr = ptr + 1
       end if
    end do
  end subroutine

#ifdef REALMAT
  subroutine MACRO_ADD_PRFX(_calc_center_radius)&
       (left, right, center, radius)
    implicit none
    COMPLEX_TYPE, intent(in) :: left
    REAL_TYPE, intent(in) :: right
    COMPLEX_TYPE, intent(out) :: center
    REAL_TYPE, intent(out) :: radius

    center = cmplx((real(left,kind(ZERO_R)) + right)/(ONE_R+ONE_R), aimag(cmplx(left)),kind(ZERO_R))
    radius = (right - real(left,kind(ZERO_R))) / (ONE_R+ONE_R) 
  end subroutine
#endif

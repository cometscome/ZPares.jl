! A and B are real symmetric and B is positive definite
program main
  use zpares
  use zpares_mumps
  implicit none  
#ifdef MPI
  include 'mpif.h'
#endif
  integer, parameter :: mat_size = 500
  integer :: num_ev, info, i, j, L, M, N, Lmax, ncv, ierr, nnzA, nnzB
  double precision :: emin, emax
  integer, allocatable :: rowptrA(:), rowptrB(:), colindA(:), colindB(:)
  double precision, allocatable :: res(:), eigval(:)
  double precision, allocatable :: valA(:), valB(:), X(:,:)
  type(zpares_prm) :: prm
  
#ifdef MPI
  call MPI_INIT(ierr)
#endif

  nnzA = mat_size + (mat_size-1) ! nnz of upper triangular part of A
  nnzB = mat_size ! B is the identity matrix
  allocate(rowptrA(mat_size+1), rowptrB(mat_size+1), colindA(nnzA), colindB(nnzB))
  allocate(valA(nnzA), valB(nnzB))
  
  ! make upper triangular part of A
  j = 1
  rowptrA(1) = 1
  do i = 1, mat_size     
     rowptrA(i+1) = rowptrA(i)
     colindA(j) = i
     valA(j) = 2d0 
     rowptrA(i+1) = rowptrA(i+1) + 1
     j = j + 1
     if ( i /= mat_size ) then
        colindA(j) = i+1
        valA(j) = 1d0
        rowptrA(i+1) = rowptrA(i+1) + 1
        j = j + 1     
     end if
  end do
  
  ! make B
  do i = 1, mat_size
     rowptrB(i) = i
     colindB(i) = i
     valB(i) = 1d0
  end do
  rowptrB(mat_size+1) = mat_size+1
  
  L = 8
  N = 32
  M = 16
  Lmax = 32
  
  call zpares_init(prm)
  prm%L = L
  prm%N = N
  prm%M = M
  prm%Lmax = Lmax

  ncv = zpares_get_ncv(prm)
  allocate(eigval(ncv), X(mat_size, ncv), res(ncv))
  emin = 3.99d0
  emax = 4.01d0
    
  call zpares_dmpssygv(prm, mat_size, rowptrA, colindA, valA, rowptrB, colindB, valB &
       , emin, emax, num_ev, eigval, X, res, info)  
  call zpares_finalize(prm)
  do i = 1, num_ev
     write(*,*) i, eigval(i), res(i)
  end do
  
#ifdef MPI  
  call MPI_FINALIZE(ierr)
#endif
  
  deallocate(rowptrA, colindA, valA, rowptrB, colindB, valB, eigval, X, res)
end program main

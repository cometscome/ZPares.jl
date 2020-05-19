program main
  use zpares
  implicit none
  include 'mpif.h'
  integer, parameter :: nrow_local = 500
  integer :: num_ev, info, i, j, L, M, N, Lmax, ncv, ierr, rank
  double precision :: right
  double precision, allocatable :: res(:)
  complex(kind(0d0)) :: left
  complex(kind(0d0)), allocatable :: A(:,:), B(:,:), eigval(:), X(:,:)
  type(zpares_prm) :: prm
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  allocate(A(nrow_local,nrow_local), B(nrow_local,nrow_local)) 
    
  A = (0d0,0d0)
  do i = 1, nrow_local
     do j = 1, nrow_local
        if ( i == j ) then
           A(i,j) = (2d0,0d0)
        else if ( abs(i-j) == 1 ) then
           A(i,j) = (1d0,0d0)
        end if
     end do
  end do
  
  B = (0d0,0d0)
  do i = 1, nrow_local
     B(i,i) = (1d0,0d0)
  end do
  
  L = 8
  N = 32
  M = 16
  Lmax = 32

  call zpares_init(prm)
  prm%L = L
  prm%N = N
  prm%M = M
  prm%Lmax = Lmax

  prm%high_comm = MPI_COMM_SELF
  prm%low_comm = MPI_COMM_WORLD

  ncv = zpares_get_ncv(prm)
  allocate(eigval(ncv), X(nrow_local, ncv), res(ncv))
  left = (3.99d0,0d0)
  right = 4.01d0

  !! When nproc = 2 
  !! solve eigenvalue problem
  !! [A 0] x = lambda [B 0] x
  !! [0 A]            [0 B]

  !! When nproc = 3
  !! solve eigenvalue problem
  !! [A 0 0] x = lambda [B 0 0] x
  !! [0 A 0]   = lambda [0 B 0] 
  !! [0 0 A]   = lambda [0 0 B] 

  call zpares_zdnsgegv(prm, nrow_local, A, nrow_local, B, nrow_local, left, right, num_ev, eigval, X, res, info)  
  call zpares_finalize(prm)
  if ( rank == 0 ) then
     do i = 1, num_ev
        write(*,*) i, eigval(i), res(i)
     end do
  end if
  
  call MPI_FINALIZE(ierr)
  
  deallocate(A, B, eigval, X, res)
end program main

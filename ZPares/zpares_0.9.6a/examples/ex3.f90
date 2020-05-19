program main
  use zpares
  implicit none
  include 'mpif.h'
  integer, parameter :: mat_size = 500
  integer :: num_ev, info, i, j, L, M, N, Lmax, ncv, ierr, rank
  double precision :: right, time_s, time_e
  double precision, allocatable :: res(:)
  complex(kind(0d0)) :: left
  complex(kind(0d0)), allocatable :: A(:,:), B(:,:), eigval(:), X(:,:)
  type(zpares_prm) :: prm

  allocate(A(mat_size,mat_size), B(mat_size,mat_size)) 
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  
  A = (0d0,0d0)
  do i = 1, mat_size
     do j = 1, mat_size
        if ( i == j ) then
           A(i,j) = (2d0,0d0)
        else if ( abs(i-j) == 1 ) then
           A(i,j) = (1d0,0d0)
        end if
     end do
  end do
  
  B = (0d0,0d0)
  do i = 1, mat_size
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

  prm%high_comm = MPI_COMM_WORLD
  prm%low_comm = MPI_COMM_SELF

  ncv = zpares_get_ncv(prm)
  allocate(eigval(ncv), X(mat_size, ncv), res(ncv))
  left = (3.99d0,0d0)
  right = 4.01d0

  time_s = MPI_WTIME()
  call zpares_zdnsgegv(prm, mat_size,  A, mat_size, B, mat_size, left, right, num_ev, eigval, X, res, info)  
  call zpares_finalize(prm)
  time_e = MPI_WTIME()
  if ( rank == 0 ) then
     do i = 1, num_ev
        write(*,*) i, eigval(i), res(i)
     end do
     write(*,*) 'time = ', time_e - time_s
  end if

  call MPI_FINALIZE(ierr)
  
  deallocate(A, B, eigval, X, res)
end program main

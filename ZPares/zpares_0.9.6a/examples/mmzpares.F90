program main
  use zpares
#ifdef MUMPS
  use zpares_mumps
#endif
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  type option
     character :: id*20
     character :: usage*20
     character :: description*512
     type(option), pointer :: next
  end type option
  
  type(zpares_prm) :: prm
  character :: argv*512, str*512, filepathA*512, filepathB*512, sign*1
  character :: rep*10, field*7, symm*19
  integer :: num_ev, info, i, j, L, M, Lmax, ncv, mat_size, idx, fileunit, iarg, rank, nprocs, ierr
  integer :: rowsA, rowsB, cols, entries, nnzA, nnzB, nnzmax
  integer, allocatable :: indx(:), jndx(:), ival(:)
  integer, pointer :: rowptrA(:), colindA(:), rowptrB(:), colindB(:)
  double precision :: right
  double precision :: darg, rleft, ileft
  real, allocatable :: ress(:)
  real, allocatable :: eigvals(:), XS(:,:)
  real, pointer :: AS(:,:), BS(:,:), svalA(:), svalB(:)
  double precision, allocatable :: resd(:), rval(:)
  double precision, allocatable :: eigvald(:), XD(:,:)
  double precision, pointer :: AD(:,:), BD(:,:), dvalA(:), dvalB(:)
  complex(kind(0d0)), allocatable :: cval(:)
  complex(kind(0d0)) :: left
  complex, allocatable :: eigvalc(:), XC(:,:)
  complex, pointer :: AC(:,:), BC(:,:), cvalA(:), cvalB(:)
  complex(kind(0d0)), allocatable :: eigvalz(:), XZ(:,:)
  complex(kind(0d0)), pointer :: AZ(:,:), BZ(:,:), zvalA(:), zvalB(:)
  logical :: generalized, A_given, B_given, realmat, hhpd, larg, single, dense, standard, sigval
  logical :: option_found, dummy, error = .false.
  integer :: iargc, itc, t_rate, t_max
  type(option), pointer :: optlist
  integer, parameter :: root = 0
  character, parameter :: UPLO = 'L' ! U or L
  logical, parameter :: UPLOtest = .false.
  double precision :: total_time = 0d0, dtc

  rank = root
  
#ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
#endif

  allocate(optlist)
  nullify(optlist%next)
  
  generalized = .false.
  A_given = .false.
  B_given = .false.
  realmat = .false.
  hhpd = .false.
  single = .false.
  dense = .false.
  standard = .false.
  fileunit = 100
  rleft = -1d0
  ileft = 0d0
  right = 1d0
  sigval = .false.

  call zpares_init(prm)
  
  do i = 1, iargc() ! this function is a compilar extension
     option_found = .false.
     call getarg(i, argv) ! this function is a compilar extension
     
     if ( get_str(argv, '-fileA', '-fileA=<file>', 'Filepath for Matrix Market file of matrix A', filepathA) ) then
        
        A_given = .true.
        
     end if
     if ( get_str(argv, '-fileB', '-fileB=<file>', 'Filepath for Matrix Market file of matrix B', filepathB) ) then

        B_given = .true.
        
     end if
     if ( get_int(argv, '-N', '-N=<integer>', 'Number of quadrature points | N > 0', iarg) ) then
        
        prm%N = iarg
        
     end if
     if ( get_int(argv, '-M', '-M=<integer>', 'Number of complex moments | M > 0', iarg) ) then
        
        prm%M = iarg
        
     end if
     if ( get_int(argv, '-L', '-L=<integer>', 'Number of columns of the source matrix | 0 < L < Lmax', iarg) ) then
        
        prm%L = iarg
        
     end if
     if ( get_int(argv, '-Lmax', '-Lmax=<integer>', 'Maximum of L | 0 < L < Lmax',iarg) ) then
        
        prm%Lmax = iarg
        
     end if
     if ( get_int(argv, '-imax', '-imax=<integer>', 'Maximum number of iterative refinement | 0 < imax', iarg) ) then
        
        prm%imax = iarg
        
     end if
     if ( get_int(argv, '-n_orth', '-n_orth=<integer>' &
          , 'Maximum number of iteration for orthogonalization | 0 < n_orth', iarg) ) then
        
        prm%n_orth = iarg
        
     end if
     if ( get_double(argv, '-emin', '-emin=<Real*8>' &
          , 'Lower bound of the interval : use -rleft and -ileft for complex', darg) ) then
        
        rleft = darg
        
     end if
     if ( get_double(argv, '-emax', '-emax=<Real*8>' &
          , 'Upper bound of the interval', darg) ) then
        
        right = darg
        
     end if
     if ( get_double(argv, '-rleft', '-rleft=<Real*8>', 'Real part of left edge of circle', darg) ) then
        
        rleft = darg
        
     end if
     if ( get_double(argv, '-ileft', '-ileft=<Real*8>', 'Imaginary part of left edge of circle',darg) ) then
        
        ileft = darg
        
     end if
     if ( get_double(argv, '-right', '-right=<Real*8>', 'Real part of right edge of circle | right > ileft ', darg) ) then
        
        right = darg
        
     end if
     if ( get_double(argv, '-delta', '-delta=<Real*8>', 'Threshold for numerical rank | delta > 0', darg) ) then
        
        prm%delta = darg
        
     end if
     if ( get_double(argv, '-asp_ratio', '-asp_ratio=<Real*8>' &
          , 'Aspect ratio of the ellipse that defines the contour path | asp_ratio > 0', darg) ) then
        
        prm%asp_ratio = darg
        
     end if
     if ( get_double(argv, '-tol', '-tol=<Real*8>','Tolerance for residuals | tol >= 0', darg) ) then
        
        prm%tol = darg
        
     end if
     if ( get_logical(argv, '-EM', '-EM', 'Use this option if you want to use explicit moment mode') ) then
        
        prm%extract = ZPARES_EXTRACT_EM
        
     end if
     if ( get_logical(argv, '-realmat', '-realmat', 'Use this option if your matrices are real') ) then
        
        realmat = .true.
        
     end if
     if ( get_logical(argv, '-hhpd', '-hhpd' &
          , 'Use this option if A is Hermitian and B is Hermitian and positive definite' ) ) then
        
        hhpd = .true.
        
     end if
     if ( get_logical(argv, '-sspd', '-sspd' &
          , 'Use this option if A is real symmetric and B is real symmetric and positive definite' ) ) then
        
        hhpd = .true.
        realmat = .true.
        
     end if
     if ( get_logical(argv, '-single', '-single', 'Use this option if you would like to use single precision') ) then
        
        single = .true.
        
     end if
     if ( get_logical(argv, '-dense', '-dense', 'Use this option if you would like to use a dense routine') ) then
        
        dense = .true.
        
     end if
     if ( get_logical(argv, '-outside', '-outside', 'Show eigenvalues outside of the ellipse') ) then
        
        prm%trim_out = .false.
        
     end if
     if ( get_logical(argv, '-spurious', '-spurious', 'Show spurious eigenvalues') ) then
        
        prm%trim_spu = .false.
        
     end if
     if ( get_double(argv, '-spu_thres', '-sputhres=<Real*8>', 'Threshold for spurious eigenvalues', darg) ) then
        
        prm%spu_thres = darg
        
     end if
     if ( get_logical(argv, '-sigval', '-sigval', 'Show singular values') ) then
        
        sigval = .true.
        
     end if
     if ( get_logical(argv, '-inner', '-inner', 'Force to do inner iteration') ) then
        
        prm%force_inner = .true.
        
     end if
     
     !!!!! Add new opiton here
     
     if ( get_logical(argv, '-help', '-help', 'Display this information')) then
        call show_optlist()
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if
     
     if ( .not. option_found ) then
        write(*,*) 'mmzpares Error : unrecognized option "' // trim(argv) // '"'
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if
        
  end do

#ifndef MUMPS
  dense = .true.
#endif
     
  if ( A_given ) then
     if ( rank == root ) then
        open(fileunit, file=trim(filepathA), status='old')
        call mminfo(fileunit, rep, field, symm, rowsA, cols, entries)
        if ( rowsA /= cols ) then
           write(*,*) 'mmzpares Error : Input matrix A must be square'
           error = .true.
        end if
        if ( ( .not. dense ) .and. (symm == 'symmetric' .or. symm == 'hermitian') .and. ( .not. hhpd ) ) then
           nnzmax = 2*entries
        else
           nnzmax = entries
        end if
     end if

     if ( rank == root .and. .not. error ) then
        if ( ( .not. dense ) .and. ( .not. (field == 'real' .and. symm == 'symmetric') ) .and. ( realmat .and. hhpd ) ) then
           write(*,*) 'mmzpares Error : Matrix Market file of A must be real symmetric type ', &
                'if -sspd option or -realmat and -hhpd option are specifined.'
           error = .true.
        end if
        if ( ( .not. dense ) .and. ( .not. (field == 'complex' .and. symm == 'hermitian') ) &
             .and. ( ( .not. realmat ) .and. hhpd) ) then
           write(*,*) 'mmzpares Error : Matrix Market file of A must be complex hermitian type if -hhpd option is specifined.'
           error = .true.
        end if
     end if

#ifdef MPI
     call MPI_BCAST(error, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     if ( error ) then
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if

#ifdef MPI
     call MPI_BCAST(nnzmax, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rowsA, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(field, 7, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(symm, 19, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(entries, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
#endif

     allocate(indx(nnzmax), jndx(nnzmax), rval(nnzmax), cval(nnzmax), ival(nnzmax))

     if ( rank == root ) then
        call mmread(fileunit, rep, field, symm, rowsA, cols, entries, entries, indx, jndx, ival, rval, cval)
        close(fileunit)
     end if

#ifdef MPI
     call MPI_BCAST(indx, entries, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(jndx, entries, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rval, entries, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(cval, entries, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
#endif
     
     call set_matrix(single, realmat, dense, hhpd, rowsA, field, symm, entries, nnzA, indx, jndx, rval, cval &
          ,rowptrA, colindA, svalA, cvalA, dvalA, zvalA, AS, AC, AD, AZ)
     deallocate(indx, jndx, rval, cval, ival)

  else

     if ( rank == root ) then
        write(*,*) 'mmzpares :'
        write(*,*) 'You should specify Matrix Market file of A'
        write(*,*) 'Use option "-fileA=<file>"'
        write(*,*) 'Use option "-help" for details'
        error = .true.
     end if

#ifdef MPI
     call MPI_BCAST(error, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     if ( error ) then
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if
  end if
     
  if ( B_given ) then
     if ( rank == root ) then
        open(fileunit, file=trim(filepathB), status='old')
        call mminfo(fileunit, rep, field, symm, rowsB, cols, entries)
        if ( rowsB /= cols ) then
           write(*,*) 'mmzpares Error : Input matrix B must be square'
           error = .true.
        end if
        if ( rowsA /= rowsB ) then
           write(*,*) 'mmzpares Error : Size of A and B must agree'
           error = .true.
        end if
        if ( ( .not. dense ) .and. (symm == 'symmetric' .or. symm == 'hermitian') .and. ( .not. hhpd ) ) then         
           nnzmax = 2*entries
        else
           nnzmax = entries
        end if
     end if

     if (  rank == root .and. .not. error ) then
        if ( ( .not. dense ) .and. ( .not. (field == 'real' .and. symm == 'symmetric') ) .and. ( realmat .and. hhpd ) ) then
           write(*,*) 'mmzpares Error : Matrix Market file of B must be real symmetric type ', &
                'if -sspd option or -realmat and -hhpd option are specifined.'
           error = .true.
        end if
        if ( ( .not. dense ) .and. ( .not. (field == 'complex' .and. symm == 'hermitian') ) &
             .and. ( ( .not. realmat ) .and. hhpd) ) then
           write(*,*) 'mmzpares Error : Matrix Market file of B must be complex hermitian type if -hhpd option is specifined.'
           error = .true.
        end if
     end if
     
#ifdef MPI
     call MPI_BCAST(error, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     if ( error ) then
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if
        
#ifdef MPI
     call MPI_BCAST(nnzmax, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rowsB, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(field, 7, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(symm, 19, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(entries, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
#endif

     allocate(indx(nnzmax), jndx(nnzmax), rval(nnzmax), cval(nnzmax), ival(nnzmax))
     
     if ( rank == root ) then
        call mmread(fileunit, rep, field, symm, rowsB, cols, entries, entries, indx, jndx, ival, rval, cval)
        close(fileunit)
     end if

#ifdef MPI
     call MPI_BCAST(indx, entries, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(jndx, entries, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rval, entries, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(cval, entries, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
#endif

     call set_matrix(single, realmat, dense, hhpd, rowsA, field, symm, entries, nnzB, indx, jndx, rval, cval &
          ,rowptrB, colindB, svalB, cvalB, dvalB, zvalB, BS, BC, BD, BZ)
#ifdef MPI
     call MPI_BCAST(error, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     if ( error ) then
#ifdef MPI
        call MPI_FINALIZE(ierr)
#endif
        stop
     end if

     deallocate(indx, jndx, rval, cval, ival)

  else

     standard = .true.

  end if
  
#ifdef MPI
  prm%high_comm = MPI_COMM_WORLD
#endif
  
  mat_size = rowsA
  left = cmplx(rleft,ileft,kind(0d0))     
  L = prm%L
  M = prm%M
  Lmax = prm%Lmax

  ncv = zpares_get_ncv(prm)

  if ( single ) then
     if ( realmat ) then
        allocate(XS(mat_size, ncv))
     else
        allocate(XC(mat_size, ncv))
     end if
     allocate(eigvals(ncv), eigvalc(ncv), ress(ncv))
  else
     if ( realmat ) then
        allocate(XD(mat_size, ncv))
     else
        allocate(XZ(mat_size, ncv))
     end if
     allocate(eigvald(ncv), eigvalz(ncv), resd(ncv))
  end if

  if ( dense .and. hhpd .and. UPLOtest ) then
     if ( single ) then
        if ( realmat ) then
           call s_zero_UPLO(mat_size, UPLO, AS)
           call s_zero_UPLO(mat_size, UPLO, BS)
        else
           call c_zero_UPLO(mat_size, UPLO, AC)
           call c_zero_UPLO(mat_size, UPLO, BC)
        end if
     else
        if ( realmat ) then
           call d_zero_UPLO(mat_size, UPLO, AD)
           call d_zero_UPLO(mat_size, UPLO, BD)
        else
           call z_zero_UPLO(mat_size, UPLO, AZ)
           CALL z_zero_UPLO(mat_size, UPLO, BZ)
        end if
     end if
  end if

#ifdef MPI
  total_time = mpi_wtime()
#else
  call system_clock(itc)
  total_time = dble(itc)
#endif

  if ( rank == root ) prm%verbose = 1

  if ( dense ) then
     if ( realmat ) then
        if ( hhpd ) then
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_sdnssyev ...'
                 call zpares_sdnssyev(prm, UPLO, mat_size, AS, mat_size, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XS, ress, info)        
              else
                 if ( rank == root ) write(*,*) 'start zpares_sdnssygv ...'
                 call zpares_sdnssygv(prm, UPLO, mat_size, AS, mat_size, BS, mat_size, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XS, ress, info)        
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_ddnssyev ...'
                 call zpares_ddnssyev(prm, UPLO, mat_size, AD, mat_size, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XD, resd, info)        
              else
                 if ( rank == root ) write(*,*) 'start zpares_ddnssygv ...'
                 call zpares_ddnssygv(prm, UPLO, mat_size, AD, mat_size, BD, mat_size, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XD, resd, info)        
              end if
           end if
        else
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_sdnsgeev ...'
                 call zpares_sdnsgeev(prm, mat_size, AS, mat_size, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XS, ress, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_sdnsgegv ...'
                 call zpares_sdnsgegv(prm, mat_size, AS, mat_size, BS, mat_size, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XS, ress, info)
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_ddnsgeev ...'
                 call zpares_ddnsgeev(prm, mat_size, AD, mat_size, &
                      left, right, num_ev, eigvalz, XD, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_ddnsgegv ...'
                 call zpares_ddnsgegv(prm, mat_size, AD, mat_size, BD, mat_size, &
                      left, right, num_ev, eigvalz, XD, resd, info)                 
              end if
           end if
        end if
     else
        if ( hhpd ) then
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_cdnsheev ...'
                 call zpares_cdnsheev(prm, UPLO, mat_size, AC, mat_size, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XC, ress, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_cdnshegv ...'
                 call zpares_cdnshegv(prm, UPLO, mat_size, AC, mat_size, BC, mat_size, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XC, ress, info)
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_zdnsheev ...'
                 call zpares_zdnsheev(prm, UPLO, mat_size, AZ, mat_size, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XZ, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_zdnshegv ...'
                 call zpares_zdnshegv(prm, UPLO, mat_size, AZ, mat_size, BZ, mat_size, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XZ, resd, info)
              end if
           end if
        else
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_cdnsgeev ...'
                 call zpares_cdnsgeev(prm, mat_size, AC, mat_size, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XC, ress, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_cdnsgegv ...'
                 call zpares_cdnsgegv(prm, mat_size, AC, mat_size, BC, mat_size, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XC, ress, info)
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_zdnsgeev ...'
                 call zpares_zdnsgeev(prm, mat_size, AZ, mat_size, &
                      left, right, num_ev, eigvalz, XZ, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_zdnsgegv ...'
                 call zpares_zdnsgegv(prm, mat_size, AZ, mat_size, BZ, mat_size, &
                      left, right, num_ev, eigvalz, XZ, resd, info)
              end if
           end if
        end if
     end if
  else
#ifdef MUMPS
     if ( realmat ) then
        if ( hhpd ) then
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_smpssyev ...'
                 call zpares_smpssyev(prm, mat_size, rowptrA, colindA, svalA, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XS, ress, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_smpssygv ...'
                 call zpares_smpssygv(prm, mat_size, rowptrA, colindA, svalA, rowptrB, colindB, svalB, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XS, ress, info)
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_dmpssyev ...'
                 call zpares_dmpssyev(prm, mat_size, rowptrA, colindA, dvalA, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XD, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_dmpssygv ...'
                 call zpares_dmpssygv(prm, mat_size, rowptrA, colindA, dvalA, rowptrB, colindB, dvalB, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XD, resd, info)
              end if
           end if
        else
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_smpsgeev ...'
                 call zpares_smpsgeev(prm, mat_size, rowptrA, colindA, svalA, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XS, ress, info)      
              else
                 if ( rank == root ) write(*,*) 'start zpares_smpsgegv ...'
                 call zpares_smpsgegv(prm, mat_size, rowptrA, colindA, svalA, rowptrB, colindB, svalB, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XS, ress, info)      
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_dmpsgeev ...'
                 call zpares_dmpsgeev(prm, mat_size, rowptrA, colindA, dvalA, &
                      left, right, num_ev, eigvalz, XD, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_dmpsgegv ...'
                 call zpares_dmpsgegv(prm, mat_size, rowptrA, colindA, dvalA, rowptrB, colindB, dvalB, &
                      left, right, num_ev, eigvalz, XD, resd, info)
              end if
           end if
        end if
     else
        if ( hhpd ) then
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_cmpsheev ...'
                 call zpares_cmpsheev(prm, mat_size, rowptrA, colindA, cvalA, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XC, ress, info)        
              else
                 if ( rank == root ) write(*,*) 'start zpares_cmpshegv ...'
                 call zpares_cmpshegv(prm, mat_size, rowptrA, colindA, cvalA, rowptrB, colindB, cvalB, &
                      real(left,kind(0.0)), real(right,kind(0.0)), num_ev, eigvals, XC, ress, info)        
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_zmpsheev ...'
                 call zpares_zmpsheev(prm, mat_size, rowptrA, colindA, zvalA, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XZ, resd, info)        
              else
                 if ( rank == root ) write(*,*) 'start zpares_zmpshegv ...'
                 call zpares_zmpshegv(prm, mat_size, rowptrA, colindA, zvalA, rowptrB, colindB, zvalB, &
                      real(left,kind(0d0)), right, num_ev, eigvald, XZ, resd, info)        
              end if
           end if
        else
           if ( single ) then
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_cmpsgeev ...'
                 call zpares_cmpsgeev(prm, mat_size, rowptrA, colindA, cvalA, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XC, ress, info)        
              else
                 if ( rank == root ) write(*,*) 'start zpares_cmpsgegv ...'
                 call zpares_cmpsgegv(prm, mat_size, rowptrA, colindA, cvalA, rowptrB, colindB, cvalB, &
                      cmplx(dble(left),aimag(left),kind(0.0)), real(right,kind(0.0)), num_ev, eigvalc, XC, ress, info)        
              end if
           else
              if ( standard ) then
                 if ( rank == root ) write(*,*) 'start zpares_zmpsgeev ...'
                 call zpares_zmpsgeev(prm, mat_size, rowptrA, colindA, zvalA, &
                      left, right, num_ev, eigvalz, XZ, resd, info)
              else
                 if ( rank == root ) write(*,*) 'start zpares_zmpsgegv ...'
                 call zpares_zmpsgegv(prm, mat_size, rowptrA, colindA, zvalA, rowptrB, colindB, zvalB, &
                      left, right, num_ev, eigvalz, XZ, resd, info)
              end if
           end if
        end if
     end if
#endif
  end if

#ifdef MPI
  if ( rank == root ) then
     write(*,'(1x,a,I8)') 'Number of MPI processes = ', nprocs
  end if
  total_time = mpi_wtime() - total_time !!!!
#else
  call system_clock(itc, t_rate, t_max)
  dtc = dble(itc)
  if ( dtc < total_time ) then
     total_time = dble(t_max) - total_time + dtc
  else
     total_time = dtc - total_time
  end if
  total_time = total_time / dble(t_rate)
#endif

  if ( rank == root ) then
     write(*,*)
     write(*,*) 'Parameters ----------------------------------------'
     sign = '-'
     if ( aimag(left) >= 0 ) sign = '+'
     if ( hhpd ) then
        write(*,'(1x,a,1x,1pe10.3)') 'emin = ', dble(left)
        write(*,'(1x,a,1pe10.3)') 'emax = ', right
     else
        write(*,'(1x,a,1x,1pe10.3,1x,a1,1x,1pe10.3," i")') 'left = ', dble(left), sign, abs(aimag(left))
        write(*,'(1x,a,1pe10.3)') 'right = ', right
     end if
     write(*,'(1x,a,I5)') 'N = ', prm%N
     write(*,'(1x,a,I5)') 'M = ', prm%M
     write(*,'(1x,a,I5)') 'initial L = ', L
     write(*,'(1x,a,I5)') 'L = ', prm%L
     write(*,'(1x,a,I5)') 'L_max = ', prm%Lmax
     write(*,'(1x,a,I5)') 'imax = ', prm%imax
     if ( prm%extract == ZPARES_EXTRACT_RR ) then
        write(*,'(1x,a,a)') 'extract = Rayleigh-Ritz'
     else
        write(*,'(1x,a,a)') 'extract = Explicit moment'
     end if
     write(*,'(1x,a,1pe10.3)') 'delta = ', prm%delta
     write(*,'(1x,a,1pe10.3)') 'asp_ratio = ', prm%asp_ratio
     write(*,'(1x,a,1pe10.3)') 'tol = ', prm%tol
     
     write(*,*)
     write(*,*) 'Results    ----------------------------------------'
     write(*,'(1x,a,I5)') 'info = ', info
     write(*,'(1x,a,I5)') 'num_basis = ', prm%num_basis
     write(*,'(1x,a,I5)') 'iter = ', prm%iter
     do i = 1, prm%iter
        write(*,'(1x,a,I2,a,I2)') 'iter_info(', i, ') = ', prm%iter_info(i)
     end do
     write(*,'(1x,a,I10)') 'Estimated eigenvalue count', prm%num_ev_est
     write(*,'(1x,a,I10)') 'Number of approximated eigenpairs', num_ev
     
     if ( hhpd ) then
        write(*,'(a6,3x,"|",1x,a10,15x,"|",1x,a8,"|",1x,a9)') 'index','eigenvalue','residual','spu. ind.'
        do i = 1, num_ev
           if ( single ) then
              write(*,'(I6,4x,1pe23.16,4x,1pe10.3,4x,1pe10.3)') i, eigvals(i), ress(i), prm%indi_spu(i)
           else
              write(*,'(I6,4x,1pe23.16,4x,1pe10.3,4x,1pe10.3)') i, eigvald(i), resd(i), prm%indi_spu(i)
           end if
        end do
     else
        write(*,'(a6,3x,"|",1x,a10,43x,"|",1x,a8,"|",1x,a9)') 'index','eigenvalue','residual','spu. ind.' 
        do i = 1, num_ev
           sign = '-'
           if ( single ) then
              if ( aimag(eigvalc(i)) >= 0 ) sign = '+'
              write(*,'(I6,4x,1pe14.7,1x,a1,1x,1pe14.7," i",4x,1pe10.3,4x,1pe10.3)') &
                   i, dble(eigvalc(i)), sign, abs(aimag(eigvalc(i))), ress(i), prm%indi_spu(i)
           else
              if ( aimag(eigvalz(i)) >= 0 ) sign = '+'
              write(*,'(I6,4x,1pe23.16,1x,a1,1x,1pe23.16," i",4x,1pe10.3,4x,1pe10.3)') &
                   i, dble(eigvalz(i)), sign, abs(aimag(eigvalz(i))), resd(i), prm%indi_spu(i)
           end if
        end do
     end if
     if ( sigval ) then
        write(*,*) 'Singular values'
        write(*,'(a6,3x,"|",1x,a10,15x,"|",1x,a10)') 'index','absolute','relative'
        do i = 1, prm%L*prm%M
           write(*,'(I6,4x,1pe23.16,4x,1pe23.16)') i, prm%sig_val(i), prm%sig_val(i) / prm%sig_val(1)
        end do
     end if

     write(*,*) ' Timing data               [sec] -------------------'
     write(*,'(1x,a,1x,f10.3)') 'Total               ', total_time
     write(*,'(1x,a,1x,f10.3)') 'Generate rand       ', prm%timer_rand%t
     write(*,'(1x,a,1x,f10.3)') 'Factorization       ', prm%timer_fact%t
     write(*,'(1x,a,1x,f10.3)') 'Solve linear system ', prm%timer_solve%t
     if ( .not. realmat .and. hhpd ) then
        write(*,'(1x,a,1x,f10.3)') 'Fact (z*B-A)^H      ', prm%timer_fact_H%t
        write(*,'(1x,a,1x,f10.3)') 'Solve (z*B-A)^H     ', prm%timer_solve_H%t
     end if
     write(*,'(1x,a,1x,f10.3)') 'Integral sum        ', prm%timer_sum%t
     write(*,'(1x,a,1x,f10.3)') 'Orthonormalizaiton  ', prm%timer_orth%t
     write(*,'(1x,a,1x,f10.3)') 'SVD                 ', prm%timer_serial_SVD%t
     write(*,'(1x,a,1x,f10.3)') 'Diag reduced problem', prm%timer_reduced_eig%t
     write(*,'(1x,a,1x,f10.3)') 'Mult A              ', prm%timer_mult_A%t
     write(*,'(1x,a,1x,f10.3)') 'Mult B              ', prm%timer_mult_B%t
     write(*,'(1x,a,1x,f10.3)') 'Gemm reduce         ', prm%timer_gemm_reduce%t
     write(*,'(1x,a,1x,f10.3)') 'Rotation            ', prm%timer_sub_rot%t
     write(*,'(1x,a,1x,f10.3)') 'Calc residual norm  ', prm%timer_res_norm%t
     
  end if

  call zpares_finalize(prm)

#ifdef MPI
  call MPI_FINALIZE(ierr)
#endif

  contains
    
    subroutine append_to_optlist(id, usage, desc)
      implicit none
      character(*), intent(in) :: id, usage, desc
      type(option), pointer :: newopt=>NULL(), p=>NULL(), q=>NULL()
      
      p => optlist
      do while(associated(p))
         if ( trim(p%id) == trim(id) ) then
            return
         end if
         q => p
         p => p%next
      end do
      
      allocate(newopt)
      newopt%id = id
      newopt%usage = usage
      newopt%description = desc
      nullify(newopt%next)
      q%next => newopt
    end subroutine append_to_optlist

    subroutine show_optlist()
      implicit none
      type(option), pointer :: p
      
      write(*,'(a)') 'Options:'
      p => optlist%next
      do while(associated(p))
         write(*,'(2x,a,a)') p%usage, trim(p%description)
         p => p%next
      end do
    end subroutine show_optlist

    logical function get_int(argv, id, usage, desc, ret)
      implicit none
      character(*), intent(in) :: argv, id, usage, desc
      integer, intent(out) :: ret
      character :: tmp*512
      integer :: idx

      call append_to_optlist(id, usage, desc)
      idx = index(argv,id // '=')
      if ( idx == 1 ) then
         option_found = .true.
         tmp = argv(len(id // '=')+1:len(argv))
         read(tmp,*) ret
      end if
      get_int = idx == 1
    end function get_int
    
    logical function get_double(argv, id, usage, desc, ret)
      implicit none
      character(*), intent(in) :: argv, id, usage, desc
      double precision, intent(out) :: ret
      character :: tmp*512
      integer :: idx

      call append_to_optlist(id, usage, desc)
      idx = index(argv,id // '=')
      if ( idx == 1 ) then
         option_found = .true.
         tmp = argv(len(id // '=')+1:len(argv))
         read(tmp,*) ret
      end if
      get_double = idx == 1
    end function get_double
    
    logical function get_str(argv, id, usage, desc, ret)
      implicit none
      character(*), intent(in) :: argv, id, usage, desc
      character(*), intent(out) :: ret
      character :: tmp*512
      integer :: idx
      
      call append_to_optlist(id, usage, desc)
      idx = index(argv,id // '=')
      if ( idx == 1 ) then
         option_found = .true.
         ret = trim(argv(len(id // '=')+1:len(argv)))
      end if
      get_str = idx == 1
    end function get_str
    
    logical function get_logical(argv, id, usage, desc)
      implicit none
      character(*), intent(in) :: argv, id, usage, desc
      
      call append_to_optlist(id, usage, desc)
      get_logical = index(argv,id) == 1
      if ( get_logical ) option_found = .true.
    end function get_logical

    subroutine make_half2full(symm, entries, indx, jndx, rval, cval)
      implicit none
      character(*), intent(in) :: symm
      integer, intent(in) :: entries
      integer, intent(inout) :: indx(:), jndx(:)
      double precision, intent(inout) :: rval(:)
      complex(kind(0d0)), intent(inout) :: cval(:)
      
      integer :: i, j, count
      count = 0
      do i = 1, entries
         if ( indx(i) /= jndx(i) ) then
            count = count + 1
            indx(entries+count) = jndx(i)
            jndx(entries+count) = indx(i)
            rval(entries+count) = rval(i)
            cval(entries+count) = conjg(cval(i))
         end if
      end do
    end subroutine make_half2full

    subroutine set_matrix(single, realmat, dense, hhpd, mat_size, field, symm, entries, nnz, indx, jndx, rcoo, ccoo &
         , rowptr, colind, sval, cval, dval, zval, sF, cF, dF, zF)
      implicit none
      logical, intent(in) :: single, realmat, dense, hhpd
      character(*), intent(in) :: field, symm
      integer, intent(in) :: mat_size, entries
      integer, intent(inout) :: indx(:), jndx(:)
      integer, intent(out) :: nnz
      integer, pointer :: rowptr(:), colind(:)
      real, pointer :: sval(:), sF(:,:)
      complex, pointer :: cval(:), cF(:,:)
      double precision, intent(inout) :: rcoo(:)
      double precision, pointer :: dval(:), dF(:,:)
      complex(kind(0d0)), intent(inout) :: ccoo(:)
      complex(kind(0d0)), pointer :: zval(:), zF(:,:)
      integer :: k

      if ( ( .not. dense ) .and. (symm == 'symmetric' .or. symm == 'hermitian') .and. ( .not. hhpd) ) then
         nnz = 2*entries
         do k = 1, entries
            if ( indx(k) == jndx(k) ) then
               nnz = nnz - 1
            end if
         end do
         call make_half2full(symm, entries, indx, jndx, rcoo, ccoo)
      else
         nnz = entries
      end if
      if ( .not. dense ) then
         allocate(rowptr(mat_size+1), colind(nnz))
      end if

      if ( single ) then
         if ( realmat ) then
            if ( dense ) then
               allocate(sF(mat_size,mat_size))
               sF = 0.0
            else
               allocate(sval(nnz))
               call s_coo2csr(mat_size, nnz, indx, jndx, rcoo, ccoo, rowptr, colind, sval)
            end if
         else
            if ( dense ) then
               allocate(cF(mat_size,mat_size))
               cF = (0.0,0.0)
            else
               allocate(cval(nnz))
               call c_coo2csr(mat_size, nnz, indx, jndx, rcoo, ccoo, rowptr, colind, cval)
            end if
         end if
      else
         if ( realmat ) then
            if ( dense ) then
               allocate(dF(mat_size,mat_size))
               dF = 0d0
            else
               allocate(dval(nnz))
               call d_coo2csr(mat_size, nnz, indx, jndx, rcoo, ccoo, rowptr, colind, dval)
            end if
         else
            if ( dense ) then
               allocate(zF(mat_size,mat_size))
               zF = (0d0,0d0)
            else
               allocate(zval(nnz))
               call z_coo2csr(mat_size, nnz, indx, jndx, rcoo, ccoo, rowptr, colind, zval)
            end if
         end if
      end if
      if ( dense ) then
         if ( field == 'complex' ) then
            if ( single ) then
               if ( realmat ) then
                  do j = 1, nnz
                     sF(indx(j),jndx(j)) = ccoo(j)
                  end do
               else
                  do j = 1, nnz
                     cF(indx(j),jndx(j)) = ccoo(j)
                  end do
               end if
            else
               if ( realmat ) then
                  do j = 1, nnz
                     dF(indx(j),jndx(j)) = ccoo(j)
                  end do
               else
                  do j = 1, nnz
                     zF(indx(j),jndx(j)) = ccoo(j)
                  end do
               end if
            end if
         else
            if ( single ) then
               if ( realmat ) then
                  do j = 1, nnz
                     sF(indx(j),jndx(j)) = cmplx(rcoo(j),0.0,kind(0.0))
                  end do
               else
                  do j = 1, nnz
                     cF(indx(j),jndx(j)) = cmplx(rcoo(j),0.0,kind(0.0))           
                  end do
               end if
            else
               if ( realmat ) then
                  do j = 1, nnz
                     dF(indx(j),jndx(j)) = cmplx(rcoo(j),0d0,kind(0d0))           
                  end do
               else
                  do j = 1, nnz
                     zF(indx(j),jndx(j)) = cmplx(rcoo(j),0d0,kind(0d0))           
                  end do
               end if
            end if
         end if
         if ( symm == 'symmetric' ) then
            if ( single ) then
               if ( realmat ) then
                  do j = 1, nnz
                     sF(jndx(j),indx(j)) = sF(indx(j),jndx(j))
                  end do
               else
                  do j = 1, nnz
                     cF(jndx(j),indx(j)) = cF(indx(j),jndx(j))
                  end do
               end if
            else
               if ( realmat ) then
                  do j = 1, nnz
                     dF(jndx(j),indx(j)) = dF(indx(j),jndx(j))
                  end do
               else
                  do j = 1, nnz
                     zF(jndx(j),indx(j)) = zF(indx(j),jndx(j))
                  end do
               end if               
            end if
         else if ( symm == 'hermitian' ) then
            if ( single ) then
               if ( realmat ) then
                  do j = 1, nnz
                     sF(jndx(j),indx(j)) = sF(indx(j),jndx(j))
                  end do
               else
                  do j = 1, nnz
                     cF(jndx(j),indx(j)) = conjg(cF(indx(j),jndx(j)))
                  end do
               end if
            else
               if ( realmat ) then
                  do j = 1, nnz
                     dF(jndx(j),indx(j)) = dF(indx(j),jndx(j))
                  end do
               else
                  do j = 1, nnz
                     zF(jndx(j),indx(j)) = conjg(zF(indx(j),jndx(j)))
                  end do
               end if
            end if
         end if
      end if
    end subroutine set_matrix
    
! include double complex routines
#include "mmzpares_inc.f90" 
#define REALMAT
! include double precision routines
#include "mmzpares_inc.f90" 
#define SINGLE
! include real routines
#include "mmzpares_inc.f90" 
#undef REALMAT
! include complex routines
#include "mmzpares_inc.f90" 

  end program main

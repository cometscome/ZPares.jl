module zpares
  implicit none

  integer, parameter :: ZPARES_TASK_NONE = -1
  integer, parameter :: ZPARES_TASK_FACTO = 1
  integer, parameter :: ZPARES_TASK_FACTO_H = 2
  integer, parameter :: ZPARES_TASK_SOLVE = 3
  integer, parameter :: ZPARES_TASK_SOLVE_H = 4
  integer, parameter :: ZPARES_TASK_MULT_A = 5
  integer, parameter :: ZPARES_TASK_MULT_B = 6
  integer, parameter :: ZPARES_TASK_FINISH = 0
  
  integer, parameter :: ZPARES_STATE_INIT1 = 0
  integer, parameter :: ZPARES_STATE_INIT2 = 1
  integer, parameter :: ZPARES_STATE_LINSOLVE = 2
  integer, parameter :: ZPARES_STATE_ALLREDUCE = 3
  integer, parameter :: ZPARES_STATE_STC = 4
  integer, parameter :: ZPARES_STATE_LADD = 5
  integer, parameter :: ZPARES_STATE_SS_SVD = 6
  integer, parameter :: ZPARES_STATE_INIT_REDUCED_EIG = 7
  integer, parameter :: ZPARES_STATE_FORM_REDUCED_EIG = 8
  integer, parameter :: ZPARES_STATE_SOLVE_REDUCED_EIG = 9
  integer, parameter :: ZPARES_STATE_CALC_RES0 = 10
  integer, parameter :: ZPARES_STATE_CALC_RES1 = 11
  integer, parameter :: ZPARES_STATE_CALC_RES2 = 12
  integer, parameter :: ZPARES_STATE_CHECK_OUTER = 13
  integer, parameter :: ZPARES_STATE_FINISH = 14
  integer, parameter :: ZPARES_STATE_FIN_REDUCED_EIG = 15
  
  integer, parameter :: ZPARES_QUAD_USER = 0
  integer, parameter :: ZPARES_QUAD_ELL_TRAP = 1
  
  integer, parameter :: ZPARES_QU_MODE_RULE = 0
  integer, parameter :: ZPARES_QU_MODE_BACK = 1
  
  integer, parameter :: ZPARES_EXTRACT_RR = 0
  integer, parameter :: ZPARES_EXTRACT_EM = 1

  integer, parameter :: ZPARES_MODE_STC = 0
  integer, parameter :: ZPARES_MODE_LADD = 1
  integer, parameter :: ZPARES_MODE_SS = 2

  integer, parameter :: ZPARES_INFO_UNDETERMINED = -1
  integer, parameter :: ZPARES_INFO_SUCCESS = 0
  integer, parameter :: ZPARES_INFO_INVALID_PARAM = 1
  integer, parameter :: ZPARES_INFO_LAPACK_ERROR = 2
  integer, parameter :: ZPARES_INFO_MPI_ERROR = 3
  integer, parameter :: ZPARES_INFO_NOSUF_USER_SRC = 4

  type zpares_timer
     double precision :: start_time
     double precision :: t = 0d0
  end type zpares_timer

  ! this derived type is only used in this module
  type pointers
     real, pointer :: sprojAB(:,:) => NULL(), sMu(:,:) => NULL()
     complex, pointer :: cprojAB(:,:) => NULL(), cMu(:,:) => NULL()
     double precision, pointer :: dprojAB(:,:) => NULL(), dMu(:,:) => NULL()
     complex(kind(0d0)), pointer :: zprojAB(:,:) => NULL(), zMu(:,:) => NULL()
  end type pointers

  !> This derived type is used to store optional parameters and context for zpares routines
  type zpares_prm
     character(8) :: version = '0.9.6a' !> Version number

     ! input
     integer :: N = 32 !< Number of integral points
     integer :: M = 16 !< Degree of moment
     integer :: Lstep = 8
     integer :: Lmax = 64 !< Maximum number of columns of the source matrix
     integer :: quad_type = ZPARES_QUAD_ELL_TRAP !< Type of quadrature rule
     integer :: extract = ZPARES_EXTRACT_RR !< Extract method
     integer :: imax = 0 !< Maximum number of iteration
     integer :: n_orth = 3 !< Number of itaration for orthonormalization
     logical :: calc_res = .true. !< Flag for calculating residuals
     logical :: user_source = .false. !< User specified source matrix is given
     logical :: Hermitian = .false. !< Input matrices A and B are Hermitian
     logical :: B_pos_def = .true. !< Input matrix B is positive definite
     logical :: standard = .false. 
     logical :: trim_out = .true. !< Discard eigenvalues located outside of the ellipse
     logical :: trim_res = .false. 
     logical :: trim_spu = .true. !< Discard spurious eigenvalues
     logical :: user_path = .false.
     logical :: sym_contour = .false.
     logical :: force_inner = .false.
     double precision :: delta = 1d-12 !< Threshold value for numerical rank
     double precision :: asp_ratio = 1d0 !< Aspect ratin of the ellipse
     double precision :: tol = 1d-14 !< Tolerance for residual
     double precision :: spu_thres = 1d-6 !< Thresfold for trimming away the spirious eigenpairs

     integer :: low_comm !< Lower level communicator
     integer :: high_comm !< Higher level communicator   
     integer :: write_unit = 6 
     integer :: verbose = 0 !< Verbose level

     ! output
     integer :: itask = ZPARES_TASK_NONE !< Reverse communication task
     integer :: num_ev_est !< Estimation of eigenvalue count
     integer :: ws !< starting columns of work arrays
     integer :: xs !< starting columns of X array
     integer :: nc !< number of columns for solve or matvec
     integer :: quad_idx !< Index of quadrature points
     integer :: iter = 0 !< Number of iteration
     integer :: num_basis !< Number of basis vectors
     integer, pointer :: iter_info(:) => NULL()
     double precision, pointer :: sig_val(:) => NULL()
     double precision, pointer :: indi_spu(:) => NULL()

     integer :: imisc(64)
     double precision :: dmisc(64)

     type(zpares_timer) :: timer_rand
     type(zpares_timer) :: timer_fact
     type(zpares_timer) :: timer_solve
     type(zpares_timer) :: timer_fact_H
     type(zpares_timer) :: timer_solve_H
     type(zpares_timer) :: timer_orth
     type(zpares_timer) :: timer_serial_SVD
     type(zpares_timer) :: timer_reduced_eig
     type(zpares_timer) :: timer_sub_rot
     type(zpares_timer) :: timer_mult_A
     type(zpares_timer) :: timer_mult_B
     type(zpares_timer) :: timer_res_norm
     type(zpares_timer) :: timer_gemm_reduce
     type(zpares_timer) :: timer_sum

     ! input/output
     integer :: L = 16 !< Number of columns of the source matrix

     ! private
     integer :: state = ZPARES_STATE_INIT1 !< User should not touch this variable
     integer :: mode = ZPARES_MODE_STC !< User should not touch this variable
     integer :: x_offset !< User should not touch this variable
     integer :: x_ncol !< User should not touch this variable

     type(pointers) :: ptrs !< User should not touch this variable
  end type zpares_prm

  character(len=15), private, parameter  :: errpfx = 'zpares error : '

  ! integer, private, parameter :: INFO_SUCCESS = 0
  ! integer, private, parameter :: INFO_UNDETERMINED = -99
  ! integer, private, parameter :: INFO_INVALID_PARAM = 1
  ! integer, private, parameter :: INFO_NOT_ENOUGH_USER_SOURCE = 2

  private :: check_inputs, para_range, get_rank_and_size, start_timer, stop_timer

  interface zpares_rci_sub_get_projAB
     module procedure zpares_srci_sub_get_projAB
     module procedure zpares_crci_sub_get_projAB
     module procedure zpares_drci_sub_get_projAB
     module procedure zpares_zrci_sub_get_projAB
  end interface zpares_rci_sub_get_projAB

  interface zpares_rci_sub_get_Mu
     module procedure zpares_srci_sub_get_Mu
     module procedure zpares_crci_sub_get_Mu
     module procedure zpares_drci_sub_get_Mu
     module procedure zpares_zrci_sub_get_Mu
  end interface zpares_rci_sub_get_Mu

  interface zpares_rci_sub_linsolve
     module procedure zpares_srci_sub_linsolve
     module procedure zpares_crci_sub_linsolve
     module procedure zpares_drci_sub_linsolve
     module procedure zpares_zrci_sub_linsolve
  end interface zpares_rci_sub_linsolve

  interface zpares_rci_sub_svd
     module procedure zpares_srci_sub_svd
     module procedure zpares_crci_sub_svd
     module procedure zpares_drci_sub_svd
     module procedure zpares_zrci_sub_svd
  end interface zpares_rci_sub_svd

  interface zpares_rci_sub_Rayleigh_Ritz
     module procedure zpares_srci_sub_Rayleigh_Ritz
     module procedure zpares_crci_sub_Rayleigh_Ritz
     module procedure zpares_drci_sub_Rayleigh_Ritz
     module procedure zpares_zrci_sub_Rayleigh_Ritz
  end interface zpares_rci_sub_Rayleigh_Ritz

  interface zpares_rci_sub_Hankel_method
     module procedure zpares_srci_sub_Hankel_method
     module procedure zpares_crci_sub_Hankel_method
     module procedure zpares_drci_sub_Hankel_method
     module procedure zpares_zrci_sub_Hankel_method
  end interface zpares_rci_sub_Hankel_method

  interface zpares_rci_sub_calc_res
     module procedure zpares_srci_sub_calc_res
     module procedure zpares_crci_sub_calc_res
     module procedure zpares_drci_sub_calc_res
     module procedure zpares_zrci_sub_calc_res
  end interface zpares_rci_sub_calc_res

  interface zpares_rci_sub_get_sym_contour
     module procedure zpares_crci_sub_get_sym_contour
     module procedure zpares_zrci_sub_get_sym_contour
  end interface zpares_rci_sub_get_sym_contour

contains

  subroutine zpares_init(prm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(zpares_prm) :: prm, new_prm
    integer :: LmaxM
    prm = new_prm
#ifdef MPI
    prm%high_comm = MPI_COMM_SELF
    prm%low_comm = MPI_COMM_SELF    
#endif
  end subroutine zpares_init

  integer function zpares_get_ncv(prm)
    implicit none
    type(zpares_prm) :: prm
    zpares_get_ncv = prm%Lmax*prm%M
  end function zpares_get_ncv

  subroutine zpares_finalize(prm)
    implicit none
    type(zpares_prm) :: prm

    if ( associated(prm%sig_val) ) then
       deallocate(prm%sig_val)
    end if
    if ( associated(prm%indi_spu) ) then
       deallocate(prm%indi_spu)
    end if
    if ( associated(prm%iter_info) ) then
       deallocate(prm%iter_info)
    end if

    if ( associated(prm%ptrs%sMu) ) then
       deallocate(prm%ptrs%sMu)
    end if
    if ( associated(prm%ptrs%cMu) ) then
       deallocate(prm%ptrs%cMu)
    end if
    if ( associated(prm%ptrs%dMu) ) then
       deallocate(prm%ptrs%dMu)
    end if
    if ( associated(prm%ptrs%zMu) ) then
       deallocate(prm%ptrs%zMu)
    end if

    if ( associated(prm%ptrs%sprojAB) ) then
       deallocate(prm%ptrs%sprojAB)
    end if
    if ( associated(prm%ptrs%cprojAB) ) then
       deallocate(prm%ptrs%cprojAB)
    end if
    if ( associated(prm%ptrs%dprojAB) ) then
       deallocate(prm%ptrs%dprojAB)
    end if
    if ( associated(prm%ptrs%zprojAB) ) then
       deallocate(prm%ptrs%zprojAB)
    end if    
  end subroutine zpares_finalize

! include double complex routines
#include "zpares_rci_inc.f90" 
#include "zpares_rci_sub_inc.f90" 
#include "zpares_dense_inc.f90" 
#define REALMAT
! include double precision routines
#include "zpares_rci_inc.f90" 
#include "zpares_rci_sub_inc.f90" 
#include "zpares_dense_inc.f90" 
#define SINGLE
! include real routines
#include "zpares_rci_inc.f90" 
#include "zpares_rci_sub_inc.f90" 
#include "zpares_dense_inc.f90" 
#undef REALMAT
! include complex routines
#include "zpares_rci_inc.f90" 
#include "zpares_rci_sub_inc.f90" 
#include "zpares_dense_inc.f90" 
  
  logical function check_inputs(prm)
    type(zpares_prm), intent(inout) :: prm
    integer :: unit
    logical :: v

    unit = prm%write_unit
    check_inputs = .false.
    v = prm%verbose >= 1

    if ( not_posint(prm%N, 'N', unit, v) ) then
       return
    end if

    if ( mod(prm%N, 2) == 1 ) then
       if (v) write(unit,*) errpfx, 'N must be an even number.'
       return
    end if

    if ( not_posint(prm%M, 'M', unit, v) ) then
       return
    end if

    if ( prm%M > prm%N ) then
       if (v) write(unit,*) errpfx, 'M must be less than or equal to the number of quadrature points N.'
       return
    end if

    if ( not_posint(prm%L, 'L', unit, v) ) then
       return
    end if

    if ( prm%delta < 0d0 ) then
       if (v) write(unit,*) errpfx, 'delta must be positive.'
       return
    end if
    
    if ( prm%extract /= 0 .and. prm%extract /= 1 ) then
       if (v) write(unit,*) errpfx, 'extract must be 0 or 1.'
       return
    end if

    if ( prm%asp_ratio <= 0d0) then
       if (v) write(unit,*) errpfx, 'asp_ratio must be strictly positive.'
       return
    end if

    if ( prm%tol < 0d0 ) then
       if (v) write(unit,*) errpfx, 'tol must be positive.'
       return
    end if

    if ( not_posint(prm%Lmax, 'Lmax', unit, v) ) then
       return
    end if

    if ( prm%L > prm%Lmax ) then
       if (v) write(unit,*) errpfx, 'L must be less than or equal to Lmax.'
       return
    end if

    if ( prm%imax < 0 ) then
       if (v) write(unit,*) errpfx, 'imax must be positive.'
       return
    end if

    if ( prm%n_orth < 0 ) then
       if (v) write(unit,*) errpfx, 'n_orth must be positive.'
       return
    end if

    if ( prm%spu_thres < 0d0 ) then
       if (v) write(unit,*) errpfx, 'spu_thres must be positive.'
       return
    end if
    
    check_inputs = .true.

    contains
      logical function not_posint(intvar, desc, unit, verbose_flg)
        implicit none
        integer, intent(in) :: intvar, unit
        character(*) :: desc
        logical, intent(in) :: verbose_flg
        
        if ( intvar <= 0 ) then
           if ( verbose_flg) write(unit,*) errpfx, desc, ' must be positive.'
           not_posint = .true.
        else
           not_posint = .false.
        end if
      end function not_posint      
  end function check_inputs

  subroutine zpares_show_result(prm)
    implicit none
    type(zpares_prm) :: prm        
  end subroutine zpares_show_result
    
  subroutine para_range(global_count, nprocs, rank, start, local_count)
    implicit none
    integer, intent(in) :: global_count, nprocs, rank
    integer, intent(out) :: start, local_count
    start = rank*(global_count/nprocs) + 1 + min(rank, mod(global_count, nprocs))
    local_count = global_count/nprocs
    if ( mod(global_count, nprocs) > rank ) local_count = local_count + 1
  end subroutine para_range

  subroutine get_rank_and_size(mpi_comm, rank, size)
    implicit none
#ifdef MPI
    integer :: ierr
#endif    
    integer, intent(in) :: mpi_comm
    integer, intent(out) :: rank, size
#ifdef MPI
    call MPI_COMM_RANK(mpi_comm, rank, ierr)
    call MPI_COMM_SIZE(mpi_comm, size, ierr)
#else
    rank = 0
    size = 1
#endif
  end subroutine get_rank_and_size

  subroutine start_timer(t)
    implicit none
    type(zpares_timer), intent(inout) :: t
#ifndef MPI
    integer :: itc
    call system_clock(itc)
    t%start_time = dble(itc)
#else
    include 'mpif.h'
    t%start_time = mpi_wtime()
#endif
  end subroutine start_timer

  subroutine stop_timer(t)
    implicit none
    type(zpares_timer), intent(inout) :: t 
#ifdef MPI  
    include 'mpif.h'
    t%t = t%t + mpi_wtime() - t%start_time
#else
    integer :: itc, t_rate, t_max
    double precision :: dtc, diff
    call system_clock(itc, t_rate, t_max)
    dtc = dble(itc)
    if ( dtc < t%start_time ) then
       diff = dble(t_max) - t%start_time + dtc
    else
       diff = dtc - t%start_time
    end if
    diff = diff / dble(t_rate)
    t%t = t%t + diff
#endif
  end subroutine stop_timer

  subroutine clear_timer(t)
    implicit none
    type(zpares_timer), intent(inout) :: t
    t%t = 0d0
  end subroutine clear_timer

end module zpares

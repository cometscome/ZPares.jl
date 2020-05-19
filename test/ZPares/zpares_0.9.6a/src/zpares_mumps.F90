module zpares_mumps

  private ! all interfaces, subroutines and module variables are private by default

  interface val_COO_zB_A
     module procedure s_val_COO_zB_A
     module procedure d_val_COO_zB_A
     module procedure c_val_COO_zB_A
     module procedure z_val_COO_zB_A
  end interface

  interface val_COO_zB_A_h2f
     module procedure c_val_COO_zB_A_h2f
     module procedure z_val_COO_zB_A_h2f
  end interface

  interface val_COO_zI_A
     module procedure s_val_COO_zI_A
     module procedure d_val_COO_zI_A
     module procedure c_val_COO_zI_A
     module procedure z_val_COO_zI_A
  end interface

  interface val_COO_zI_A_h2f
     module procedure c_val_COO_zI_A_h2f
     module procedure z_val_COO_zI_A_h2f
  end interface

  interface matvec_csr
     module procedure s_matvec_csr
     module procedure d_matvec_csr
     module procedure c_matvec_csr
     module procedure z_matvec_csr     
  end interface

  interface matvec_csr_hermite
     module procedure s_matvec_csr_hermite
     module procedure d_matvec_csr_hermite
     module procedure c_matvec_csr_hermite
     module procedure z_matvec_csr_hermite
  end interface

  interface MUMPS_wrap
     module procedure c_MUMPS_wrap
     module procedure z_MUMPS_wrap
  end interface

  public :: zpares_smpsgegv, zpares_cmpsgegv, zpares_dmpsgegv, zpares_zmpsgegv 
  public :: zpares_smpssygv, zpares_cmpshegv, zpares_dmpssygv, zpares_zmpshegv 
  public :: zpares_smpsgeev, zpares_cmpsgeev, zpares_dmpsgeev, zpares_zmpsgeev 
  public :: zpares_smpssyev, zpares_cmpsheev, zpares_dmpssyev, zpares_zmpsheev 

  contains

    integer function count_nnz_zB_A(mat_size, rowptrA, colindA, rowptrB, colindB)
      implicit none
      integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
      integer :: i, j, jA, jB, posA, posB, posAend, posBend
          
      count_nnz_zB_A = 0
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
               posA = posA + 1
               posB = posB + 1
            else if ( jB < jA ) then
               posB = posB + 1
            else if ( jA < jB ) then
               posA = posA + 1
            end if
            count_nnz_zB_A = count_nnz_zB_A + 1
            
            if ( posA > posAend .and. posB > posBend ) then
               exit
            end if
         end do
      end do
    end function count_nnz_zB_A

    integer function count_nnz_zB_A_h2f(mat_size, rowptrA, colindA, rowptrB, colindB)
      implicit none
      integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
      integer :: i, j, jA, jB, posA, posB, posAend, posBend
    
      count_nnz_zB_A_h2f = 0
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
               posA = posA + 1
               posB = posB + 1
               if ( jA /= i ) count_nnz_zB_A_h2f = count_nnz_zB_A_h2f + 1
            else if ( jB < jA ) then
               posB = posB + 1
               if ( jB /= i ) count_nnz_zB_A_h2f = count_nnz_zB_A_h2f + 1
            else if ( jA < jB ) then
               posA = posA + 1
               if ( jA /= i ) count_nnz_zB_A_h2f = count_nnz_zB_A_h2f + 1
            end if
            count_nnz_zB_A_h2f = count_nnz_zB_A_h2f + 1
            
            if ( posA > posAend .and. posB > posBend ) then
               exit
            end if
         end do
      end do
    end function count_nnz_zB_A_h2f


    integer function count_nnz_zI_A(mat_size, rowptrA, diagptrA)
      implicit none
      integer, intent(in) :: mat_size, rowptrA(*), diagptrA(*)
      integer :: i
          
      count_nnz_zI_A = rowptrA(mat_size+1) - 1 ! nnz of A
      do i = 1, mat_size
         if ( diagptrA(i) == rowptrA(i+1) ) then ! means A(i,i) == 0
            count_nnz_zI_A = count_nnz_zI_A + 1
         end if
      end do
    end function count_nnz_zI_A


    integer function count_nnz_zI_A_h2f(mat_size, rowptrA, diagptrA)
      implicit none
      integer, intent(in) :: mat_size, rowptrA(*), diagptrA(*)
      integer :: i
    
      count_nnz_zI_A_h2f = 2*(rowptrA(mat_size+1) - 1)
      do i = 1, mat_size
         if ( diagptrA(i) == rowptrA(i+1) ) then ! means A(i,i) == 0
            count_nnz_zI_A_h2f = count_nnz_zI_A_h2f + 1
         else
            count_nnz_zI_A_h2f = count_nnz_zI_A_h2f - 1
         end if
      end do
    end function count_nnz_zI_A_h2f

  subroutine idx_COO_zB_A &
       (mat_size, rowptrA, colindA, rowptrB, colindB, rowout, colout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    integer, intent(out) :: rowout(*), colout(*)

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
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jA
             posA = posA + 1
             posB = posB + 1
          else if ( jB < jA ) then
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jB
             posB = posB + 1
          else if ( jA < jB ) then
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jA
             posA = posA + 1
          end if
          pos_zB_A = pos_zB_A + 1
                    
          if ( posA > posAend .and. posB > posBend ) then
             exit
          end if          
       end do
    end do    
  end subroutine idx_COO_zB_A

  subroutine idx_COO_zB_A_h2f &
       (mat_size, rowptrA, colindA, rowptrB, colindB, rowout, colout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), rowptrB(*), colindB(*)
    integer, intent(out) :: rowout(*), colout(*)

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
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jA
             if ( jA /= i ) then
                pos_zB_A = pos_zB_A + 1
                rowout(pos_zB_A) = jA
                colout(pos_zB_A) = i
             end if
             posA = posA + 1
             posB = posB + 1
          else if ( jB < jA ) then
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jB
             if ( jB /= i ) then
                pos_zB_A = pos_zB_A + 1
                rowout(pos_zB_A) = jB
                colout(pos_zB_A) = i
             end if
             posB = posB + 1
          else if ( jA < jB ) then
             rowout(pos_zB_A) = i
             colout(pos_zB_A) = jA
             if ( jA /= i ) then
                pos_zB_A = pos_zB_A + 1
                rowout(pos_zB_A) = jA
                colout(pos_zB_A) = i
             end if
             posA = posA + 1
          end if
          pos_zB_A = pos_zB_A + 1

          if ( posA > posAend .and. posB > posBend ) then
             exit
          end if          
       end do
    end do
  end subroutine idx_COO_zB_A_h2f

  subroutine idx_COO_zI_A &
       (mat_size, rowptrA, colindA, diagptrA, rowout, colout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), diagptrA(*)
    integer, intent(out) :: rowout(*), colout(*)
    
    integer :: i, j, pos
    
    pos = 1
    do i = 1, mat_size
       do j = rowptrA(i), diagptrA(i) - 1
          rowout(pos) = i
          colout(pos) = colindA(j)
          pos = pos + 1
       end do
       do j = diagptrA(i) + 1, rowptrA(i+1) - 1
          rowout(pos) = i
          colout(pos) = colindA(j)
          pos = pos + 1
       end do
       if ( diagptrA(i) /= rowptrA(i+1) ) then ! there is diagonal entry
          rowout(pos) = i
          colout(pos) = i
       else
          rowout(pos) = i
          colout(pos) = i
       end if
       pos = pos + 1          
    end do    
  end subroutine idx_COO_zI_A

  subroutine idx_COO_zI_A_h2f(mat_size, rowptrA, colindA, diagptrA, rowout, colout)
    implicit none
    integer, intent(in) :: mat_size, rowptrA(*), colindA(*), diagptrA(*)
    integer, intent(out) :: rowout(*), colout(*)
    
    integer :: i, j, pos
    
    pos = 1
    do i = 1, mat_size
       do j = rowptrA(i), diagptrA(i) - 1
          rowout(pos) = i
          colout(pos) = colindA(j)
          pos = pos + 1
          rowout(pos) = colindA(j)
          colout(pos) = i
          pos = pos + 1
       end do
       do j = diagptrA(i) + 1, rowptrA(i+1) - 1
          rowout(pos) = i
          colout(pos) = colindA(j)
          pos = pos + 1
          rowout(pos) = colindA(j)
          colout(pos) = i
          pos = pos + 1
       end do
       if ( diagptrA(i) /= rowptrA(i+1) ) then ! there is diagonal entry
          rowout(pos) = i
          colout(pos) = i
       else
          rowout(pos) = i
          colout(pos) = i
       end if
       pos = pos + 1          
    end do
  end subroutine idx_COO_zI_A_h2f

    subroutine get_diagptr(mat_size, rowptr, colind, diagptr)
      implicit none
      integer, intent(in) :: mat_size, rowptr(*), colind(*)
      integer, intent(out) :: diagptr(*)
      
      integer :: i, left, right, m

      do i = 1, mat_size
         diagptr(i) = rowptr(i+1) ! when there is no diagonal entry
         ! binary search
         left = rowptr(i)
         right = rowptr(i+1) - 1
         do while ( left <= right )
            m = (left + right) / 2
            if ( colind(m) == i ) then
               diagptr(i) = m
               exit
            else if ( colind(m) < i ) then
               left = m + 1
            else
               right = m - 1
            end if
         end do
      end do
    end subroutine get_diagptr

    
! include double complex routines
#include "zpares_mumps_inc.f90"
#define REALMAT
! include double precision routines
#include "zpares_mumps_inc.f90"
#define SINGLE
! include real routines
#include "zpares_mumps_inc.f90"
#undef REALMAT
! include complex routines
#include "zpares_mumps_inc.f90"

    subroutine set_mumps_icntl(icntl)
      implicit none
      integer :: icntl(:)
      !-- Options for MUMPS.
      ICNTL( 1) = -1 ! Output stream for error messages.
      ICNTL( 2) = -1 ! Output stream for diagonostic printing, statistics, and warning messages.
      ICNTL( 3) = -1 ! Output stream for global information, collected on the host.
      ICNTL( 5) =  0 ! Input matrix must be given in assembled format.
      ICNTL( 6) =  7 ! Automatic choice
      ICNTL( 7) =  7 ! Automatic ( for sequential analysis. 3: SCOTCH, 4: PORD, 5: MeTiS, 7: Automatic choice )
      ICNTL( 8) = 77 ! Automatic choice of scaling option may be performed.
      ICNTL( 9) =  1 ! Ax = b is solved.
      ICNTL(10) =  0 ! Number of iterative refinement ( only single right-hand side ).
      ICNTL(11) =  0
      ICNTL(12) =  0
      ICNTL(13) =  0
      ICNTL(14) = 50
      ICNTL(18) =  0 ! The inputted matrix is centralized on the host.
      ICNTL(19) =  0
      ICNTL(20) =  0 ! Centralized dense right-hand side
      ICNTL(21) =  0 ! Centralized dense solution
      ICNTL(22) =  0 ! Option for Out-Of-Core ( 0: Disable, 1: Enable )
      ICNTL(23) =  0
      ICNTL(24) =  0
      ICNTL(25) =  0
      ICNTL(26) =  0
      ICNTL(27) = -8
      ICNTL(28) =  0 ! Parallel analysis   ( 0: Automatic, 1: Sequential, 2: Parallel )
      ICNTL(29) =  0 ! Reordering strategy ( 0: Automatic, 1: PT-SCOTCH , 2: ParMeTiS )
      ICNTL(30) =  0
      ICNTL(31) =  0
      ICNTL(33) =  0    
  end subroutine set_mumps_icntl

  subroutine start_timer(t0)
    implicit none
    double precision, intent(out) :: t0
#ifndef MPI
    integer :: itc
    call system_clock(itc)
    t0 = dble(itc)
#else
    include 'mpif.h'
    t0 = mpi_wtime()
#endif
  end subroutine start_timer

  subroutine stop_timer(t0, t)
    implicit none
    double precision, intent(in) :: t0
    double precision, intent(inout) :: t
#ifdef MPI  
    include 'mpif.h'
    t = t + mpi_wtime() - t0
#else
    integer :: itc, t_rate, t_max
    double precision :: dtc, diff
    call system_clock(itc, t_rate, t_max)
    dtc = dble(itc)
    if ( dtc < t0 ) then
       diff = dble(t_max) - t0 + dtc
    else
       diff = dtc - t0
    end if
    diff = diff / dble(t_rate)
    t = t + diff
#endif
  end subroutine stop_timer

end module zpares_mumps

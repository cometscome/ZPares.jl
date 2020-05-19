#undef MATRIX_TYPE
#undef ZERO

#ifdef REALMAT
#ifdef SINGLE
#define MATRIX_TYPE real
#define ZERO 0.0
#else
#define MATRIX_TYPE double precision
#define ZERO 0d0
#endif
#else
#ifdef SINGLE
#define MATRIX_TYPE complex
#define ZERO (0.0,0.0)
#else
#define MATRIX_TYPE complex(kind(0d0))
#define ZERO (0d0,0d0)
#endif
#endif


#ifdef REALMAT
#ifdef SINGLE
  subroutine s_zero_UPLO &
#else
  subroutine d_zero_UPLO &
#endif
#else
#ifdef SINGLE
  subroutine c_zero_UPLO &
#else
  subroutine z_zero_UPLO &
#endif
#endif
  (mat_size, UPLO, Mat)
    implicit none
    integer, intent(in) :: mat_size
    character, intent(in) :: UPLO
    MATRIX_TYPE, intent(out) :: Mat(mat_size,*)
    integer :: i, j

    if ( UPLO == 'U' .or. UPLO == 'u' ) then
       do i = 1, mat_size
          do j = 1, i-1
             Mat(i,j) = ZERO
          end do
       end do
    elseif ( UPLO == 'L' .or. UPLO == 'l' ) then
       do i = 1, mat_size
          do j = 1, i-1
             Mat(j,i) = ZERO
          end do
       end do
    else
       ! error
    end if
  end subroutine


#ifdef REALMAT
#ifdef SINGLE
  subroutine s_coo2csr &
#else
  subroutine d_coo2csr &
#endif
#else
#ifdef SINGLE
  subroutine c_coo2csr &
#else
  subroutine z_coo2csr &
#endif
#endif
  (mat_size, nnz, indx, jndx, rcoo, ccoo, rowptr, colind, val)
    implicit none
    integer, intent(in) :: mat_size, nnz, jndx(*)
    integer, intent(inout) ::  indx(*)
    integer, intent(out) :: rowptr(*), colind(*)
    double precision, intent(in) :: rcoo(*)
    complex(kind(0d0)), intent(in) :: ccoo(*)
    MATRIX_TYPE, intent(out) :: val(*)

    integer :: k, kold, j, i, m
    MATRIX_TYPE :: t, tnext

    rowptr(1:mat_size+1) = 0
    do k = 1, nnz
       rowptr(indx(k)) = rowptr(indx(k)) + 1
    end do

    k = 1
    do j = 1, mat_size+1
       kold = rowptr(j)
       rowptr(j) = k
       k = k + kold
    end do

    do k = 1, nnz
       i = indx(k)
       j = jndx(k)
       m = rowptr(i)
       if ( field == 'complex' ) then
          val(m) = ccoo(k)
       else
          val(m) = rcoo(k)
       end if
       colind(m) = j
       rowptr(i) = m + 1
    end do

    do j = mat_size, 1, -1
       rowptr(j+1) = rowptr(j)
    end do
    rowptr(1) = 1

#ifdef REALMAT
#ifdef SINGLE
  call s_sort_csr &
#else
  call d_sort_csr &
#endif
#else
#ifdef SINGLE
  call c_sort_csr &
#else
  call z_sort_csr &
#endif
#endif
  (mat_size, rowptr, colind, val)
end subroutine


#ifdef REALMAT
#ifdef SINGLE
  subroutine s_sort_csr &
#else
  subroutine d_sort_csr &
#endif
#else
#ifdef SINGLE
  subroutine c_sort_csr &
#else
  subroutine z_sort_csr &
#endif
#endif
  (mat_size, rowptr, colind, val)  
  implicit none
  integer, intent(in) :: mat_size, rowptr(*)
  integer, intent(inout) :: colind(*)
  MATRIX_TYPE, intent(inout) :: val(*)
  
  integer :: i, st, en, j
  
  do i = 1, mat_size     
     st = rowptr(i)
     en = rowptr(i+1) - 1
     if ( st == en - 1 ) then
        cycle
     end if
#ifdef REALMAT
#ifdef SINGLE
     call s_heap_sort &
#else
     call d_heap_sort &
#endif
#else
#ifdef SINGLE
     call c_heap_sort &
#else
     call z_heap_sort &
#endif
#endif
     (en-st+1, colind(st:en), val(st:en))
  end do
end subroutine


#ifdef REALMAT
#ifdef SINGLE
subroutine s_heap_sort &
#else
subroutine d_heap_sort &
#endif
#else
#ifdef SINGLE
subroutine c_heap_sort &
#else
subroutine z_heap_sort &
#endif
#endif
  (length, colind_slice, val_slice)
  implicit none
  integer, intent(in) :: length
  integer, intent(inout) :: colind_slice(:)
  MATRIX_TYPE, intent(inout) :: val_slice(:)
  
  integer :: i, j, ir, l, n, itmp
  MATRIX_TYPE :: xtmp
  
  if ( length <= 1 ) then
     return
  end if

  l = ishft(length, -1) + 1
  ir = length
  
  do
     if ( l > 1 ) then
        l = l - 1
        itmp = colind_slice(l)
        xtmp = val_slice(l)
     else
        itmp = colind_slice(ir)
        xtmp = val_slice(ir)        
        colind_slice(ir) = colind_slice(1)
        val_slice(ir) = val_slice(1)
        ir = ir - 1
        if ( ir == 1 ) then
           colind_slice(1) = itmp
           val_slice(1) = xtmp
           return
        end if
     end if
     i = l
     j = ishft(l, 1)
     do while ( j <= ir )
        if ( j < ir ) then
           if ( colind_slice(j) < colind_slice(j+1) ) j = j + 1
        end if
        if ( itmp < colind_slice(j) ) then
           colind_slice(i) = colind_slice(j)
           val_slice(i) = val_slice(j)
           i = j
           j = j + i
        else
           j = ir + 1
        end if
     end do
     colind_slice(i) = itmp
     val_slice(i) = xtmp
  end do
end subroutine
  

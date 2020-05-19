module zpares_aux
  implicit none

  interface quad_ell_trap
     module procedure s_quad_ell_trap
     module procedure d_quad_ell_trap
  end interface

  interface ALLREDUCE_SUM_1D
     module procedure s_ALLREDUCE_SUM_1D
     module procedure d_ALLREDUCE_SUM_1D
     module procedure c_ALLREDUCE_SUM_1D
     module procedure z_ALLREDUCE_SUM_1D
  end interface

  interface ALLREDUCE_SUM_2D
     module procedure s_ALLREDUCE_SUM_2D
     module procedure d_ALLREDUCE_SUM_2D
     module procedure c_ALLREDUCE_SUM_2D
     module procedure z_ALLREDUCE_SUM_2D
  end interface

  interface norm2_blk
     module procedure s_norm2_blk
     module procedure d_norm2_blk
     module procedure c_norm2_blk
     module procedure z_norm2_blk
  end interface

  interface GEMM_ALLREDUCE
     module procedure SGEMM_ALLREDUCE
     module procedure DGEMM_ALLREDUCE
     module procedure CGEMM_ALLREDUCE
     module procedure ZGEMM_ALLREDUCE
  end interface
  
  interface create_hutch_samples
     module procedure s_create_hutch_samples
     module procedure c_create_hutch_samples
     module procedure d_create_hutch_samples
     module procedure z_create_hutch_samples
  end interface

  interface create_rand_matrix
     module procedure s_create_rand_matrix
     module procedure d_create_rand_matrix
     module procedure c_create_rand_matrix
     module procedure z_create_rand_matrix
  end interface

  interface orth_SVD
     module procedure s_orth_SVD
     module procedure d_orth_SVD
     module procedure c_orth_SVD
     module procedure z_orth_SVD
  end interface

  interface DOT_ALLREDUCE
     module procedure SDOT_ALLREDUCE
     module procedure DDOT_ALLREDUCE
     module procedure CDOT_ALLREDUCE
     module procedure ZDOT_ALLREDUCE
  end interface

  interface block_Hankel
     module procedure s_block_Hankel
     module procedure d_block_Hankel
     module procedure c_block_Hankel
     module procedure z_block_Hankel
  end interface

  interface serial_SVD
     module procedure s_serial_SVD
     module procedure d_serial_SVD
     module procedure c_serial_SVD
     module procedure z_serial_SVD
  end interface

  interface LAPACK_QR
     module procedure s_LAPACK_QR
     module procedure d_LAPACK_QR
     module procedure c_LAPACK_QR
     module procedure z_LAPACK_QR
  end interface

  interface GEGV_reduced_eig
     module procedure SGEGV_reduced_eig
     module procedure DGEGV_reduced_eig
     module procedure CGEGV_reduced_eig
     module procedure ZGEGV_reduced_eig
  end interface

  interface GEEV_reduced_eig
     module procedure SGEEV_reduced_eig
     module procedure DGEEV_reduced_eig
     module procedure CGEEV_reduced_eig
     module procedure ZGEEV_reduced_eig
  end interface

  interface HEGV_reduced_eig
     module procedure SSYGV_reduced_eig
     module procedure DSYGV_reduced_eig
     module procedure CHEGV_reduced_eig
     module procedure ZHEGV_reduced_eig
  end interface

  interface HEEV_reduced_eig
     module procedure SSYEV_reduced_eig
     module procedure DSYEV_reduced_eig
     module procedure CHEEV_reduced_eig
     module procedure ZHEEV_reduced_eig
  end interface

  interface basis_rotation
     module procedure s_basis_rotation
     module procedure d_basis_rotation
     module procedure c_basis_rotation
     module procedure z_basis_rotation
  end interface

  interface inside_ellipse
     module procedure s_inside_ellipse
     module procedure d_inside_ellipse
  end interface

  interface packing
     module procedure s_packing
     module procedure d_packing
     module procedure c_packing
     module procedure z_packing
  end interface

  interface calc_center_radius
     module procedure s_calc_center_radius
     module procedure d_calc_center_radius
  end interface

contains

! include double complex routines
#include "zpares_aux_inc.f90" 
#define REALMAT
! include double precision routines
#include "zpares_aux_inc.f90" 
#define SINGLE
! include real routines
#include "zpares_aux_inc.f90" 
#undef REALMAT
! include complex routines
#include "zpares_aux_inc.f90" 

end module zpares_aux

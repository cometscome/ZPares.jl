#undef REAL_TYPE
#undef COMPLEX_TYPE
#undef ZERO_R
#undef ZERO_C
#undef ONE_R
#undef ONE_C
#undef MATRIX_TYPE
#undef ZERO_M
#undef ONE_M
#undef MPI_TYPE
#undef PREFIX
#undef R_PREFIX
#undef C_PREFIX

#ifdef SINGLE
#define REAL_TYPE real
#define COMPLEX_TYPE complex
#define ZERO_R 0.0
#define ZERO_C (0.0,0.0)
#define ONE_R 1.0
#define ONE_C (1.0,0.0)
#define R_PREFIX s
#define C_PREFIX c
#ifdef REALMAT
#define MATRIX_TYPE real
#define ZERO_M 0.0
#define ONE_M 1.0
#define MPI_TYPE MPI_REAL
#define PREFIX s
#else
#define MATRIX_TYPE complex
#define ZERO_M (0.0,0.0)
#define ONE_M (1.0,0.0)
#define MPI_TYPE MPI_COMPLEX
#define PREFIX c
#endif
#else
#define REAL_TYPE double precision
#define COMPLEX_TYPE complex(kind(0d0))
#define ZERO_R 0d0
#define ZERO_C (0d0,0d0)
#define ONE_R 1d0
#define ONE_C (1d0,0d0)
#define R_PREFIX d
#define C_PREFIX z
#ifdef REALMAT
#define MATRIX_TYPE double precision
#define ZERO_M 0d0
#define ONE_M 1d0
#define MPI_TYPE MPI_DOUBLE_PRECISION
#define PREFIX d
#else
#define MATRIX_TYPE complex(kind(0d0))
#define ZERO_M (0d0,0d0)
#define ONE_M (1d0,0d0)
#define MPI_TYPE MPI_DOUBLE_COMPLEX
#define PREFIX z
#endif
#endif

#define MACRO_PASTE(A) A
#define MACRO_ADD_PRFX(STR) MACRO_PASTE(PREFIX)STR
#define MACRO_ADD_R_PRFX(STR) MACRO_PASTE(R_PREFIX)STR
#define MACRO_ADD_C_PRFX(STR) MACRO_PASTE(C_PREFIX)STR
#define MACRO_INSERT_PRFX(FIRST,LAST) MACRO_PASTE(FIRST)MACRO_PASTE(PREFIX)LAST

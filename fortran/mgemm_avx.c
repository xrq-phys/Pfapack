#include <mmintrin.h>
#include <xmmintrin.h>  // SSE
#include <pmmintrin.h>  // SSE2
#include <emmintrin.h>  // SSE3
#include <immintrin.h>  // AVX

/* Create macros so that the matrices are stored in column-major order */

#define A(i,j) a[ (j)*lda + (i) ]
#define B(i,j) b[ (j)*ldb + (i) ]
#define C(i,j) c[ (j)*ldc + (i) ]

typedef union
{
  __m256d blk;
  __m128d vec[2];
  double d[4];
} v4df_t;

#ifndef __COL_BASE_4x4

/*
 * Block-based avx 4x4 starts here.
 */
void adddott4x4_( int *k_, double *alpha_, 
        double *a, int *lda_,  double *b, int *ldb_, double *c, int *ldc_ )
{
  int p,
      k = *k_,
      lda = *lda_,
      ldb = *ldb_,
      ldc = *ldc_;

  v4df_t
    alpha_vreg,
    c_blk_00_vreg,    c_blk_01_vreg,
    c_blk_10_vreg,    c_blk_11_vreg,
    a_blk_0p_vreg,    a_blk_1p_vreg,
    b_blk_0p_vreg,    b_blk_1p_vreg;

  /* expand constant alpha */
  alpha_vreg.vec[0] = _mm_loaddup_pd( (double *) alpha_ );
  alpha_vreg.vec[1] = alpha_vreg.vec[0];

  c_blk_00_vreg.blk = _mm256_setzero_pd();   
  c_blk_10_vreg.blk = _mm256_setzero_pd();
  c_blk_01_vreg.blk = _mm256_setzero_pd(); 
  c_blk_11_vreg.blk = _mm256_setzero_pd(); 

  for ( p=0; p<k; p++ ){
    a_blk_0p_vreg.vec[0] = _mm_load_pd( (double*) a ); /* 1:2, p */
    a_blk_1p_vreg.vec[0] = _mm_load_pd( (double*) (a+2) ); /* 3:4, p */
    a_blk_0p_vreg.vec[1] = a_blk_0p_vreg.vec[0];
    a_blk_1p_vreg.vec[1] = a_blk_1p_vreg.vec[0];
    a += lda;

    b_blk_0p_vreg.vec[0] = _mm_loaddup_pd( (double*) b ); /* 1, p */
    b_blk_0p_vreg.vec[1] = _mm_loaddup_pd( (double*) (b+1) ); /* 2, p */
    b_blk_1p_vreg.vec[0] = _mm_loaddup_pd( (double*) (b+2) ); /* 3, p */
    b_blk_1p_vreg.vec[1] = _mm_loaddup_pd( (double*) (b+3) ); /* 4, p */
    b += ldb;

    c_blk_00_vreg.blk += a_blk_0p_vreg.blk * b_blk_0p_vreg.blk;
    c_blk_10_vreg.blk += a_blk_1p_vreg.blk * b_blk_0p_vreg.blk;
    c_blk_01_vreg.blk += a_blk_0p_vreg.blk * b_blk_1p_vreg.blk;
    c_blk_11_vreg.blk += a_blk_1p_vreg.blk * b_blk_1p_vreg.blk;
  }
  /* alpha factor */
  c_blk_00_vreg.blk *= alpha_vreg.blk;
  c_blk_10_vreg.blk *= alpha_vreg.blk;
  c_blk_01_vreg.blk *= alpha_vreg.blk;
  c_blk_11_vreg.blk *= alpha_vreg.blk;

  /* write back to memory */
  C( 0, 0 ) += c_blk_00_vreg.d[0];  C( 0, 1 ) += c_blk_00_vreg.d[2];
  C( 1, 0 ) += c_blk_00_vreg.d[1];  C( 1, 1 ) += c_blk_00_vreg.d[3];

  C( 2, 0 ) += c_blk_10_vreg.d[0];  C( 2, 1 ) += c_blk_10_vreg.d[2];
  C( 3, 0 ) += c_blk_10_vreg.d[1];  C( 3, 1 ) += c_blk_10_vreg.d[3];

  C( 0, 2 ) += c_blk_01_vreg.d[0];  C( 0, 3 ) += c_blk_01_vreg.d[2];
  C( 1, 2 ) += c_blk_01_vreg.d[1];  C( 1, 3 ) += c_blk_01_vreg.d[3];

  C( 2, 2 ) += c_blk_11_vreg.d[0];  C( 2, 3 ) += c_blk_11_vreg.d[2];
  C( 3, 2 ) += c_blk_11_vreg.d[1];  C( 3, 3 ) += c_blk_11_vreg.d[3];
}

#else

/* 
 * This is the column-based avx implementation.
 * Dot without transposing is not used thus not included here.
 */
void adddott4x4_( int *k_, double *alpha_, 
        double *a, int *lda_,  double *b, int *ldb_, double *c, int *ldc_ )
{
  int p,
      k = *k_,
      lda = *lda_,
      ldb = *ldb_,
      ldc = *ldc_;

  v4df_t
    alpha_vreg,
    c_col_0_vreg,    c_col_1_vreg,
    c_col_2_vreg,    c_col_3_vreg,
    a_col_p_vreg,
    b_p0_vreg, b_p1_vreg, b_p2_vreg, b_p3_vreg; 

  /* expand constant alpha */
  alpha_vreg.vec[0] = _mm_loaddup_pd( (double *) alpha_ );
  alpha_vreg.vec[1] = alpha_vreg.vec[0];

  c_col_0_vreg.blk = _mm256_setzero_pd();   
  c_col_1_vreg.blk = _mm256_setzero_pd();
  c_col_2_vreg.blk = _mm256_setzero_pd(); 
  c_col_3_vreg.blk = _mm256_setzero_pd(); 

  for ( p=0; p<k; p++ ){
    /* load row */
    a_col_p_vreg.blk = _mm256_load_pd( (double *) a );
    a += lda;

    /* load and duplicate columns */
    b_p0_vreg.vec[0] = _mm_loaddup_pd( (double *) b );
    b_p1_vreg.vec[0] = _mm_loaddup_pd( (double *) ( b+1 ) );
    b_p2_vreg.vec[0] = _mm_loaddup_pd( (double *) ( b+2 ) );
    b_p3_vreg.vec[0] = _mm_loaddup_pd( (double *) ( b+3 ) );
    b_p0_vreg.vec[1] = b_p0_vreg.vec[0];
    b_p1_vreg.vec[1] = b_p1_vreg.vec[0];
    b_p2_vreg.vec[1] = b_p2_vreg.vec[0];
    b_p3_vreg.vec[1] = b_p3_vreg.vec[0];
    b += ldb;

    c_col_0_vreg.blk += a_col_p_vreg.blk * b_p0_vreg.blk;
    c_col_1_vreg.blk += a_col_p_vreg.blk * b_p1_vreg.blk;
    c_col_2_vreg.blk += a_col_p_vreg.blk * b_p2_vreg.blk;
    c_col_3_vreg.blk += a_col_p_vreg.blk * b_p3_vreg.blk;
  }
  /* alpha factor */
  c_col_0_vreg.blk *= alpha_vreg.blk;
  c_col_1_vreg.blk *= alpha_vreg.blk;
  c_col_2_vreg.blk *= alpha_vreg.blk;
  c_col_3_vreg.blk *= alpha_vreg.blk;

  /* write back to memory */
  C( 0, 0 ) += c_col_0_vreg.d[0];  C( 0, 1 ) += c_col_1_vreg.d[0];
  C( 1, 0 ) += c_col_0_vreg.d[1];  C( 1, 1 ) += c_col_1_vreg.d[1];

  C( 2, 0 ) += c_col_0_vreg.d[2];  C( 2, 1 ) += c_col_1_vreg.d[2];
  C( 3, 0 ) += c_col_0_vreg.d[3];  C( 3, 1 ) += c_col_1_vreg.d[3];

  C( 0, 2 ) += c_col_2_vreg.d[0];  C( 0, 3 ) += c_col_3_vreg.d[0];
  C( 1, 2 ) += c_col_2_vreg.d[1];  C( 1, 3 ) += c_col_3_vreg.d[1];

  C( 2, 2 ) += c_col_2_vreg.d[2];  C( 2, 3 ) += c_col_3_vreg.d[2];
  C( 3, 2 ) += c_col_2_vreg.d[3];  C( 3, 3 ) += c_col_3_vreg.d[3];
}

#endif


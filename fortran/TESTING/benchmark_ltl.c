extern void dcopy_( const int *n, double *dx, const int *incx, double *dy, const int *incy );
extern void daxpy_( const int *n, const double *da, double *dx, const int *incx, double *dy, const int *incy );
extern void dlarnv_( const int *idist, int *iseed4, const int *n, double *x );
extern double dlange_( const char *norm, const int *m, const int *n, double *a, const int *lda, double *work );		
extern void dsktrf_( const char *uplo, const char *mode, const int *n, double *a, const int *lda, int *ipiv, double *work, const int *lwork, int *info );

int dsktrf( const char uplo, const char mode, const int n, double *a, const int lda, int *ipiv, double *work, const int lwork ) {
  int info = 0;
  dsktrf_( &uplo, &mode, &n, a, &lda, ipiv, work, &lwork, &info );

  return info;
}


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct timespec timespec_t;

int64_t timediff_ns(timespec_t start, timespec_t end)
{
  timespec_t temp;
  if ((end.tv_nsec - start.tv_nsec) < 0) {
    temp.tv_sec = end.tv_sec - start.tv_sec - 1;
    temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  }
  return (int64_t)temp.tv_sec * 1000000000 + temp.tv_nsec;
}

int main(const int argc, const char *argv[]) {
  int nmax, nmin, nstep;
  bool blocked = false;
  if ( argc < 4 ) {
    fprintf(stderr, "./.x nmax nmin nstep [blocked=False]\n");
    return -1;
  }
  nmax  = atoi(argv[1]);
  nmin  = atoi(argv[2]);
  nstep = atoi(argv[3]);
  if ( argc >= 5 ) {
    if ( argv[4][0] == 'T' || argv[4][0] == 't' ) {
      blocked = true;
    }
  }
  int nexp = 10;

  double *A = (double *)malloc(sizeof(double) * nmax * nmax * 3);
  double *A2 = A + nmax * nmax;
  double *A3 = A2 + nmax * nmax;
  double *W = (double *)malloc(sizeof(double) * nmax * nmax);
  int *iPiv2 = (int *)malloc(sizeof(int) * nmax);
  int *iPiv3 = (int *)malloc(sizeof(int) * nmax);

  {
    int idist = 3;
    int iseed[4] = { 1, 1234, 4567, 123 };
    int siz = nmax * nmax;
    dlarnv_( &idist, iseed, &siz, A );
  }
  fprintf(stdout, "#n, right- vs. left-looking diff., right-looking msec, left-looking msec, ipiv mismatch\n");

  for ( int n = nmax; n >= nmin; n -= nstep ) {
    int lwork_r, lwork_l;
    int siz = n * n;
    int inc = 1;
    int npm = 0;
    double dmone = -1.0;
    dcopy_( &siz, A, &inc, A2, &inc );
    dcopy_( &siz, A, &inc, A3, &inc );

    if ( !blocked ) {
      lwork_r = 1;
      lwork_l = n;
    } else {
      lwork_r = nmax * nmax;
      lwork_l = nmax * nmax;
    }

    dsktrf( 'l', 'n', n, A2, n, iPiv2, W, lwork_r );
    dsktrf( 'l', 'l', n, A3, n, iPiv3, W, lwork_l );

    daxpy_( &siz, &dmone, A2, &inc, A3, &inc );
    for ( int i = 0; i < n; ++i ) { A3[i + i * n] = 0; }

    double diff_norm = dlange_( "F", &n, &n, A3, &n, W );
    for ( int i = 0; i < n; ++i ) { npm += ( iPiv2[i] != iPiv3[i] ); }

    timespec_t time1, time2, time3;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for ( int i = 0; i < nexp; ++i ) {
      dsktrf( 'l', 'l', n, A3, n, iPiv3, W, lwork_l );
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    for ( int i = 0; i < nexp; ++i ) {
      dsktrf( 'l', 'n', n, A2, n, iPiv2, W, lwork_r );
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);

    fprintf(stdout, "%8d %18.14e %18.4e %18.4e %d\n", n, diff_norm,
        (double)timediff_ns(time2, time3) / nexp / 1000000,
        (double)timediff_ns(time1, time2) / nexp / 1000000,
        npm);
  }
  // Zero placeholder
  fprintf(stdout, "%8d %18.14e %18.4e %18.4e %d\n", 0, 0.0, 1e-20, 1e-20, 0);

  free(A);
  free(W);
  return 0;
}


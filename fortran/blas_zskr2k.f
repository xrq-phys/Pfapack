      SUBROUTINE ZSKR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      DOUBLE COMPLEX MALPHA
      DOUBLE COMPLEX BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  ZSKR2K  performs one of the skew-symmetric rank 2k operations
*
*     C := alpha*A*B^T - alpha*B*A^T + beta*C,
*
*  or
*
*     C := alpha*A^T*B - alpha*B^T*A + beta*C,
*
*  where  alpha and beta are scalars,  C is an  n by n
*  skew-symmetric matrix and  A and B  are  n by k matrices in the first case
*  and  k by n  matrices in the second case.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'    C := alpha*A*B^T   -
*                                         alpha*B*A^T   +
*                                         beta*C.
*
*              TRANS = 'T' or 't'    C := alpha*A^T*B   -
*                                         alpha*B^T*A   +
*                                         beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
*           matrices  A and B.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE COMPLEX  array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE COMPLEX   array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE COMPLEX          .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE COMPLEX  array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  skew-symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  skew-symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the diagonal elements need
*           not be set,  they are assumed to be zero,  and on exit they
*           are set to zero.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
      EXTERNAL ZGEMMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP1,TEMP2
      INTEGER I,INFO,J,L,NROWA
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      IF (LSAME(TRANS,'N')) THEN
          NROWA = N
      ELSE
          NROWA = K
      END IF
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
          INFO = 1
      ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND.
     +         (.NOT.LSAME(TRANS,'T'))) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (K.LT.0) THEN
          INFO = 4
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
          INFO = 7
      ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
          INFO = 9
      ELSE IF (LDC.LT.MAX(1,N)) THEN
          INFO = 12
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZSKR2K',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR.
     +    (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
          IF (UPPER) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 20 J = 1,N
                      DO 10 I = 1,J
                          C(I,J) = ZERO
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 J = 1,N
                      DO 30 I = 1,J - 1
                          C(I,J) = BETA*C(I,J)
   30                 CONTINUE
                      C(J,J) = ZERO
   40             CONTINUE
              END IF
          ELSE
              IF (BETA.EQ.ZERO) THEN
                  DO 60 J = 1,N
                      DO 50 I = J,N
                          C(I,J) = ZERO
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 J = 1,N
                      C(J,J) = ZERO
                      DO 70 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
*
*     Negate ALPHA.
*
      MALPHA = -ALPHA
*
*     Call GEMMT.
*
      IF (LSAME(TRANS,'N')) THEN
         CALL ZGEMMT(UPLO,'N','T',N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         CALL ZGEMMT(UPLO,'N','T',N,K,MALPHA,B,LDB,A,LDA,ONE,C,LDC)
      ELSE
         CALL ZGEMMT(UPLO,'T','N',N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
         CALL ZGEMMT(UPLO,'T','N',N,K,MALPHA,B,LDB,A,LDA,ONE,C,LDC)
      END IF
*
      RETURN
*
*     End of ZSKR2K.
*
      END

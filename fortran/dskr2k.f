      SUBROUTINE DSKR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION BETA
      INTEGER K,LDA,LDB,LDC,N
      CHARACTER TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  DSKR2K  performs one of the skew-symmetric rank 2k operations
*
*     C := alpha*A*B^T - alpha*B*A^T + beta*C,
*
*  or
*
*     C := alpha*A^T*B - alpha*B^T*A + beta*C,
*
*  where  alpha and beta  are scalars,  C is an  n by n
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
*  ALPHA  - DOUBLE PRECISION         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION       array of DIMENSION ( LDA, ka ), where ka is
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
*  B      - DOUBLE PRECISION       array of DIMENSION ( LDB, kb ), where kb is
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
*  BETA   - DOUBLE PRECISION          .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION          array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  skew-symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           --- After 2020 patch:
*             Lower triangular part serves as scratchpad.
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
*  Level 3 Blas routine.
*
*  -- Written on 10/22/2010
*     Michael Wimmer, Universiteit Leiden
*     Based on ZHER2K from BLAS (www.netlib.org)
*  -- Patched on 03/07/2020
*     RuQing Xu, The University of Tokyo
*     Do blocking for better performance (ref. github.com/michael-lehn/ulmBLAS).
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP1,TEMP2,TBETA
      INTEGER I,INFO,J,L,NROWA,MI,MJ,ML,NBLK,NBLK_,KBLK,KBLK_,
     +        LENI,LENJ,LENL,I_,J_,L_
      LOGICAL UPPER
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE= 1.0D+0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO= 0.0D+0)
      INTEGER MBLK
*     ..
*
*     Prepare scratchpad space.
*
      IF (LDC.GE.66) THEN
*         SPMA =>C(3,1)
*         SPMB =>C(3,2)
          MBLK = 8
      ELSE
*         Do not do blocking.
          MBLK = N
      END IF
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
          CALL XERBLA('DSKR2K',INFO)
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
*     Blocking scheme.
*
      NBLK = N/MBLK
      NBLK_ = MOD(N,MBLK)
      IF (NBLK_.GT.0) THEN
          NBLK = NBLK+1
      END IF
      KBLK = K/MBLK
      KBLK_ = MOD(K,MBLK)
      IF (KBLK_.GT.0) THEN
          KBLK = KBLK+1
      END IF
*
*     Start the operations.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  C := alpha*A* B^T - alpha*B*A^T + beta*
*                   C.
*
          IF (UPPER) THEN
*
*           Blocked algorithm starts here.
*
            DO 131 MJ = 1,NBLK
              IF (MJ.EQ.NBLK.AND.NBLK_.GT.0) THEN
                LENJ = NBLK_
              ELSE
                LENJ = MBLK
              END IF
              DO 121 ML = 1,KBLK
                IF (ML.EQ.KBLK.AND.KBLK_.GT.0) THEN
                  LENL = KBLK_
                ELSE
                  LENL = MBLK
                END IF
*
*               Only multiplies BETA at the first block.
*
                IF (1.EQ.ML) THEN
                  TBETA = BETA
                ELSE
                  TBETA = 1
                END IF
*
*               Packs memory for fast execution.
*
                IF (NBLK.NE.1) 
     +            CALL DMPACK('T',LENL,LENJ,ALPHA,
     +                        B(MJ*MBLK-MBLK+1,ML*MBLK-MBLK+1),LDB,
     +                        C(3,2))
*
                DO 111 MI = 1,MJ
                  IF (MI.EQ.NBLK.AND.NBLK_.GT.0) THEN
                    LENI = NBLK_
                  ELSE
                    LENI = MBLK
                  END IF
                  IF (NBLK.NE.1)
     +              CALL DMPACK('N',LENI,LENL,ALPHA,
     +                          A(MI*MBLK-MBLK+1,ML*MBLK-MBLK+1),LDA,
     +                          C(3,1))
*
                  IF (MI.EQ.MJ) THEN
*
*                   Vanilla kernel along the diagonal.
*
                    CALL DMSKR2K(LENI,LENJ,LENL,
     +                           C(3,1),C(3,2),TBETA,
     +                           C(MI*MBLK-MBLK+1,MJ*MBLK-MBLK+1),LDC)
                  ELSE
*
*                   Use GEMM core here. 
*
                    CALL DMGEMM(LENI,LENJ,LENL,.TRUE.,
     +                          C(3,1),C(3,2),TBETA,
     +                          C(MI*MBLK-MBLK+1,MJ*MBLK-MBLK+1),LDC)
                  END IF
  111           CONTINUE
  121         CONTINUE
  131       CONTINUE
*
            IF (NBLK.GE.2) THEN
*
              TBETA = 1
*
              DO 132 MJ = 2,NBLK
                IF (MJ.EQ.NBLK.AND.NBLK_.GT.0) THEN
                  LENJ = NBLK_
                ELSE
                  LENJ = MBLK
                END IF
                DO 122 ML = 1,KBLK
                  IF (ML.EQ.KBLK.AND.KBLK_.GT.0) THEN
                    LENL = KBLK_
                  ELSE
                    LENL = MBLK
                  END IF
                  IF (NBLK.NE.1) 
     +              CALL DMPACK('T',LENL,LENJ,ALPHA,
     +                          A(MJ*MBLK-MBLK+1,ML*MBLK-MBLK+1),LDA,
     +                          C(3,2))
*
                  DO 112 MI = 1,MJ-1
                    IF (MI.EQ.NBLK.AND.NBLK_.GT.0) THEN
                      LENI = NBLK_
                    ELSE
                      LENI = MBLK
                    END IF
                    IF (NBLK.NE.1)
     +                CALL DMPACK('N',LENI,LENL,ALPHA,
     +                            B(MI*MBLK-MBLK+1,ML*MBLK-MBLK+1),LDB,
     +                            C(3,1))
*
                      CALL DMGEMM(LENI,LENJ,LENL,.FALSE.,
     +                            C(3,1),C(3,2),TBETA,
     +                            C(MI*MBLK-MBLK+1,MJ*MBLK-MBLK+1),LDC)
  112             CONTINUE
  122           CONTINUE
  132         CONTINUE
            END IF
          ELSE
              DO 180 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 140 I = J,N
                          C(I,J) = ZERO
  140                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 150 I = J + 1,N
                          C(I,J) = BETA*C(I,J)
  150                 CONTINUE
                      C(J,J) = ZERO
                  ELSE
                      C(J,J) = ZERO
                  END IF
                  DO 170 L = 1,K
                      IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                          TEMP1 = ALPHA*B(J,L)
                          TEMP2 = ALPHA*A(J,L)
                          DO 160 I = J + 1,N
                              C(I,J) = C(I,J) + A(I,L)*TEMP1 -
     +                                 B(I,L)*TEMP2
  160                     CONTINUE
                          C(J,J) = ZERO
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
*
*        Form  C := alpha*A^T*B - alpha*B^T*A + beta*
*                   C.
*
          IF (UPPER) THEN
              DO 210 J = 1,N
                  DO 200 I = 1,J
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 190 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  190                 CONTINUE
                      IF (I.EQ.J) THEN
                         C(J,J) = ZERO
                      ELSE
                          IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP1 - ALPHA*TEMP2
                          ELSE
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 -
     +                                 ALPHA*TEMP2
                          END IF
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 J = 1,N
                  DO 230 I = J,N
                      TEMP1 = ZERO
                      TEMP2 = ZERO
                      DO 220 L = 1,K
                          TEMP1 = TEMP1 + A(L,I)*B(L,J)
                          TEMP2 = TEMP2 + B(L,I)*A(L,J)
  220                 CONTINUE
                      IF (I.EQ.J) THEN
                         C(J,J) = ZERO
                      ELSE
                          IF (BETA.EQ.ZERO) THEN
                              C(I,J) = ALPHA*TEMP1 - ALPHA*TEMP2
                          ELSE
                              C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 -
     +                                 ALPHA*TEMP2
                          END IF
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DSKR2K.
*
      END
*
*     Auxilliary core MGEMM
*
      SUBROUTINE DMGEMM(LENI,LENJ,LENL,PM,A,B,BETA,C,LDC)
*     .. Arguments ..
      INTEGER LENI,LENJ,LENL,LDC
      LOGICAL PM
      DOUBLE PRECISION BETA
      DOUBLE PRECISION C(LDC,*)
      DOUBLE PRECISION A(LENI,*),B(LENJ,*)
*     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE= 1.0D+0)
*     .. Local variables ..
      DOUBLE PRECISION TEMP,TEMP2
      INTEGER I,J,L
*
*     Execution of microkernel.
*
      IF (LENI.EQ.8.AND.LENJ.EQ.8) THEN
*         .. Beta Update ..
          C(1:8,1:8) = C(1:8,1:8)*BETA
*
*         Inner part needs manual expanding.
*         Julia is usually better at this.
*
          IF (PM) THEN
              DO 302 L = 1,LENL
*                 .. Column 1 ..
                  C(1,1) = C(1,1) + A(1,L) * B(L,1);
                  C(2,1) = C(2,1) + A(2,L) * B(L,1);
                  C(3,1) = C(3,1) + A(3,L) * B(L,1);
                  C(4,1) = C(4,1) + A(4,L) * B(L,1);
                  C(5,1) = C(5,1) + A(5,L) * B(L,1);
                  C(6,1) = C(6,1) + A(6,L) * B(L,1);
                  C(7,1) = C(7,1) + A(7,L) * B(L,1);
                  C(8,1) = C(8,1) + A(8,L) * B(L,1);
*                 .. Column 2 ..
                  C(1,2) = C(1,2) + A(1,L) * B(L,2);
                  C(2,2) = C(2,2) + A(2,L) * B(L,2);
                  C(3,2) = C(3,2) + A(3,L) * B(L,2);
                  C(4,2) = C(4,2) + A(4,L) * B(L,2);
                  C(5,2) = C(5,2) + A(5,L) * B(L,2);
                  C(6,2) = C(6,2) + A(6,L) * B(L,2);
                  C(7,2) = C(7,2) + A(7,L) * B(L,2);
                  C(8,2) = C(8,2) + A(8,L) * B(L,2);
*                 .. Column 3 ..
                  C(1,3) = C(1,3) + A(1,L) * B(L,3);
                  C(2,3) = C(2,3) + A(2,L) * B(L,3);
                  C(3,3) = C(3,3) + A(3,L) * B(L,3);
                  C(4,3) = C(4,3) + A(4,L) * B(L,3);
                  C(5,3) = C(5,3) + A(5,L) * B(L,3);
                  C(6,3) = C(6,3) + A(6,L) * B(L,3);
                  C(7,3) = C(7,3) + A(7,L) * B(L,3);
                  C(8,3) = C(8,3) + A(8,L) * B(L,3);
*                 .. Column 4 ..
                  C(1,4) = C(1,4) + A(1,L) * B(L,4);
                  C(2,4) = C(2,4) + A(2,L) * B(L,4);
                  C(3,4) = C(3,4) + A(3,L) * B(L,4);
                  C(4,4) = C(4,4) + A(4,L) * B(L,4);
                  C(5,4) = C(5,4) + A(5,L) * B(L,4);
                  C(6,4) = C(6,4) + A(6,L) * B(L,4);
                  C(7,4) = C(7,4) + A(7,L) * B(L,4);
                  C(8,4) = C(8,4) + A(8,L) * B(L,4);
*                 .. Column 5 ..
                  C(1,5) = C(1,5) + A(1,L) * B(L,5);
                  C(2,5) = C(2,5) + A(2,L) * B(L,5);
                  C(3,5) = C(3,5) + A(3,L) * B(L,5);
                  C(4,5) = C(4,5) + A(4,L) * B(L,5);
                  C(5,5) = C(5,5) + A(5,L) * B(L,5);
                  C(6,5) = C(6,5) + A(6,L) * B(L,5);
                  C(7,5) = C(7,5) + A(7,L) * B(L,5);
                  C(8,5) = C(8,5) + A(8,L) * B(L,5);
*                 .. Column 6 ..
                  C(1,6) = C(1,6) + A(1,L) * B(L,6);
                  C(2,6) = C(2,6) + A(2,L) * B(L,6);
                  C(3,6) = C(3,6) + A(3,L) * B(L,6);
                  C(4,6) = C(4,6) + A(4,L) * B(L,6);
                  C(5,6) = C(5,6) + A(5,L) * B(L,6);
                  C(6,6) = C(6,6) + A(6,L) * B(L,6);
                  C(7,6) = C(7,6) + A(7,L) * B(L,6);
                  C(8,6) = C(8,6) + A(8,L) * B(L,6);
*                 .. Column 7 ..
                  C(1,7) = C(1,7) + A(1,L) * B(L,7);
                  C(2,7) = C(2,7) + A(2,L) * B(L,7);
                  C(3,7) = C(3,7) + A(3,L) * B(L,7);
                  C(4,7) = C(4,7) + A(4,L) * B(L,7);
                  C(5,7) = C(5,7) + A(5,L) * B(L,7);
                  C(6,7) = C(6,7) + A(6,L) * B(L,7);
                  C(7,7) = C(7,7) + A(7,L) * B(L,7);
                  C(8,7) = C(8,7) + A(8,L) * B(L,7);
*                 .. Column 8 ..
                  C(1,8) = C(1,8) + A(1,L) * B(L,8);
                  C(2,8) = C(2,8) + A(2,L) * B(L,8);
                  C(3,8) = C(3,8) + A(3,L) * B(L,8);
                  C(4,8) = C(4,8) + A(4,L) * B(L,8);
                  C(5,8) = C(5,8) + A(5,L) * B(L,8);
                  C(6,8) = C(6,8) + A(6,L) * B(L,8);
                  C(7,8) = C(7,8) + A(7,L) * B(L,8);
                  C(8,8) = C(8,8) + A(8,L) * B(L,8);
  302         CONTINUE
          ELSE
              DO 312 L = 1,LENL
*                 .. Column 1 ..
                  C(1,1) = C(1,1) - A(1,L) * B(L,1);
                  C(2,1) = C(2,1) - A(2,L) * B(L,1);
                  C(3,1) = C(3,1) - A(3,L) * B(L,1);
                  C(4,1) = C(4,1) - A(4,L) * B(L,1);
                  C(5,1) = C(5,1) - A(5,L) * B(L,1);
                  C(6,1) = C(6,1) - A(6,L) * B(L,1);
                  C(7,1) = C(7,1) - A(7,L) * B(L,1);
                  C(8,1) = C(8,1) - A(8,L) * B(L,1);
*                 .. Column 2 ..                   
                  C(1,2) = C(1,2) - A(1,L) * B(L,2);
                  C(2,2) = C(2,2) - A(2,L) * B(L,2);
                  C(3,2) = C(3,2) - A(3,L) * B(L,2);
                  C(4,2) = C(4,2) - A(4,L) * B(L,2);
                  C(5,2) = C(5,2) - A(5,L) * B(L,2);
                  C(6,2) = C(6,2) - A(6,L) * B(L,2);
                  C(7,2) = C(7,2) - A(7,L) * B(L,2);
                  C(8,2) = C(8,2) - A(8,L) * B(L,2);
*                 .. Column 3 ..                   
                  C(1,3) = C(1,3) - A(1,L) * B(L,3);
                  C(2,3) = C(2,3) - A(2,L) * B(L,3);
                  C(3,3) = C(3,3) - A(3,L) * B(L,3);
                  C(4,3) = C(4,3) - A(4,L) * B(L,3);
                  C(5,3) = C(5,3) - A(5,L) * B(L,3);
                  C(6,3) = C(6,3) - A(6,L) * B(L,3);
                  C(7,3) = C(7,3) - A(7,L) * B(L,3);
                  C(8,3) = C(8,3) - A(8,L) * B(L,3);
*                 .. Column 4 ..                   
                  C(1,4) = C(1,4) - A(1,L) * B(L,4);
                  C(2,4) = C(2,4) - A(2,L) * B(L,4);
                  C(3,4) = C(3,4) - A(3,L) * B(L,4);
                  C(4,4) = C(4,4) - A(4,L) * B(L,4);
                  C(5,4) = C(5,4) - A(5,L) * B(L,4);
                  C(6,4) = C(6,4) - A(6,L) * B(L,4);
                  C(7,4) = C(7,4) - A(7,L) * B(L,4);
                  C(8,4) = C(8,4) - A(8,L) * B(L,4);
*                 .. Column 5 ..                   
                  C(1,5) = C(1,5) - A(1,L) * B(L,5);
                  C(2,5) = C(2,5) - A(2,L) * B(L,5);
                  C(3,5) = C(3,5) - A(3,L) * B(L,5);
                  C(4,5) = C(4,5) - A(4,L) * B(L,5);
                  C(5,5) = C(5,5) - A(5,L) * B(L,5);
                  C(6,5) = C(6,5) - A(6,L) * B(L,5);
                  C(7,5) = C(7,5) - A(7,L) * B(L,5);
                  C(8,5) = C(8,5) - A(8,L) * B(L,5);
*                 .. Column 6 ..                   
                  C(1,6) = C(1,6) - A(1,L) * B(L,6);
                  C(2,6) = C(2,6) - A(2,L) * B(L,6);
                  C(3,6) = C(3,6) - A(3,L) * B(L,6);
                  C(4,6) = C(4,6) - A(4,L) * B(L,6);
                  C(5,6) = C(5,6) - A(5,L) * B(L,6);
                  C(6,6) = C(6,6) - A(6,L) * B(L,6);
                  C(7,6) = C(7,6) - A(7,L) * B(L,6);
                  C(8,6) = C(8,6) - A(8,L) * B(L,6);
*                 .. Column 7 ..                   
                  C(1,7) = C(1,7) - A(1,L) * B(L,7);
                  C(2,7) = C(2,7) - A(2,L) * B(L,7);
                  C(3,7) = C(3,7) - A(3,L) * B(L,7);
                  C(4,7) = C(4,7) - A(4,L) * B(L,7);
                  C(5,7) = C(5,7) - A(5,L) * B(L,7);
                  C(6,7) = C(6,7) - A(6,L) * B(L,7);
                  C(7,7) = C(7,7) - A(7,L) * B(L,7);
                  C(8,7) = C(8,7) - A(8,L) * B(L,7);
*                 .. Column 8 ..                   
                  C(1,8) = C(1,8) - A(1,L) * B(L,8);
                  C(2,8) = C(2,8) - A(2,L) * B(L,8);
                  C(3,8) = C(3,8) - A(3,L) * B(L,8);
                  C(4,8) = C(4,8) - A(4,L) * B(L,8);
                  C(5,8) = C(5,8) - A(5,L) * B(L,8);
                  C(6,8) = C(6,8) - A(6,L) * B(L,8);
                  C(7,8) = C(7,8) - A(7,L) * B(L,8);
                  C(8,8) = C(8,8) - A(8,L) * B(L,8);
  312         CONTINUE
          END IF
*
      ELSE
*
        DO 300 J = 1,LENJ
          IF (BETA.NE.ONE) THEN
              DO I = 1,LENI
                  C(I,J) = C(I,J)*BETA
              END DO
          END IF
          IF (PM) THEN
              DO 320 L = 1,LENL
                  TEMP = B(L,J)
                  DO 330 I = 1,LENI
                      C(I,J) = C(I,J) + TEMP*A(I,L)
  330             CONTINUE
  320         CONTINUE
          ELSE
              DO 321 L = 1,LENL
                  TEMP = B(L,J)
                  DO 331 I = 1,LENI
                      C(I,J) = C(I,J) - TEMP*A(I,L)
  331             CONTINUE
  321         CONTINUE
          END IF
  300   CONTINUE
      END IF
*
      RETURN
      END
*
*     Auxilliary core MSKR2K('N','U').
*
      SUBROUTINE DMSKR2K(LENI,LENJ,LENL,A,B,BETA,C,LDC)
*     .. Arguments ..
      INTEGER LENI,LENJ,LENL,LDC
      DOUBLE PRECISION BETA
      DOUBLE PRECISION C(LDC,*)
      DOUBLE PRECISION A(LENI,*),B(LENJ,*)
*     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE= 1.0D+0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO= 0.0D+0)
*     .. Local variables ..
      DOUBLE PRECISION TEMP1,TEMP2
      INTEGER I,J,L
*
*     Execution of microkernel.
*
      DO 130 J = 1,LENJ
          IF (BETA.EQ.ZERO) THEN
              DO 90 I = 1,J
                  C(I,J) = ZERO
   90         CONTINUE
          ELSE IF (BETA.NE.ONE) THEN
              DO 100 I = 1,J - 1
                  C(I,J) = BETA*C(I,J)
  100         CONTINUE
              C(J,J) = ZERO
          ELSE
              C(J,J) = ZERO
          END IF
          DO 120 L = 1,LENL
              IF ((A(J,L).NE.ZERO)  .OR. 
     +            (B(L,J).NE.ZERO)) THEN
                  TEMP1 = B(L,J)
                  TEMP2 = A(J,L)
                  DO 110 I = 1,J - 1
                      C(I,J) = C(I,J) + A(I,L)*TEMP1 - B(L,I)*TEMP2
  110             CONTINUE
                  C(J,J) = ZERO
              END IF
  120     CONTINUE
  130 CONTINUE
*
      RETURN
      END
*
*     Memory packing.
*
      SUBROUTINE DMPACK(TRANS,M,N,ALPHA,A,LDA,T)
*     .. Arguments ..
      CHARACTER TRANS
      INTEGER M,N,LDA,I
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION A(LDA,*),T(M,*)
      DOUBLE PRECISION ONE
      PARAMETER (ONE= 1.0D+0)
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*
*     Pack to tight memory.
*
      IF (LSAME(TRANS,'N')) THEN
          IF (ALPHA.EQ.ONE) THEN
              T(1:M,1:N) = A(1:M,1:N)
          ELSE
              T(1:M,1:N) = ALPHA*A(1:M,1:N)
          END IF
      ELSE
*
*         A transpose packing which could cause memory stall.
*
          IF (ALPHA.EQ.ONE) THEN
!GCC$ UNROLL 8
!DIR$ UNROLL =8
            DO 980 I = 1,M
              T(I,1:N) = A(1:N,I)
  980       CONTINUE
          ELSE
!GCC$ UNROLL 8
!DIR$ UNROLL =8
            DO 981 I = 1,M
              T(I,1:N) = ALPHA*A(1:N,I)
  981       CONTINUE
          END IF
      END IF
      END


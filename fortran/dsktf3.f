      SUBROUTINE DSKTF3( UPLO, MODE, N, NB, A, LDA, IPIV, W, LDW,
     $                   INIT, GLOBALK, INFO )
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, LDW, N, NB, GLOBALK
      LOGICAL            INIT
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION   W( LDW, * )
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     .. Local Scalars ..
      LOGICAL            UPPER, NORMAL
      INTEGER            K, KK, KP
      DOUBLE PRECISION   COLMAX, DDOT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX, DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSKR2, DGEMV, DGEMMT, XERBLA
      EXTERNAL           DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX

*     .. Placeholder for left-looking variant: Always FALSE here ..
      NORMAL = LSAME( MODE, 'N' )
      UPPER = LSAME( UPLO, 'U' )
      INFO = 0

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MOD(N,2).EQ.1 ) THEN
*     We need an even-dimensional matrix
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( NB.LT.3 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKTF3', -INFO )
         RETURN
      END IF

*     Quick return if possible
      IF( N .EQ. 0 ) RETURN

      IF( UPPER ) THEN
        INFO = -1
        CALL XERBLA( 'DSKTF3', -INFO )
        RETURN
      ELSE
*     Factorize A as L * T * L^T using the lower triangle of A

         IF ( INIT ) IPIV( 1 ) = 1

         DO 20 K=1, MIN(N-1, NB), 1

*     Pivoting for one column
*     Find the pivot
            KP = K + IDAMAX(N-K, A( K+1, K ), 1)
            COLMAX = ABS( A( KP, K ) )

            IF( COLMAX.EQ.ZERO ) THEN
*     The column is completely zero - do nothing
               IF( INFO.EQ.0 ) THEN
                  INFO = K
               END IF
               KP = K+1
            END IF

*     swap rows and columns K+1 and IMAX in the
*     full matrix A(1:N,1:N)
            KK = K+1

            IF( KP .NE. KK ) THEN
               IF( KP.LT.N ) THEN
                  CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ),1 )
               END IF

               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1,
     $                              A( KP, KK+1 ), LDA )

               CALL DSWAP( K+GLOBALK-1, A( KK, 2-GLOBALK), LDA,
     $                                  A( KP, 2-GLOBALK), LDA)

               CALL DSCAL(KP-KK, -ONE, A(KK+1, KK), 1)
               CALL DSCAL(KP-KK-1, -ONE, A(KP, KK+1), LDA)
            END IF
*     Store Pivot
            IPIV( K+1 ) = KP

            CALL DSCAL(N-K-1, ONE/A( K+1, K ), A( K+2, K ), 1)

***********************************************************************
*     Left-looking updates
***********************************************************************
            IF( (K.GT.1 .OR. .NOT.INIT) .AND. K.LE.N-2 ) THEN
               IF( NB.GT.1 .AND. K.EQ.NB ) THEN
*     .. reusing KP as loop variable ..
                  IF( INIT ) THEN
                     DO KP=1, N-K-1, 1
                        W( 1, KP ) = -A( 3, 2 ) * A( K+KP, 2 )
                        DO KK=2, K-2, 1
                           W( KK, KP ) = A( KK+1, KK ) * A( K+KP, KK-1 )
     $                                  -A(KK+2, KK+1) * A( K+KP, KK+1 )
                        END DO
                        IF( KP.EQ.1 ) THEN
                           W( K-1, KP ) = A( K, K-1 ) * A( K+KP, K-2 )
     $                                   -A( K+1, K ) * ONE
                        ELSE
                           W( K-1, KP ) = A( K, K-1 ) * A( K+KP, K-2 )
     $                                   -A( K+1, K ) * A( K+KP, K )
                        END IF
                        W( K, KP ) = A( K+1, K ) * A( K+KP, K-1 )
                     END DO

                     CALL DGEMMT( 'L', 'N', 'N', N-K-1, K, -ONE,
     $                            A( K+2, 1 ), LDA, W, LDW, ONE,
     $                            A( K+2, K+1 ), LDA )
                  ELSE
                     DO KP=1, N-K-1, 1
                        W( 1, KP ) = -A( 2, 1 ) * A( K+KP, 1 )
                        DO KK=1, K-2, 1
                           W(KK+1, KP) = A( KK+1, KK ) * A( K+KP, KK-1 )
     $                                  -A(KK+2, KK+1) * A( K+KP, KK+1 )
                        END DO
                        IF( KP.EQ.1 ) THEN
                           W( K, KP ) = A( K, K-1 ) * A( K+KP, K-2 )
     $                                 -A( K+1, K ) * ONE
                        ELSE
                           W( K, KP ) = A( K, K-1 ) * A( K+KP, K-2 )
     $                                 -A( K+1, K ) * A( K+KP, K )
                        END IF
                        W( K+1, KP ) = A( K+1, K ) * A( K+KP, K-1 )
                     END DO

                     CALL DGEMMT( 'L', 'N', 'N', N-K-1, K+1, -ONE,
     $                            A( K+2, 0 ), LDA, W, LDW, ONE,
     $                            A( K+2, K+1 ), LDA )
                  END IF

               ELSE
                  IF( K.EQ.1 ) THEN
                     W( 1, 1 ) = -A( 2, 1 )
                     W( 2, 1 ) =  A( 2, 1 ) * A( K+1, 0 )
                  ELSE IF( K.EQ.2 .AND. INIT ) THEN
                     W( 1, 1 ) = -A( 3, 2 )
                     W( 2, 1 ) =  A( 3, 2 ) * A( K+1, 1 )
                  ELSE IF( INIT ) THEN
                     W( 1, 1 ) = -A( 3, 2 ) * A( K+1, 2 )
                     DO KK=2, K-2, 1
                        W( KK, 1 ) = A( KK+1, KK ) * A( K+1, KK-1 )
     $                              -A( KK+2, KK+1 ) * A( K+1, KK+1 )
                     END DO
                     W( K-1, 1 ) = A( K, K-1 ) * A( K+1, K-2 )
     $                            -A( K+1, K ) * ONE
                     W( K, 1 ) = A( K+1, K ) * A( K+1, K-1 )
                  ELSE
                     W( 1, 1 ) = -A( 2, 1 ) * A( K+1, 1 )
                     DO KK=1, K-2, 1
                        W( KK+1, 1 ) = A( KK+1, KK ) * A( K+1, KK-1 )
     $                                -A( KK+2, KK+1 ) * A( K+1, KK+1 )
                     END DO
                     W( K, 1 ) = A( K, K-1 ) * A( K+1, K-2 )
     $                          -A( K+1, K ) * ONE
                     W( K+1, 1 ) = A( K+1, K ) * A( K+1, K-1 )
                  END IF

                  IF( INIT ) THEN
                     CALL DGEMV( 'N', N-K-1, K, -ONE, A( K+2, 1 ), LDA,
     $                            W, 1, ONE, A( K+2, K+1 ), 1 )
                  ELSE
                     CALL DGEMV( 'N', N-K-1, K+1,-ONE, A( K+2, 0 ), LDA,
     $                            W, 1, ONE, A( K+2, K+1 ), 1 )
                  END IF
               END IF
            END IF
***********************************************************************

 20      CONTINUE

      END IF

      END

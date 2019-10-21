      SUBROUTINE gemm(N,A,B,C)
      INTEGER N
      REAL*8  A(N,N), B(N,N), C(N,N)
      
      INTEGER I,J,K

      DO 10,J=1,N
         DO 20,K=1,N
            DO 30,I=1,N
               C(I,J) = C(I,J)+A(I,K)*B(K,J)
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE

      END

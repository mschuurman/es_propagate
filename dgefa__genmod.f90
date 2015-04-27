        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar  4 15:33:10 2014
        MODULE DGEFA__genmod
          INTERFACE 
            SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,N)
              INTEGER(KIND=4) :: IPVT(N)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEFA
          END INTERFACE 
        END MODULE DGEFA__genmod

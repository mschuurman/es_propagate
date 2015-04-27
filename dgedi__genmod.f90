        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar  4 15:33:10 2014
        MODULE DGEDI__genmod
          INTERFACE 
            SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,N)
              INTEGER(KIND=4) :: IPVT(N)
              REAL(KIND=8) :: DET(2)
              REAL(KIND=8) :: WORK(N)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DGEDI
          END INTERFACE 
        END MODULE DGEDI__genmod

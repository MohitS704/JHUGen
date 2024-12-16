C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     -(P(-1,3)*Gamma5(-2,1)*Gamma(-1,2,-3)*Gamma(3,-3,-2)) +
C      P(-1,3)*Gamma5(-2,1)*Gamma(-1,-3,-2)*Gamma(3,2,-3)
C     
      SUBROUTINE FFV2_1(F2, V3, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(6)
      COMPLEX*16 F2(*)
      REAL*8 M1
      REAL*8 P1(0:3)
      REAL*8 P3(0:3)
      COMPLEX*16 V3(*)
      REAL*8 W1
      COMPLEX*16 DENOM
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      F1(1) = +F2(1)+V3(1)
      F1(2) = +F2(2)+V3(2)
      P1(0) = -DBLE(F1(1))
      P1(1) = -DBLE(F1(2))
      P1(2) = -DIMAG(F1(2))
      P1(3) = -DIMAG(F1(1))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(3)= DENOM*2D0 * CI*(F2(5)*(P1(1)*(P3(0)*(-V3(4)+CI*(V3(5)))
     $ +(P3(1)*(V3(3)+V3(6))+(P3(2)*(-1D0)*(+CI*(V3(3)+V3(6)))+P3(3)*(
     $ -V3(4)+CI*(V3(5))))))+(P1(2)*(P3(0)*(-1D0)*(+CI*(V3(4))+V3(5))
     $ +(P3(1)*(+CI*(V3(3)+V3(6)))+(P3(2)*(V3(3)+V3(6))-P3(3)*(+CI
     $ *(V3(4))+V3(5)))))+(P1(0)*(-P3(0)*V3(6)-CI*(P3(1)*V3(5))+CI
     $ *(P3(2)*V3(4))+P3(3)*V3(3))+P1(3)*(-P3(0)*V3(6)-CI*(P3(1)*V3(5))
     $ +CI*(P3(2)*V3(4))+P3(3)*V3(3)))))+(F2(6)*(P1(0)*(P3(0)*(-1D0)
     $ *(V3(4)+CI*(V3(5)))+(P3(1)*(V3(3)-V3(6))+(P3(2)*(+CI*(V3(3))-CI
     $ *(V3(6)))+P3(3)*(V3(4)+CI*(V3(5))))))+(P1(3)*(P3(0)*(-1D0)
     $ *(V3(4)+CI*(V3(5)))+(P3(1)*(V3(3)-V3(6))+(P3(2)*(+CI*(V3(3))-CI
     $ *(V3(6)))+P3(3)*(V3(4)+CI*(V3(5))))))+(P1(1)*(P3(0)*V3(6)+CI
     $ *(P3(1)*V3(5))-CI*(P3(2)*V3(4))-P3(3)*V3(3))+P1(2)*(+CI*(P3(0)
     $ *V3(6))-P3(1)*V3(5)+P3(2)*V3(4)-CI*(P3(3)*V3(3))))))+M1*(F2(4)
     $ *(P3(0)*(V3(4)+CI*(V3(5)))+(P3(1)*(-1D0)*(V3(3)+V3(6))+(P3(2)*(
     $ -1D0)*(+CI*(V3(3)+V3(6)))+P3(3)*(V3(4)+CI*(V3(5))))))+F2(3)
     $ *(P3(0)*V3(6)-CI*(P3(1)*V3(5))+CI*(P3(2)*V3(4))-P3(3)*V3(3)))))
      F1(4)= DENOM*2D0 * CI*(F2(5)*(P1(0)*(P3(0)*(-V3(4)+CI*(V3(5)))
     $ +(P3(1)*(V3(3)+V3(6))+(P3(2)*(-1D0)*(+CI*(V3(3)+V3(6)))+P3(3)*(
     $ -V3(4)+CI*(V3(5))))))+(P1(3)*(P3(0)*(V3(4)-CI*(V3(5)))+(P3(1)*(
     $ -1D0)*(V3(3)+V3(6))+(P3(2)*(+CI*(V3(3)+V3(6)))+P3(3)*(V3(4)-CI
     $ *(V3(5))))))+(P1(1)*(-P3(0)*V3(6)-CI*(P3(1)*V3(5))+CI*(P3(2)
     $ *V3(4))+P3(3)*V3(3))+P1(2)*(+CI*(P3(0)*V3(6))-P3(1)*V3(5)+P3(2)
     $ *V3(4)-CI*(P3(3)*V3(3))))))+(F2(6)*(P1(1)*(P3(0)*(-1D0)*(V3(4)
     $ +CI*(V3(5)))+(P3(1)*(V3(3)-V3(6))+(P3(2)*(+CI*(V3(3))-CI*(V3(6))
     $ )+P3(3)*(V3(4)+CI*(V3(5))))))+(P1(2)*(P3(0)*(+CI*(V3(4))-V3(5))
     $ +(P3(1)*(-CI*(V3(3))+CI*(V3(6)))+(P3(2)*(V3(3)-V3(6))+P3(3)*(
     $ -CI*(V3(4))+V3(5)))))+(P1(0)*(P3(0)*V3(6)+CI*(P3(1)*V3(5))-CI
     $ *(P3(2)*V3(4))-P3(3)*V3(3))+P1(3)*(-P3(0)*V3(6)-CI*(P3(1)*V3(5))
     $ +CI*(P3(2)*V3(4))+P3(3)*V3(3)))))+M1*(F2(3)*(P3(0)*(V3(4)-CI
     $ *(V3(5)))+(P3(1)*(-V3(3)+V3(6))+(P3(2)*(+CI*(V3(3))-CI*(V3(6)))
     $ +P3(3)*(-V3(4)+CI*(V3(5))))))+F2(4)*(-P3(0)*V3(6)+CI*(P3(1)
     $ *V3(5))-CI*(P3(2)*V3(4))+P3(3)*V3(3)))))
      F1(5)= DENOM*(-2D0 * CI)*(F2(3)*(P1(1)*(P3(0)*(-V3(4)+CI*(V3(5)))
     $ +(P3(1)*(V3(3)-V3(6))+(P3(2)*(-CI*(V3(3))+CI*(V3(6)))+P3(3)
     $ *(V3(4)-CI*(V3(5))))))+(P1(2)*(P3(0)*(-1D0)*(+CI*(V3(4))+V3(5))
     $ +(P3(1)*(+CI*(V3(3))-CI*(V3(6)))+(P3(2)*(V3(3)-V3(6))+P3(3)*(
     $ +CI*(V3(4))+V3(5)))))+(P1(0)*(P3(0)*V3(6)-CI*(P3(1)*V3(5))+CI
     $ *(P3(2)*V3(4))-P3(3)*V3(3))+P1(3)*(-P3(0)*V3(6)+CI*(P3(1)*V3(5))
     $ -CI*(P3(2)*V3(4))+P3(3)*V3(3)))))+(F2(4)*(P1(0)*(P3(0)*(V3(4)
     $ +CI*(V3(5)))+(P3(1)*(-1D0)*(V3(3)+V3(6))+(P3(2)*(-1D0)*(+CI
     $ *(V3(3)+V3(6)))+P3(3)*(V3(4)+CI*(V3(5))))))+(P1(3)*(P3(0)*(-1D0)
     $ *(V3(4)+CI*(V3(5)))+(P3(1)*(V3(3)+V3(6))+(P3(2)*(+CI*(V3(3)
     $ +V3(6)))-P3(3)*(V3(4)+CI*(V3(5))))))+(P1(1)*(P3(0)*V3(6)-CI
     $ *(P3(1)*V3(5))+CI*(P3(2)*V3(4))-P3(3)*V3(3))+P1(2)*(+CI*(P3(0)
     $ *V3(6))+P3(1)*V3(5)-P3(2)*V3(4)-CI*(P3(3)*V3(3))))))+M1*(F2(6)
     $ *(P3(0)*(-1D0)*(V3(4)+CI*(V3(5)))+(P3(1)*(V3(3)-V3(6))+(P3(2)*(
     $ +CI*(V3(3))-CI*(V3(6)))+P3(3)*(V3(4)+CI*(V3(5))))))+F2(5)*(
     $ -P3(0)*V3(6)-CI*(P3(1)*V3(5))+CI*(P3(2)*V3(4))+P3(3)*V3(3)))))
      F1(6)= DENOM*(-2D0 * CI)*(F2(3)*(P1(0)*(P3(0)*(V3(4)-CI*(V3(5)))
     $ +(P3(1)*(-V3(3)+V3(6))+(P3(2)*(+CI*(V3(3))-CI*(V3(6)))+P3(3)*(
     $ -V3(4)+CI*(V3(5))))))+(P1(3)*(P3(0)*(V3(4)-CI*(V3(5)))+(P3(1)*(
     $ -V3(3)+V3(6))+(P3(2)*(+CI*(V3(3))-CI*(V3(6)))+P3(3)*(-V3(4)+CI
     $ *(V3(5))))))+(P1(1)*(-P3(0)*V3(6)+CI*(P3(1)*V3(5))-CI*(P3(2)
     $ *V3(4))+P3(3)*V3(3))+P1(2)*(+CI*(P3(0)*V3(6))+P3(1)*V3(5)-P3(2)
     $ *V3(4)-CI*(P3(3)*V3(3))))))+(F2(4)*(P1(1)*(P3(0)*(-1D0)*(V3(4)
     $ +CI*(V3(5)))+(P3(1)*(V3(3)+V3(6))+(P3(2)*(+CI*(V3(3)+V3(6)))
     $ -P3(3)*(V3(4)+CI*(V3(5))))))+(P1(2)*(P3(0)*(+CI*(V3(4))-V3(5))
     $ +(P3(1)*(-1D0)*(+CI*(V3(3)+V3(6)))+(P3(2)*(V3(3)+V3(6))+P3(3)*(
     $ +CI*(V3(4))-V3(5)))))+(P1(0)*(-P3(0)*V3(6)+CI*(P3(1)*V3(5))-CI
     $ *(P3(2)*V3(4))+P3(3)*V3(3))+P1(3)*(-P3(0)*V3(6)+CI*(P3(1)*V3(5))
     $ -CI*(P3(2)*V3(4))+P3(3)*V3(3)))))+M1*(F2(5)*(P3(0)*(-V3(4)+CI
     $ *(V3(5)))+(P3(1)*(V3(3)+V3(6))+(P3(2)*(-1D0)*(+CI*(V3(3)+V3(6)))
     $ +P3(3)*(-V3(4)+CI*(V3(5))))))+F2(6)*(P3(0)*V3(6)+CI*(P3(1)*V3(5)
     $ )-CI*(P3(2)*V3(4))-P3(3)*V3(3)))))
      END



MODULE ROTATION
  ! Utilities for working with rotations of rigid objects
  ! especially in terms of quaternions
  USE GENUTIL
  USE MT19937
  IMPLICIT NONE

CONTAINS
  SUBROUTINE XYZANG2QUAT(ANG,QUAT)
    ! convert an *extrinsic* rotation around (fixed) axes x, y, z by the
    ! three angles given in ANG, to a quaternion representation
    ! input ANG is the x-y-z extrinsic Tait-Bryan angles
    ! NOTE: left-muliplying by a quaternion gives an extrinsic rotation
    
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: ANG(3)
    DOUBLE PRECISION, INTENT(OUT) :: QUAT(4)
    DOUBLE PRECISION :: QX(4), QY(4), QZ(4), QYX(4)   
    
    ! Rotation around original x
    CALL ROTQUAT(ANG(1),(/1D0,0D0,0D0/),QX)
    ! rotation around original y
    CALL ROTQUAT(ANG(2),(/0D0,1D0,0D0/),QY)
    ! rotation around original z
    CALL ROTQUAT(ANG(3),(/0D0,0D0,1D0/),QZ)   
    
    QYX = QUATMULT(QY,QX)    
    QUAT = QUATMULT(QZ,QYX)
  END SUBROUTINE XYZANG2QUAT
  
  SUBROUTINE ROTQUAT(THETA,AX,QV,DT)
    ! get the quaternion (QV) corresponding to rotation around axis AX by angle THETA
    ! optionally, also get the derivative wrt theta
    ! WARNING: AX assumed to be normalized; does not check for this!

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: THETA, AX(3)
    DOUBLE PRECISION, INTENT(OUT) :: QV(4)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DT(4)
    DOUBLE PRECISION :: CT, ST

    CT = COS(THETA/2); ST = SIN(THETA/2);
    QV(1) = CT
    QV(2:4) = ST*AX

    IF (PRESENT(DT)) THEN
       DT(1) = -ST/2
       DT(2:4) = CT/2*AX
    ENDIF
    
  END SUBROUTINE ROTQUAT
  
  SUBROUTINE QUAT2ROTMAT(Q,MAT,DMAT)
    ! convert a quaternion object to a rotation matrix and (optionally) return derivatives
    ! input: Q -> quaternion expressed as 4-vector
    ! output: MAT -> rotation matrix
    ! DMAT(i,j,k) -> derivative of MAT(i,j) wrt quaternion coordinate k
    ! NOTE: no normalization
    IMPLICIT NONE
    DOUBLE PRECISION :: Q(4)
    DOUBLE PRECISION, INTENT(OUT) :: MAT(3,3)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: DMAT(3,3,4)
    DOUBLE PRECISION :: A,B,C,D,AA,BB,CC,DD,AB,AC,AD,BC,BD,CD

    A = Q(1); B = Q(2); C = Q(3); D = Q(4)
    AA = A*A; BB = B*B; CC = C*C; DD = D*D
    AB = 2*A*B; AC = 2*A*C; AD = 2*A*D
    BC = 2*B*C; BD = 2*B*D
    CD = 2*C*D

    MAT(1,:) = (/AA+BB-CC-DD,BC-AD,AC+BD/)
    MAT(2,:) = (/AD+BC,AA-BB+CC-DD,CD-AB/)
    MAT(3,:) = (/BD-AC,AB+CD,AA-BB-CC+DD/)

    IF (PRESENT(DMAT)) THEN
       dMAT(1,:,1) = (/A,-D,C/)
       dMAT(1,:,2) = (/B,C,D/)
       dMAT(1,:,3) = (/-C,B,A/)
       dMAT(1,:,4) = (/-D,-A,B/)
       dMAT(2,:,1) = (/D,A,-B/)
       dMAT(2,:,2) = (/C,-B,-A/)
       dMAT(2,:,3) = (/B,C,D/)
       dMAT(2,:,4) = (/A,-D,C/)
       dMAT(3,:,1) = (/-C,B,A/)
       dMAT(3,:,2) = (/D,A,-B/)
       dMAT(3,:,3) = (/-A,D,-C/)
       dMAT(3,:,4) = (/B,C,D/)
       dMAT = dMAT*2
    ENDIF
  END SUBROUTINE QUAT2ROTMAT
    
  SUBROUTINE UNIFORMRANDQUAT(QUAT)
    ! Generate a random quaternion rotation,
    ! uniformly distributed in magnitude and direction

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: QUAT(4)
    INTEGER :: I
    DOUBLE PRECISION :: AX(3), THETA
    
    ! generate a uniformly random axis
    DO I = 1,3
       AX(I) = GRND()
    END DO
    CALL NORMALIZE(AX)

    ! generate a uniformly random rotation
    THETA = GRND()*2*PI

    ! convert to a quaternion
    CALL ROTQUAT(THETA,AX,QUAT)
  END SUBROUTINE UNIFORMRANDQUAT

  FUNCTION QUATMULT(P,Q)
    ! Product of 2 quaternions
    ! Note: for a body with orientation Q, this performs rotation P around extrinsic (lab-frame) axes.
    ! right-multiply by a quaternion to rotate around intrinsic axes

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: P(4),Q(4)
    DOUBLE PRECISION :: QUATMULT(4)

    QUATMULT(1) = P(1)*Q(1) - P(2)*Q(2) - P(3)*Q(3) - P(4)*Q(4)
    QUATMULT(2) = P(1)*Q(2) + P(2)*Q(1) + P(3)*Q(4) - P(4)*Q(3)
    QUATMULT(3) = P(1)*Q(3) - P(2)*Q(4) + P(3)*Q(1) + P(4)*Q(2)
    QUATMULT(4) = P(1)*Q(4) + P(2)*Q(3) - P(3)*Q(2) + P(4)*Q(1)
    
  END FUNCTION QUATMULT
END MODULE ROTATION

PROGRAM MAIN
  ! various subroutines for testing code
  USE KEYS
  USE PARTICLEUTIL
  USE GENUTIL
  IMPLICIT NONE
  TYPE(PARTICLEGROUP), POINTER :: PARTGROUP(:)
  TYPE(PARTICLEGROUP), POINTER :: PARTP

  DOUBLE PRECISION, ALLOCATABLE :: INITPOS(:,:)
  DOUBLE PRECISION :: TESTMAT(3,3)
  INTEGER :: I,J, PC = 1
  DOUBLE PRECISION :: LEBPTS(6,3)
  DOUBLE PRECISION :: TMP(3)

	!parses free format input file
!  CALL READKEY
!  
!  ALLOCATE(PARTGROUP(NTRIALS))
!  ALLOCATE(INITPOS(3,NPART))

!	DO PC = 1,NTRIALS
!	  PARTP=>PARTGROUP(PC)
!	  CALL SETUPPARTICLEGROUP(PARTP,NPART)
!	  CALL SETPARTPARAMS(PARTP)
!	  INITPOS(:,1) = (/0D0,0D0,DBLE(PC)/)
!	  CALL INITIALIZEPARTICLES(PARTP, INITPOS(:,1:NPART),RANDORIENT=.true.)
!  END DO

!	DO PC = 1,NTRIALS
!		PARTP=>PARTGROUP(PC)
!		PRINT*, PARTP%POS(:,1)
!	END DO
!  INITPOS(:,1) = (/0,0,0/)
!  INITPOS(:,2) = (/0,5,5/)
!  
!  print*, sqrt(sum((initpos(:,1)-initpos(:,2))**2))

!  CALL INITIALIZEPARTICLES(PARTP, INITPOS(:,1:NPART),RANDORIENT=.true.)

!  SELECT CASE (ACTION)
!  CASE ('BROWNDYN')
!     CALL BROWNDYNDRIVER(PARTP)
!  CASE DEFAULT
!     PRINT*, 'UNKNOWN ACTION:', ACTION
!     STOP 1
!  END SELECT
  
!  PARTP%POS(:,1) = (/1D0,1D0,1D0/)
!  PARTP%POS(:,2) = (/2D0,2D0,2D0/)
!  PRINT*, 'TESTXy: ', PARTP%POS(:,1),PARTP%POS(:,2)
!  CALL TESTPASSBYREF(PARTP%POS(:,1),PARTP%POS(:,2))
!  PRINT*, 'TESTXB: ', PARTP%POS(:,1),PARTP%POS(:,2)

  ! CALL SNAPSHOTPARTICLES(PARTP,SNAPSHOTFILE,.FALSE.)
!  CALL TESTQUATMULT
	TMP = (/1D0,1D0,1D0/)
	NORM2(TMP)
	

  !CALL TESTENERGYDERIVS(CHAINP)
  !CALL TESTBROWNDYN(PARTP)
!	DO PC = 1,NTRIALS
!		PARTP=>PARTGROUP(PC)
!	  CALL CLEANUPPARTICLEGROUP(PARTP)
!  END DO
!  DEALLOCATE(PARTGROUP)

!	OPEN(UNIT=44,FILE="/home/smogre/proj/spheretouchBD/source/lebfiles/lebedev_003.txt",STATUS="OLD",ACTION="READ") 
!	DO I = 1,6
!		READ(UNIT=44,FMT=*) LEBPTS(I,:)
!		PRINT*, LEBPTS(I,:)
!	END DO

CONTAINS
	
	SUBROUTINE TESTPASSBYREF(POS0,POS1)
		IMPLICIT NONE
		DOUBLE PRECISION, INTENT(IN) :: POS0(3)
		DOUBLE PRECISION, INTENT(OUT) :: POS1(3)
		POS1(2) = 100D0
	
	END SUBROUTINE TESTPASSBYREF
	
  SUBROUTINE TESTQUATMULT
    ! test quaternion multiplication to represent rotation
    USE ROTATION
    IMPLICIT NONE
    DOUBLE PRECISION :: Q(4), P(4), R1(4), R2(4), ROTMAT(3,3), ROTMAT2(3,3)
    INTEGER :: I

    ! quaternion for cannonical orientation
    CALL ROTQUAT(PI/2,(/0D0,0D0,1D0/),Q)
!		CALL ROTQUAT(0d0,(/0D0,0D0,1D0/),Q)
    CALL QUAT2ROTMAT(Q,ROTMAT)

	  PRINT*,'Canonical rotation:'
	  PRINT*, 'ROTATION MATRIX:'
    DO I = 1,3
	  	DO J = 1,3
       PRINT*, 'I= ',I,' J = ',J,' MATRIX ELEMENT: ', ROTMAT(I,J)
    	END DO
    ENDDO
!    PRINT*, 'QUATERNION:'
!    PRINT*, Q

    ! quaternions for 2 rotations
    CALL ROTQUAT(PI/2,(/0D0,0D0,1D0/),R1)
    CALL ROTQUAT(PI/2,(/1D0,0D0,0D0/),R2)

    ! p = R2*R1*Q
    P = QUATMULT(R1,Q)
		CALL QUAT2ROTMAT(p,ROTMAT)
    DO I = 1,3
       PRINT*, ROTMAT(I,:)
    ENDDO

    P = QUATMULT(R2,p)

    CALL QUAT2ROTMAT(P,ROTMAT2)
    DO I = 1,3
       PRINT*, ROTMAT2(I,:)
    ENDDO

    PRINT*, 'TESTING XYZ EXTRINSIC ROTATION:'

    CALL XYZANG2QUAT((/PI/2,PI/2,PI/2/),P)
    CALL QUAT2ROTMAT(P,ROTMAT2)
    DO I = 1,3
       PRINT*, ROTMAT2(I,:)
    ENDDO


  END SUBROUTINE TESTQUATMULT

  SUBROUTINE BROWNDYNDRIVER(PARTP)
    ! driver for brownian dynamics calculations
    USE BROWNDYN, ONLY : RUNBROWNDYNSIM
    USE KEYS, ONLY : BDSTEPS, DELT, KT, OUTFILE, BDPRINTEVERY, SNAPSHOTFILE,&
         & SNAPSHOTEVERY, APPENDSNAPSHOTS, DOBROWN, DUMPSNAPSHOTS
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    INTEGER :: SNAPEVERY

    IF (DUMPSNAPSHOTS) THEN
       SNAPEVERY = SNAPSHOTEVERY
    ELSE
       SNAPEVERY = HUGE(1)
    ENDIF

    CALL  RUNBROWNDYNSIM(PARTP, BDSTEPS, DELT,KT,OUTFILE,BDPRINTEVERY,&
         & SNAPSHOTFILE,SNAPEVERY,APPENDSNAPSHOTS,DOBROWN)

  END SUBROUTINE BROWNDYNDRIVER

END PROGRAM MAIN

MODULE BROWNDYN
  ! Utilities for running brownian dynamics simulations (rotational and translational motion)
  ! for a group of spherical particles
  ! ignores hydrodynamic interactions, no direct coupling of rotation / translation
  USE PARTICLEUTIL, ONLY : PARTICLEGROUP, OUTPUTSNAPSHOT, GETFORCESTORQUES, REFLECTPART, CHECKRXN, CHECKTRACKPROX
  USE MT19937, ONLY : RNORM, GRND
  USE ROTATION, ONLY : QUATMULT, XYZANG2QUAT
  IMPLICIT NONE

CONTAINS
  SUBROUTINE RUNBROWNDYNSIM(PARTP, NSTEP, DELT,KT,OUTFILE,PRINTEVERY,&
       & SNAPSHOTFILE,SNAPSHOTEVERY,APPENDSNAPSHOTS,DOBROWN,WHICHREACT, CURTIME, WHICHTRACK, TRACKHITTIME)
    ! run a brownian dynamics simulation for a group of particles
    ! NSTEP: number of steps to simulate
    ! DELT: timestep
    ! kt : temperature for brownian forces/torques
    ! OUTFILE: output file for printing out info (what do we want to keep track of?)
    ! PRINTEVERY: how often to print info into OUTFILE
    ! SNAPSHOTFILE: file for dumping chain configuration snapshots
    ! SNAPSHOTEVERY: how often to dump snapshots
    ! APPENDSNAPSHOTS: append snapshots to existing files (otherwise, start from scratch and append as we go)
    ! DOBROWN: include brownian forces
		
    USE KEYS, ONLY:DOMRAD,DOMLEN,KSP,KTETH,ATTLEN,NTRIALS, STOPHITTRACK
    
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    INTEGER, INTENT(IN) :: NSTEP
    DOUBLE PRECISION, INTENT(IN) :: DELT, kt
    CHARACTER(LEN=*), INTENT(IN) :: SNAPSHOTFILE, OUTFILE
    INTEGER, INTENT(IN) :: PRINTEVERY, SNAPSHOTEVERY
    LOGICAL, INTENT(IN) :: APPENDSNAPSHOTS, DOBROWN
    INTEGER, INTENT(OUT) :: WHICHREACT, WHICHTRACK
    DOUBLE PRECISION, INTENT(OUT) :: CURTIME,TRACKHITTIME
    
    INTEGER :: STEP
    DOUBLE PRECISION :: ENERGY, ENERGY0, FORCES(3,PARTP%NP)
    
    DOUBLE PRECISION, ALLOCATABLE :: INFO(:)

    OPEN(FILE=OUTFILE,UNIT=88,STATUS='UNKNOWN')

    CURTIME = 0D0
    
    IF(PRINTEVERY.NE.0) THEN
	    PRINT'(A27,I7,4ES15.4E2)', 'STEP, TIME, PARTICLE 1 POS:', 0, CURTIME, PARTP%POS(:,1)    
    END IF

    ! Initial snapshot
    ALLOCATE(INFO(11+PARTP%NWALK+PARTP%NTETH))
    INFO = (/CURTIME,PARTP%FRIC(1), PARTP%FRICR(1), PARTP%SRAD(1),PARTP%CRAD(1),DFLOAT(PARTP%WALKINDS),&
    			DFLOAT(PARTP%TETHINDS), PARTP%VEL,DOMRAD,DOMLEN,KSP,KTETH,ATTLEN/)
	  
    CALL OUTPUTSNAPSHOT(PARTP,SNAPSHOTFILE,INFO,APPENDSNAPSHOTS)

    WHICHTRACK = 0
    TRACKHITTIME = 0D0
    WHICHREACT = 0
    
    DO STEP = 1,NSTEP
       
       ! Check for proximity to track
       IF (WHICHTRACK.EQ.0) THEN
          CALL CHECKTRACKPROX(PARTP,WHICHTRACK)
          IF (WHICHTRACK.GT.0) THEN ! first time hitting track
             TRACKHITTIME = CURTIME
             IF (STOPHITTRACK) RETURN ! stop after first hitting a track
          ENDIF
       ENDIF

       ! check if any reactions occured
       IF (PARTP%NRXN.GT.0) THEN
          CALL CHECKRXN(PARTP, WHICHREACT)
          IF (WHICHREACT.GT.0) THEN
             ! rxn has occured. Stop running this trial             
             RETURN
          ENDIF
       ENDIF
       
       !CALL LANGEVINSTEPRK4(CHAINP,DELT,ENERGY,DOBROWN)
       !CALL EULERSTEP(PARTP,DELT,KT,DOBROWN)
       CALL RK4STEP(PARTP,DELT,KT,DOBROWN)       
       
       
       CURTIME = CURTIME+DELT

    	IF ((PRINTEVERY.NE.0).AND.MOD(STEP,PRINTEVERY).EQ.0) THEN
      	! Print before starting a new iteration
        ! Output information on particle status
        PRINT'(A27,I8,4ES15.4E2)', 'STEP, TIME, PARTICLE 1 POS:', STEP, CURTIME, PARTP%POS(:,1)
      ENDIF

      IF (MOD(STEP,SNAPSHOTEVERY).EQ.0) THEN
         ! output snapshot
         INFO(1) = CURTIME
         CALL OUTPUTSNAPSHOT(PARTP,SNAPSHOTFILE,INFO,.TRUE.)
      ENDIF
    ENDDO

    CLOSE(88)
  END SUBROUTINE RUNBROWNDYNSIM

  SUBROUTINE EULERSTEP(PARTP,DELT,KT,DOBROWN)
    ! Simple forward euler step for the particle group
    ! INPUT:
    ! PARTP = object describing group of particles
    ! DELT = timestep
    ! KT = temperature for brownian forces/torques
    ! DOBROWN = include brownian forces?
    ! OUTPUT: %pos and %quat fields in partp are updated after the step
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION, INTENT(IN) :: DELT, KT
    LOGICAL, INTENT(IN) :: DOBROWN
    DOUBLE PRECISION :: STDF(PARTP%NP), STDTQ(PARTP%NP)
    INTEGER :: PC, I
    DOUBLE PRECISION :: FBROWN(3,PARTP%NP), TQBROWN(3,PARTP%NP), TORQUEBODY(3,PARTP%NP)
    DOUBLE PRECISION :: FORCE(3,PARTP%NP), TORQUE(3,PARTP%NP)
    DOUBLE PRECISION :: ANGROT(3), QROT(4)

    ! get standard deviations for Brownian force and torque
    STDF = SQRT(2*KT*PARTP%FRIC/DELT)
    STDTQ = SQRT(2*KT*PARTP%FRICR/DELT)

    ! get the brownian forces (in lab frame)
    IF (DOBROWN) THEN
       DO PC = 1,PARTP%NP
          DO I = 1,3
             FBROWN(I,PC) = RNORM()*STDF(pc)
          ENDDO
       END DO
    ELSE
       FBROWN = 0D0
    ENDIF

    ! get the brownian torques (in particle frame)
    IF (DOBROWN) THEN
       DO PC = 1,PARTP%NP
          DO I = 1,3
             TQBROWN(I,PC) = RNORM()*STDTQ(PC)
          ENDDO
       END DO
    ELSE
       TQBROWN = 0D0
    ENDIF

    ! get the deterministic forces and torques (lab frame of ref)
    CALL GETFORCESTORQUES(PARTP,FORCE,TORQUE)

    ! torques from lab frame to particle frame
    ! ------------ NOT IMPLEMENTED YET ------
     TORQUEBODY = 0D0

    ! propagate particles
    ! implement reflecting/periodic boundary conditions
    DO PC = 1,PARTP%NP
        IF(ANY(PARTP%WALKINDS.EQ.PC)) THEN
            PARTP%POS(:,PC) = PARTP%POS(:,PC) + (/0d0,0d0,PARTP%VEL/)*DELT
        ELSE
            ! get rotation quaternion from torque
            ANGROT = (TQBROWN(:,PC)+TORQUEBODY(:,PC))/PARTP%FRICR(PC)*DELT ! xyz rotation angles
            CALL XYZANG2QUAT(ANGROT,QROT)
            ! carry out particle rotation (in body frame -> right multiply by quaternion)
            PARTP%QUAT(:,PC) = QUATMULT(PARTP%QUAT(:,PC),QROT)

            ! re-normalize quaternion for numerical stability
            PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))

            ! shift particle (in lab frame)
            PARTP%POS(:,PC) = PARTP%POS(:,PC) + (FORCE(:,PC)+FBROWN(:,PC))/PARTP%FRIC(PC)*DELT
        END IF
    ENDDO
  END SUBROUTINE EULERSTEP
  
  SUBROUTINE RK4STEP(PARTP,DELT,KT,DOBROWN)
    ! RK4 step for the particle group
    ! INPUT:
    ! PARTP = object describing group of particles
    ! DELT = timestep
    ! KT = temperature for brownian forces/torques
    ! DOBROWN = include brownian forces?
    ! OUTPUT: %pos and %quat fields in partp are updated after the step
    USE KEYS, ONLY: DOMRAD,DOMLEN,PERBOUND,NTRACK,ATTLEN,NPART,STERRAD
    USE GENUTIL, ONLY: PI
    USE ROTATION, ONLY: ROTQUAT
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION, INTENT(IN) :: DELT, KT
    LOGICAL, INTENT(IN) :: DOBROWN
    DOUBLE PRECISION :: STDF(PARTP%NP), STDTQ(PARTP%NP)
    INTEGER :: PC, I
    DOUBLE PRECISION :: FBROWN(3,PARTP%NP), TQBROWN(3,PARTP%NP)
    DOUBLE PRECISION :: FORCE(3,PARTP%NP), TORQUE(3,PARTP%NP)
    DOUBLE PRECISION :: ANGROT(3), QROT(4)
    DOUBLE PRECISION, DIMENSION(3) :: K1ANG,K2ANG,K3ANG,K4ANG,K1POS,K2POS,K3POS,K4POS
    DOUBLE PRECISION :: POS0(3), QUAT0(4)
    DOUBLE PRECISION :: DIST,TMP(3)
    DOUBLE PRECISION :: X0,Y0,XT,YT, TH, MULT = 1D0
    LOGICAL :: DISTFLAG
    INTEGER :: IND,LC,CNT
    

    ! get standard deviations for Brownian force and torque
    STDF = SQRT(2*KT*PARTP%FRIC/DELT)
    STDTQ = SQRT(2*KT*PARTP%FRICR/DELT)

    ! get the brownian forces (in lab frame)
    IF (DOBROWN) THEN
       DO PC = 1,PARTP%NP
          DO I = 1,3
             FBROWN(I,PC) = RNORM()*STDF(pc)
          ENDDO
       END DO
    ELSE
       FBROWN = 0D0
    ENDIF

    ! get the brownian torques (in particle frame)
    IF (DOBROWN) THEN
       DO PC = 1,PARTP%NP
          DO I = 1,3
             TQBROWN(I,PC) = RNORM()*STDTQ(PC)
          ENDDO
       END DO
    ELSE
       TQBROWN = 0D0
    ENDIF

    !MOVE SPRING
    DO PC = 1,PARTP%NP
		  IF(ANY(PARTP%WALKINDS.EQ.PC)) THEN
		  	PARTP%SPRPOS(3,PC) = PARTP%SPRPOS(3,PC)+PARTP%VEL*DELT
		  END IF
	  END DO

    !PROPAGATE PARTICLES
    DO PC = 1,PARTP%NP
  	  POS0 = PARTP%POS(:,PC)
	    QUAT0 = PARTP%QUAT(:,PC)
	    ATTLEN = STERRAD(PC)

            
        !-------------1ST RK4 STEP----------------
        
      !GET TORQUE AND FORCES FOR CURRENT STEP
      CALL GETFORCESTORQUES(PARTP,FORCE,TORQUE)
        
      !GET POSITION AND ANGLE FOR CURRENT STEP
      K1POS = (FORCE(:,PC)+FBROWN(:,PC))/PARTP%FRIC(PC)
      K1ANG = (TQBROWN(:,PC) + TORQUE(:,PC))/PARTP%FRICR(PC)
     
      !UPDATE POSITION AND ANGLE FOR NEXT STEP
      PARTP%POS(:,PC) = POS0 + DELT/2*K1POS

      CALL XYZANG2QUAT(DELT/2*K1ANG,QROT)
      PARTP%QUAT(:,PC) = QUATMULT(QUAT0,QROT)
      
      ! re-normalize quaternion for numerical stability
      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
      
      !-------------2ND RK4 STEP----------------
      
      !GET TORQUE AND FORCES FOR CURRENT STEP
      CALL GETFORCESTORQUES(PARTP,FORCE,TORQUE)
      
      !GET POSITION AND ANGLE FOR CURRENT STEP
      K2POS = (FORCE(:,PC)+FBROWN(:,PC))/PARTP%FRIC(PC)
      K2ANG = (TQBROWN(:,PC) + TORQUE(:,PC))/PARTP%FRICR(PC)
      
      !UPDATE POSITION AND ANGLE FOR NEXT STEP
      PARTP%POS(:,PC) = POS0 + DELT/2*K2POS

      CALL XYZANG2QUAT(DELT/2*K2ANG,QROT)
      PARTP%QUAT(:,PC) = QUATMULT(QUAT0,QROT)
      
      ! re-normalize quaternion for numerical stability
      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
      
      !-------------3RD RK4 STEP----------------
      
      !GET TORQUE AND FORCES FOR CURRENT STEP
      CALL GETFORCESTORQUES(PARTP,FORCE,TORQUE)
      
      !GET POSITION AND ANGLE FOR CURRENT STEP
      K3POS = (FORCE(:,PC)+FBROWN(:,PC))/PARTP%FRIC(PC)
      K3ANG = (TQBROWN(:,PC) + TORQUE(:,PC))/PARTP%FRICR(PC)

      !UPDATE POSITION AND ANGLE FOR NEXT STEP
      PARTP%POS(:,PC) = POS0 + DELT*K3POS

      CALL XYZANG2QUAT(DELT*K3ANG,QROT)
      PARTP%QUAT(:,PC) = QUATMULT(QUAT0,QROT)
      
      ! re-normalize quaternion for numerical stability
      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
      
      !-------------4TH RK4 STEP----------------
      
      !GET TORQUE AND FORCES FOR CURRENT STEP
      CALL GETFORCESTORQUES(PARTP,FORCE,TORQUE)
      
      !GET POSITION AND ANGLE FOR CURRENT STEP
      K4POS = (FORCE(:,PC)+FBROWN(:,PC))/PARTP%FRIC(PC)
      K4ANG = (TQBROWN(:,PC) + TORQUE(:,PC))/PARTP%FRICR(PC)
              
   		!ADD UP RK4 COMPONENTS
      PARTP%POS(:,PC) = POS0 + DELT/6*(K1POS + 2*K2POS + 2*K3POS + K4POS)
      ANGROT = DELT/6*(K1ANG+2*K2ANG+2*K3ANG+K4ANG)
      CALL XYZANG2QUAT(ANGROT,QROT)
      PARTP%QUAT(:,PC) = QUATMULT(QUAT0,QROT)
      
      ! re-normalize quaternion for numerical stability
      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
      
      !REFLECT PARTICLE IF OUT OF BOUNDS
      IF(ANY(PARTP%WALKINDS.EQ.PC)) THEN
      
	      IF((ABS(PARTP%SPRPOS(3,PC)).GT.DOMLEN).AND.PERBOUND.AND.(ABS(PARTP%POS(3,PC)).GT.DOMLEN)) THEN
	      	PARTP%POS(3,PC) = PARTP%POS(3,PC) - SIGN(1D0,PARTP%POS(3,PC))*2*DOMLEN
					PARTP%SPRPOS(3,PC) = PARTP%SPRPOS(3,PC) - SIGN(1D0,PARTP%SPRPOS(3,PC))*2*DOMLEN
					
					! REFLECT PARTICLE BACK TO ANOTHER TRACK
					DISTFLAG = .TRUE.
					MULT = 1D0
					CNT = 0
					DO WHILE (DISTFLAG)
						! CHOOSE MICROTUBULE
						IND = CEILING(GRND()*NTRACK)
						XT = PARTP%TRKPOS(1,IND)
						YT = PARTP%TRKPOS(2,IND)
						
						!ORIENT PARTICLE CENTER RELATIVE TO MICROTUBULE
						TH = 2*PI*GRND()
						X0 = XT+ATTLEN*COS(TH)
						Y0 = YT+ATTLEN*SIN(TH)
						IF((X0**2+Y0**2).GT.(DOMRAD-ATTLEN)) CYCLE
						DISTFLAG = .FALSE.
						
						DO LC = 1,NPART
							IF(LC.EQ.PC) CYCLE
							TMP = PARTP%POS(:,LC)-(/X0,Y0,PARTP%POS(3,PC)/)
							DIST = SQRT(DOT_PRODUCT(TMP,TMP))
							
							IF(DIST.LT.MULT*(STERRAD(PC)+STERRAD(LC))) THEN
								DISTFLAG = .TRUE.
								CNT = CNT+1
								IF(CNT.GT.1E4) THEN
										MULT = 0.9*MULT
										CNT = 0
								END IF
								EXIT
							END IF
							
						END DO
						
					END DO
					PARTP%POS(1,PC) = X0
					PARTP%POS(2,PC) = Y0
					PARTP%SPRPOS(1,PC) = XT
					PARTP%SPRPOS(2,PC) = YT
					CALL ROTQUAT(TH,(/0D0,0D0,1D0/),PARTP%QUAT(:,PC))
		      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
					
				END IF

      ELSE
      	!CALL REFLECTPART(POS0,PARTP%POS(:,PC))
      	IF(PERBOUND.AND.(ABS(PARTP%POS(3,PC)).GT.DOMLEN)) THEN
	      	PARTP%POS(3,PC) = PARTP%POS(3,PC) - SIGN(1D0,PARTP%POS(3,PC))*2*DOMLEN
	      
			    ! re-normalize quaternion for numerical stability
		      PARTP%QUAT(:,PC) = PARTP%QUAT(:,PC)/SQRT(SUM(PARTP%QUAT(:,PC)**2))
				END IF
			END IF
    END DO
  END SUBROUTINE RK4STEP
  
  
  ! SUBROUTINE LANGEVINSTEPRK4(CHAINP,DELT,ENERGY,DOBROWN)
  !   ! propagate chain forward in time, using a fourth-order Runge-Kutta method
  !   ! DELT: timestep
  !   ! ENERGY: returns the final energy of the chain after the step
  !   ! DOBROWN: toggle whether to include brownian forces
  !   USE MT19937, ONLY : RNORM, GRND
  !   USE CHAINUTIL, ONLY : CHAIN, GETCHAINENERGY

  !   IMPLICIT NONE
  !   TYPE(CHAIN), POINTER :: CHAINP
  !   DOUBLE PRECISION, INTENT(IN) :: DELT
  !   DOUBLE PRECISION, INTENT(OUT) :: ENERGY
  !   LOGICAL, INTENT(IN) :: DOBROWN
  !   DOUBLE PRECISION :: S2DT
  !   DOUBLE PRECISION, DIMENSION(3,CHAINP%NPT) :: POS0, K1POS, K2POS, K3POS, &
  !        & K4POS, FORCES, FBROWN
  !   INTEGER :: B, I

  !   S2DT = SQRT(2*CHAINP%KT*CHAINP%FRICT/DELT)

  !   POS0 = CHAINP%POS

  !   ! get the brownian forces
  !   IF (DOBROWN) THEN
  !      DO B = 1,CHAINP%NPT
  !         DO I = 1,3
  !            FBROWN(I,B) = RNORM()*S2DT
  !         ENDDO
  !      END DO
  !   ELSE
  !      FBROWN = 0D0
  !   ENDIF

  !   ! --------- 1ST RK STEP---------------
  !   CALL GETCHAINENERGY(CHAINP,ENERGY,FORCES,.TRUE.)
  !   K1POS = (FBROWN + FORCES)/CHAINP%FRICT
  !   IF (CHAINP%NFIX.GT.0) K1POS(:,CHAINP%FIXBEADS(1:CHAINP%NFIX)) = 0D0

  !   CHAINP%POS = POS0 + DELT/2*K1POS

  !   ! --------- 2ND RK STEP---------------
  !   CALL GETCHAINENERGY(CHAINP,ENERGY,FORCES,.TRUE.)
  !   K2POS = (FBROWN + FORCES)/CHAINP%FRICT
  !   IF (CHAINP%NFIX.GT.0) K2POS(:,CHAINP%FIXBEADS(1:CHAINP%NFIX)) = 0D0

  !   CHAINP%POS = POS0 + DELT/2*K2POS

  !   ! --------- 3RD RK STEP---------------
  !   CALL GETCHAINENERGY(CHAINP,ENERGY,FORCES,.TRUE.)
  !   K3POS = (FBROWN + FORCES)/CHAINP%FRICT
  !   IF (CHAINP%NFIX.GT.0) K3POS(:,CHAINP%FIXBEADS(1:CHAINP%NFIX)) = 0D0

  !   CHAINP%POS = POS0 + DELT*K3POS

  !   ! --------- 3RD RK STEP---------------
  !   CALL GETCHAINENERGY(CHAINP,ENERGY,FORCES,.TRUE.)
  !   K4POS = (FBROWN + FORCES)/CHAINP%FRICT
  !   IF (CHAINP%NFIX.GT.0) K4POS(:,CHAINP%FIXBEADS(1:CHAINP%NFIX)) = 0D0

  !   CHAINP%POS = POS0 + DELT/6*(K1POS + 2*K2POS + 2*K3POS + K4POS)

  !   ! calculate final energy
  !   CALL GETCHAINENERGY(CHAINP,ENERGY,FORCES,.false.)
  ! END SUBROUTINE LANGEVINSTEPR4K
END MODULE BROWNDYN

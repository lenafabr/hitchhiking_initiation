MODULE PARTICLEUTIL
  ! utilities for defining parameters and states of groups of particles

  TYPE PARTICLEGROUP
     INTEGER :: NP ! number of particles in the group

     ! particle position
     ! first index is x,y,z, second is particle number
     DOUBLE PRECISION, POINTER :: POS(:,:)
     ! particle orientation (as quaternion)
     ! first index is w,x,y,z, second is particle number
     DOUBLE PRECISION, POINTER :: QUAT(:,:)
     ! friction coefficients for the particles
     DOUBLE PRECISION, POINTER :: FRIC(:), FRICR(:)
     ! steric size of particles
     ! and contact radius of particles
     DOUBLE PRECISION, POINTER :: SRAD(:), CRAD(:)

     ! arrays have been allocated
     LOGICAL :: ARRAYSET = .FALSE.

     ! index of particles walking on microtubule. Default is 0 -> none walking
     INTEGER :: NWALK
     INTEGER, POINTER :: WALKINDS(:)
     
     ! index of particles tethered to a microtubule
     INTEGER :: NTETH
     INTEGER, POINTER :: TETHINDS(:)

     ! walk velocity
     DOUBLE PRECISION :: VEL = 0D0
      
     ! SPRING POSITION
     DOUBLE PRECISION, POINTER :: SPRPOS(:,:)
   
     ! TRACK POSITION
     INTEGER :: NTRACK
     DOUBLE PRECISION, POINTER :: TRKPOS(:,:)
     ! which track is each particle on
     INTEGER, pointer :: PARTTRACK(:)

     
     ! reactive pairs
     DOUBLE PRECISION, POINTER :: RXNPAIR(:,:)
     INTEGER :: NRXN = 0

     ! check for track proximity
     INTEGER :: CHECKTRACKPART ! which particle to check for track proximity
     DOUBLE PRECISION :: TRACKDIST ! contact distance for track
     
  END TYPE PARTICLEGROUP

CONTAINS

  SUBROUTINE CHECKTRACKPROX(PARTP,WHICHTRACK)
    ! check for proximity to a track
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    INTEGER, INTENT(OUT) :: WHICHTRACK
    INTEGER :: PC, TC
    DOUBLE PRECISION :: DIST2
    
    WHICHTRACK = 0
    IF (PARTP%CHECKTRACKPART.LE.0) RETURN
    
    PC = PARTP%CHECKTRACKPART
    DO TC = 1,PARTP%NTRACK
       DIST2 = SUM((PARTP%POS(1:2,PC)-PARTP%TRKPOS(1:2,TC))**2)
       IF (DIST2.LT.PARTP%TRACKDIST**2) THEN
          WHICHTRACK = TC
          RETURN
       ENDIF
    ENDDO
  END SUBROUTINE CHECKTRACKPROX
  
  SUBROUTINE CHECKRXN(PARTP, WHICHREACT)
    ! check if any reactive pair is touching each other
    ! returns only the first rxn found
    ! returns 0 if no rxn found
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    INTEGER, INTENT(OUT) :: WHICHREACT
    INTEGER :: RC, PC1, PC2
    DOUBLE PRECISION :: RVEC(3), DIST2, RD2
    

    WHICHREACT = 0
    
    DO RC = 1,PARTP%NRXN ! cycle through reactive pairs
       PC1 = PARTP%RXNPAIR(RC,1)
       PC2 = PARTP%RXNPAIR(RC,2)
       IF (PC1.GT.PARTP%NP.OR.PC2.GT.PARTP%NP) THEN
          PRINT*, 'ERROR IN CHECKRXN. BAD PARTICLE NUMBER', PC1, PC2, PARTP%NP
          STOP 1
       ENDIF
       RVEC = PARTP%POS(:,PC1)-PARTP%POS(:,PC2)
       DIST2 = SUM(RVEC**2)

       RD2 = (PARTP%CRAD(PC1)+PARTP%CRAD(PC2))**2
       IF (DIST2.LT.RD2) THEN
          WHICHREACT = RC
          RETURN
       ENDIF
    ENDDO
  END SUBROUTINE CHECKRXN
  
  SUBROUTINE GETFORCESTORQUES(PARTP,FORCES,TORQUES)
    ! get forces and torques associated with a given configuration of particles
!    USE KEYS, ONLY: DOMRAD,DOMLEN,KSP
		USE ROTATION, ONLY: QUAT2ROTMAT
		USE KEYS, ONLY: KSP,ATTLEN,DOMLEN,DOMRAD, KCONF, KSTER, KTETH, PERBOUND, STERRAD
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION, INTENT(OUT) :: FORCES(3,PARTP%NP), TORQUES(3,PARTP%NP)
    DOUBLE PRECISION :: E0,EN, DEL, R0
    DOUBLE PRECISION :: POS0(3), FSPRING(3), TQSPRING(3), ROTMAT(3,3), DR(3), UPOS(3)
    DOUBLE PRECISION :: FCONF(3), FSTER(3)
    DOUBLE PRECISION :: RVEC(3), XVEC(3), ATTPT(3), ZVEC(3)
    DOUBLE PRECISION :: DIST
    INTEGER :: PC, I, CC 
    
    FORCES = 0D0
    TORQUES = 0D0
    DEL = 1D-5

    DO PC = 1,PARTP%NP

       FCONF = 0D0
       FSTER = 0D0
       FSPRING = 0D0
       ATTLEN = STERRAD(PC)

       !save initial position
       POS0 = PARTP%POS(:,PC)
       !			print*, 'testxa, POS0: ',POS0
       R0 = SQRT(DOT_PRODUCT(POS0,POS0))
       UPOS = POS0/R0


       !get energy at initial position
       !			CALL GETENERGY(PARTP,E0)
       !			
       !			!find gradient of energy
       !			DO I = 1,3
       !				PARTP%POS(I,PC) = PARTP%POS(I,PC) + DEL
       !				!CALL REFLECTPART(POS0,PARTP%POS(:,PC))				
       !				CALL GETENERGY(PARTP,EN)
       !				FORCES(I,PC) = (E0-EN)/DEL
       !				PARTP%POS(I,PC) = POS0(I)
       !			END DO


       !SPRING FORCE AND TORQUE ON WALKING/TETHERED PARTICLES
       IF(ANY(PARTP%WALKINDS.EQ.PC).OR.ANY(PARTP%TETHINDS.EQ.PC)) THEN
          CALL QUAT2ROTMAT(PARTP%QUAT(:,PC),ROTMAT)
          XVEC = ROTMAT(:,1)
          RVEC = -ATTLEN*XVEC
          ATTPT = PARTP%POS(:,PC)+RVEC
          !				PRINT'(A27,I1, 3ES18.4E2)', 'TESTXA, PC, PART POSITION: ', PC, PARTP%POS(:,PC)
          DR = ATTPT-PARTP%SPRPOS(:,PC)
          IF(ANY(PARTP%WALKINDS.EQ.PC)) THEN
             FSPRING = FSPRING-KSP*DR
          ELSEIF(ANY(PARTP%TETHINDS.EQ.PC)) THEN
             FSPRING = FSPRING-KTETH*DR
          END IF
          TQSPRING = (/RVEC(2)*FSPRING(3)-RVEC(3)*FSPRING(2), &
               RVEC(3)*FSPRING(1)-RVEC(1)*FSPRING(3), &
               RVEC(1)*FSPRING(2)-RVEC(2)*FSPRING(1)/)

       END IF

       !FORCE DUE TO CONFINEMENT
       IF(ABS(DOMLEN).GT.1E-5) THEN
          !CYLINDRICAL CONFINEMENT
          !Z POSITION
          ZVEC = (/0D0,0D0,POS0(3)/)
          IF(POS0(3)+PARTP%SRAD(PC).GT.DOMLEN) THEN
             FCONF = FCONF - KCONF*(PARTP%SRAD(PC)+POS0(3)-DOMLEN)*ZVEC/POS0(3)
          END IF
          !R POSITION
          RVEC = POS0-ZVEC
          DIST = SQRT(DOT_PRODUCT(RVEC,RVEC))
          IF(DIST+PARTP%SRAD(PC).GT.DOMRAD) THEN
!             IF (PC.EQ.PARTP%NP.AND.DIST+PARTP%SRAD(PC).GT.1.2*domrad) PRINT*, 'TESTXA:', PC, DIST+PARTP%SRAD(PC)
             FCONF = FCONF - KCONF*(DIST+PARTP%SRAD(PC)-DOMRAD)*RVEC/DIST
          END IF
       ELSE
          !SPHERICAL CONFINEMENT
          IF((R0+PARTP%SRAD(PC)).GT.DOMRAD) THEN
             FCONF = FCONF - KCONF*(R0-DOMRAD)*UPOS
          END IF
       END IF

       !STERIC FORCE
       DO CC = 1,PARTP%NP
          IF (CC.EQ.PC) CYCLE
          RVEC = PARTP%POS(:,PC)-PARTP%POS(:,CC)
          DIST = SQRT(DOT_PRODUCT(RVEC,RVEC))			
          IF((DIST.LT.(PARTP%SRAD(CC)+PARTP%SRAD(PC))).AND.(DIST.NE.0)) THEN
             !					FSTER = FSTER + KSTER/DIST**13*RVEC/DIST
             FSTER = FSTER + KSTER*(DIST-(PARTP%SRAD(CC)+PARTP%SRAD(PC)))*RVEC/DIST
          END IF
       END DO

       !NO Z-CONFINEMENT IF PERIODIC BC
       IF(PERBOUND) THEN
          FCONF(3) = 0D0
       END IF

       FORCES(:,PC) = FSPRING + FSTER + FCONF
       TORQUES(:,PC) = TQSPRING			

    END DO

  END SUBROUTINE GETFORCESTORQUES
  
  SUBROUTINE GETENERGY(PARTP,ENERGY)
  	!GET ENERGY OF A SYSTEM OF PARTICLES IN A GIVEN CONFIGURATION
  	USE KEYS, ONLY: KSP, KSTER
  	USE ROTATION, ONLY: QUAT2ROTMAT
  	IMPLICIT NONE
		TYPE(PARTICLEGROUP), POINTER :: PARTP
		DOUBLE PRECISION, INTENT(OUT) :: ENERGY
		DOUBLE PRECISION :: DIST, DX
		DOUBLE PRECISION :: ROTMAT(3,3), YVEC(3), ATTPT(3)
		INTEGER :: PC, DC
		
		!ONLY STERIC ENERGY IMPLEMENTED, ADD ENERGY FROM CHAIN ETC
		ENERGY = 0D0

!		!SPRING ENERGY
!		IF(PARTP%WALKIND.NE.0) THEN
!			CALL QUAT2ROTMAT(PARTP%QUAT(:,PARTP%WALKIND),ROTMAT)
!			YVEC = ROTMAT(2,:)
!			ATTPT = PARTP%POS(:,PARTP%WALKIND)-PARTP%SRAD(PARTP%WALKIND)*YVEC
!			DX = SQRT(SUM((ATTPT-PARTP%SPRPOS)**2))
!			ENERGY = ENERGY + 0.5D0*KSP*DX**2
!		END IF		

		!STERIC ENERGY
		DO PC = 1,PARTP%NP-1
			DO DC = PC+1,PARTP%NP
				DIST = SQRT(SUM((PARTP%POS(:,PC)-PARTP%POS(:,DC))**2))
				IF((DIST.LT.(PARTP%SRAD(PC)+PARTP%SRAD(DC))).AND.(DIST.NE.0)) THEN
					!ENERGY OF CONTACT
					ENERGY = ENERGY+KSTER/DIST**12
!					PRINT*, 'TESTX, STERIC ENERGY:', KSTER/DIST**12
				END IF
			END DO
	 END DO
   
  END SUBROUTINE GETENERGY
  
  SUBROUTINE REFLECTPART(POS0,POS1)
  !SUBROUTINE TO REFLECT PARTICLES OFF THE BOUNDARY
  USE KEYS, ONLY:DOMLEN, DOMRAD, PERBOUND
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: POS0(3)
  DOUBLE PRECISION, INTENT(OUT) :: POS1(3)
  DOUBLE PRECISION :: X1,Y1,X2,Y2,X0,Y0,A,B,C,C1,M,DISC
  	
  	IF(PERBOUND) THEN
		  DO WHILE ((ABS(POS1(3)).GT.DOMLEN))
		    POS1(3) = POS1(3) - SIGN(1D0,POS1(3))*2*DOMLEN
  	  END DO
	  ELSE
		  DO WHILE ((ABS(POS1(3)).GT.DOMLEN))
		    POS1(3) = SIGN(1D0,POS1(3))*(2*DOMLEN-ABS(POS1(3)))
  	  END DO
	  END IF
        
    DO WHILE (SQRT((POS1(1)**2+POS1(2)**2)).GT.DOMRAD)
			
			!FIND CO-ORDINATES OF THE INTERSECTION POINT
			X1 = POS0(1)
			Y1 = POS0(2)
			X2 = POS1(1)
      Y2 = POS1(2)
      M = (Y2-Y1)/(X2-X1)
      C1 = Y1-M*X1
      	
      A = 1+M**2
      B = 2*M*C1
      C = C1**2-DOMRAD**2
      
      IF(B**2-4*A*C .LT.0) THEN
      	PRINT*,	'TESTXA:', B**2-4*A*C 
      	PRINT*, 'TESTXA:',POS0,POS1
      	PRINT*, 'TESTXA:',SQRT((POS0(1)**2+POS0(2)**2)),SQRT((POS1(1)**2+POS1(2)**2))
!      	STOP 1
    	END IF
      DISC = SQRT(B**2-4*A*C)
      
      	
      IF(((-B-DISC)/(2*A).GT.MIN(X1,X2)).AND.((-B-DISC)/(2*A).LT.MAX(X1,X2))) THEN
      	X0 = (-B-DISC)/(2*A)
    	ELSE
    		X0 = (-B+DISC)/(2*A)
  		END IF
  			
  		Y0 = M*X0+C1
  			
  		!GET REFLECTION OFF THE TANGENT AT INTERSECTION POINT
  		POS1(1) = X2-2*X0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
  		POS1(2) = Y2-2*Y0*(X0*X2+Y0*Y2-DOMRAD**2)/DOMRAD**2
  			
		END DO

  END SUBROUTINE REFLECTPART

  SUBROUTINE OUTPUTSNAPSHOT(PARTP,FILENAME,INFO,APPEND)
    ! Output a snapshot of current particle configuration
    ! ------------------
    ! input parameters:
    ! PARTP: object describing group of particles
    ! FILENAME: output file
    ! INFO: additional informational float (eg: time)
    ! APPEND: append to file or rewrite?
    ! -------------
    ! output format:
    ! 1st line -> number of particles.
    ! Successive lines -> one line per particle
    ! for each particle, lists INDEX, 3 coordinates of position,
    ! followed by orientation in matrix form (x axis, y axis, z axis of particle frame)
    USE ROTATION, ONLY : QUAT2ROTMAT
    USE KEYS, ONLY : DOMRAD, DOMLEN,NTRIALS
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    DOUBLE PRECISION, INTENT(IN) :: INFO(:)
    LOGICAL, INTENT(IN) :: APPEND
    INTEGER :: PC,NINFO
    DOUBLE PRECISION :: ROTMAT(3,3)
    CHARACTER(LEN=2) :: X1
    CHARACTER(LEN=19) :: FMT1

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN',POSITION='APPEND')
    ELSE
       OPEN(UNIT=99,FILE=FILENAME,STATUS='UNKNOWN')
    END IF
		
		NINFO = SIZE(INFO)
		WRITE (X1,'(I2)') NINFO
		FMT1 = '(I5,3I3,'//TRIM(X1)//'ES18.4E2)'
    WRITE(99,FMT1) NTRIALS,PARTP%NP, PARTP%NWALK, PARTP%NTETH, INFO

    DO PC = 1,PARTP%NP
       ! get rotation matrix from quaternion orientation

       CALL QUAT2ROTMAT(PARTP%QUAT(:,PC),ROTMAT)
       IF(ANY(PARTP%WALKINDS.EQ.PC).OR.ANY(PARTP%TETHINDS.EQ.PC)) THEN
	       WRITE(99,'(I3,15ES18.4E2)') PC, PARTP%POS(:,PC),ROTMAT(:,1),ROTMAT(:,2),ROTMAT(:,3), PARTP%SPRPOS(:,PC)
       ELSE
	       WRITE(99,'(I3,12ES18.4E2)') PC, PARTP%POS(:,PC),ROTMAT(:,1),ROTMAT(:,2),ROTMAT(:,3)
       END IF
    END DO

    CLOSE(99)

  END SUBROUTINE OUTPUTSNAPSHOT

  SUBROUTINE GETINITPOS(PARTP,INITPOS)
    ! get initial positions for particles
    ! walking particles start at the center of the domain
    ! diffusing particles start uniformly at the central cross section
    USE KEYS, ONLY : DOMLEN, DOMRAD, ATTLEN, RANDINITPOS, STERRAD
    USE MT19937, ONLY : GRND
    USE GENUTIL, ONLY : PI
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION, INTENT(OUT) :: INITPOS(3,PARTP%NP)
    DOUBLE PRECISION :: R,TH
    INTEGER :: PC

    DO PC = 1,PARTP%NP
    		ATTLEN = STERRAD(PC)
        IF(ANY(PARTP%WALKINDS.EQ.PC)) THEN
          INITPOS(:,PC) = (/ATTLEN,0D0,0D0/)
        ELSE IF(RANDINITPOS) THEN
          R = (DOMRAD-STERRAD(PC))*SQRT(GRND())
          TH = 2*PI*GRND()
          INITPOS(:,PC) = (/R*COS(TH),R*SIN(TH),0D0/)
        ELSE
        	INITPOS(:,PC) = (/0D0,0D0,0D0/)
        END IF

    END DO

  END SUBROUTINE

  SUBROUTINE INITIALIZEPARTICLES(PARTP,INITPOS,INITQUAT,INITSPRPOS,RANDORIENT)
    ! initialize particle positions and orientations
    ! must provide initial positions
    ! optional: provide initial quaternion
    ! default: canonical orientation, unless RANDORIENT is provided and true
    ! in which case, pick a random orientation
    USE ROTATION, ONLY : UNIFORMRANDQUAT
    USE KEYS, ONLY: ATTLEN, STERRAD
    IMPLICIT NONE

    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION, INTENT(IN) :: INITPOS(3,PARTP%NP)
    DOUBLE PRECISION, INTENT(IN),OPTIONAL :: INITQUAT(4,PARTP%NP), INITSPRPOS(3,PARTP%NP)
    LOGICAL, INTENT(IN), OPTIONAL :: RANDORIENT
!    DOUBLE PRECISION, TARGET :: TMP(3)
    LOGICAL :: DORANDORIENT
    INTEGER :: PC
    DOUBLE PRECISION :: TH

    ! set particle positions
    PARTP%POS = INITPOS
    PARTP%SPRPOS = 0D0

    IF (PRESENT(RANDORIENT)) THEN
       DORANDORIENT = RANDORIENT
    ELSE
       DORANDORIENT = .FALSE.
    ENDIF

    IF (DORANDORIENT) THEN
       ! uniformly random orientation
       DO PC = 1,PARTP%NP
          CALL UNIFORMRANDQUAT(PARTP%QUAT(:,PC))
       ENDDO
    ELSE
       IF (PRESENT(INITQUAT)) THEN
          ! input orientation
          PARTP%QUAT = INITQUAT
       ELSE
          ! cannonical orientation (zero rotation)
          DO PC = 1,PARTP%NP
             PARTP%QUAT(:,PC) = (/1D0,0D0,0D0,0D0/)
          ENDDO
       END IF
    ENDIF
    
    IF (PRESENT(INITSPRPOS)) THEN
    	PARTP%SPRPOS = INITSPRPOS
		ELSE
		  IF (PARTP%NWALK.GT.0.OR.PARTP%NTETH.GT.0) THEN
		  	DO PC = 1,PARTP%NP
		  		IF(ANY(PARTP%WALKINDS.EQ.PC).OR.ANY(PARTP%TETHINDS.EQ.PC)) THEN
	!		    	PARTP%SPRPOS(:,PC) = INITPOS(:,PC)-(/ATTLEN,0D0,0D0/)
	!	    	ELSE IF(ANY(PARTP%TETHINDS.EQ.PC)) THEN 
						ATTLEN = STERRAD(PC)
			  		TH = ATAN2(INITPOS(2,PC),INITPOS(1,PC))
			  		PARTP%SPRPOS(:,PC) = INITPOS(:,PC)-ATTLEN*(/COS(TH),SIN(TH),0D0/)
			  	END IF
				END DO
		  END IF
	  END IF

  END SUBROUTINE INITIALIZEPARTICLES

  SUBROUTINE SETPARTPARAMS(PARTP)
    ! set particle parameters based on keyword arguments
    USE GENUTIL, ONLY : PI
    USE KEYS, ONLY : STERRAD, CONTACTRAD, FRIC, FRICR,SETFRIC,&
                    & SETFRICR, NTETH, TETHINDS, VISCOSITY, NWALK, NTRACK, WALKINDS, VEL,&
                    & DOMRAD, NRXN, RXNPAIR, CHECKTRACKPART,TRACKDIST, VELBIDIR
    USE MT19937, ONLY : GRND
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    DOUBLE PRECISION :: TH,R0, U
    INTEGER :: PC

    IF (.NOT.PARTP%ARRAYSET) THEN
       PRINT*, 'ERROR IN SETPARTPARAMS: particle arrays not allocated'
       STOP 1
    ENDIF

    PARTP%SRAD = STERRAD(1:PARTP%NP)
    PARTP%CRAD = CONTACTRAD(1:PARTP%NP)
    
    
    IF (SETFRIC) THEN
       PARTP%FRIC = FRIC(1:PARTP%NP)
    ELSE
       ! default friction coefficient of sphere
       PARTP%FRIC = 6*PI*VISCOSITY*PARTP%SRAD
    ENDIF

    IF (SETFRICR) THEN
       PARTP%FRICR = FRICR(1:PARTP%NP)
    ELSE
       ! default rotational friction coefficient of sphere
       PARTP%FRICR = 8*PI*VISCOSITY*PARTP%SRAD**3
    ENDIF

    IF(NWALK.GT.0) THEN
       PARTP%WALKINDS = WALKINDS(1:NWALK)
       IF (VELBIDIR) THEN
          U = GRND()
          IF (U.LE.0.5D0) THEN
             PARTP%VEL = VEL
          ELSE
             PARTP%VEL = -VEL
          ENDIF          
       ELSE
          PARTP%VEL = VEL
       ENDIF
    END IF
    
    IF(NTETH.GT.0) THEN
       PARTP%TETHINDS = TETHINDS(1:NTETH)
    END IF

    IF(NTRACK.GT.0) THEN
       DO PC = 1,NTRACK
          TH = 2*PI*GRND()
          R0 = DOMRAD*SQRT(GRND())
          PARTP%TRKPOS(:,PC) = (/R0*COS(TH),R0*SIN(TH),0D0/)
       END DO
       PARTP%CHECKTRACKPART = CHECKTRACKPART
       PARTP%TRACKDIST = TRACKDIST
    END IF

    ! rxn pairs
    IF (NRXN.GT.0) THEN
       PARTP%RXNPAIR = RXNPAIR(1:NRXN,:)
   ENDIF
  END SUBROUTINE SETPARTPARAMS

  SUBROUTINE SETUPPARTICLEGROUP(PARTP,NP,NWALK,NTETH,NTRACK, NRXN)
    ! Initialize arrays for a group of particles
    IMPLICIT NONE

    ! allocate arrays for a particle group, assuming NP particles
    TYPE(PARTICLEGROUP), POINTER :: PARTP
    INTEGER, INTENT(IN) :: NP
    INTEGER, INTENT(IN), OPTIONAL :: NWALK, NTETH, NTRACK, NRXN 

    PARTP%NP = NP

    IF (PRESENT(NRXN)) THEN
       PARTP%NRXN = NRXN
    ELSE
       PARTP%NRXN = 0
    ENDIF
    
    IF(PRESENT(NWALK)) THEN
    	PARTP%NWALK = NWALK
  	ELSE
  		PARTP%NWALK = 0
		END IF
		
    IF(PRESENT(NTETH)) THEN
    	PARTP%NTETH = NTETH
  	ELSE
  		PARTP%NTETH = 0
		END IF
			
    ALLOCATE(PARTP%POS(3,NP),PARTP%QUAT(4,NP))
    ALLOCATE(PARTP%FRIC(NP),PARTP%FRICR(NP),PARTP%SRAD(NP),PARTP%CRAD(NP))
    ALLOCATE(PARTP%SPRPOS(3,NP))
    ALLOCATE(PARTP%TRKPOS(3,NTRACK))
    ALLOCATE(PARTP%WALKINDS(NWALK),PARTP%TETHINDS(NTETH))

    ALLOCATE(PARTP%RXNPAIR(NRXN,2))
    ALLOCATE(PARTP%PARTTRACK(NP))
    IF (NP.GT.0) PARTP%PARTTRACK = 0

    IF (PRESENT(NTRACK)) THEN
       PARTP%NTRACK = NTRACK
    ELSE
       PARTP%NTRACK = 0
    ENDIF
    
    PARTP%ARRAYSET= .TRUE.
  END SUBROUTINE SETUPPARTICLEGROUP

  SUBROUTINE CLEANUPPARTICLEGROUP(PARTP)
    ! deallocate arrays
    IMPLICIT NONE
    TYPE(PARTICLEGROUP), POINTER :: PARTP

    DEALLOCATE(PARTP%POS,PARTP%QUAT,PARTP%FRIC,PARTP%SRAD,PARTP%CRAD,PARTP%FRICR, PARTP%SPRPOS)
    DEALLOCATE(PARTP%WALKINDS,PARTP%TETHINDS, PARTP%RXNPAIR, PARTP%PARTTRACK)
    PARTP%NWALK = 0
    PARTP%NTETH = 0
    PARTP%VEL = 0
    PARTP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPPARTICLEGROUP
END MODULE PARTICLEUTIL

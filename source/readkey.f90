SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3, PC1, PC2, RC
  INTEGER :: REACTWITHALL(NRXNMAX), REACTALLCT
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! input/output
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .false. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  DOPRINT = .TRUE. !print statements to stdout
  NTRIALS = 1

  ! particle mechanics
  NPART = 1
  STERRAD = 0.1
  CONTACTRAD = 0.12
  
  NWALK = 0 !default -> no walking particle
  NTETH = 0 !DEFAULT -> NO TETHERED PARTICLE
  VEL = 0D0
  VELBIDIR = .FALSE.
  KSP = 1D4
  KTETH = 1D4
  KSTER = 1D4
  ATTLEN = 1D0

  ! brownian dynamics
  DELT = 1D-3 ! time step
  FRIC = 1D0 ! friction coefficient
  FRICR = 1D0 ! rotational friction coefficient
  VISCOSITY = 1D0 ! fluid viscosity
  KT = 1D0 ! temperature
  BDSTEPS = 1000 ! number of brownian steps to run
  BDPRINTEVERY = 1 ! how often to print output
  DOBROWN = .TRUE. ! include brownian forces

  ! environment properties
  DOMLEN = 100*MAXVAL(STERRAD) ! length of walking domain
  DOMRAD = 100*MAXVAL(STERRAD) ! radius of walking domain
  NTRACK = 1 !number of microtubule tracks
  PERBOUND = .TRUE.
  KCONF = 1D3

  ! reactive pairs
  NRXN = 0
  RXNPAIR = 0
  RXNFILE = '*.rxn.out'
  REACTALLCT = 0

  CHECKTRACKPART = 0
  TRACKDIST = 0.1D0
  STOPHITTRACK = .FALSE.
  
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1
		 IF(DOPRINT) THEN
	     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     END IF
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
        	CALL READA(ACTION, CASESET=1)
      	CASE('ATTLEN')
      		CALL READF(ATTLEN)
        CASE('BDPRINTEVERY')
        	CALL READI(BDPRINTEVERY)
        CASE('BDSTEPS')
           CALL READI(BDSTEPS)
        CASE('CHECKTRACK')
           CALL READI(CHECKTRACKPART)
           CALL READF(TRACKDIST)
           IF (NITEMS.GT.3) CALL READO(STOPHITTRACK)
        CASE('CONTACTRAD')
        	CALL READF(CONTACTRAD(1))
        	CONTACTRAD = CONTACTRAD(1)
        	 DO I = 2,MIN(MAXNPART,NITEMS-1)
        	   CALL READF(CONTACTRAD(I))
        	ENDDO
      	CASE('CYLCONF')
      		DOMLEN = HUGE(1D0)
        CASE('DELT')
        	CALL READF(DELT)
        CASE('DOMLEN')
        	CALL READF(DOMLEN)
        CASE('DOMRAD')
          CALL READF(DOMRAD)
        CASE('FRIC')
          CALL READF(FRIC(1))
          FRIC = FRIC(1)
          DO I = 2,MIN(MAXNPART,NITEMS-1)
          	CALL READF(FRIC(I))
          ENDDO
          SETFRIC = .TRUE.
      	CASE('FRICR')
          CALL READF(FRICR(1))
          FRICR = FRICR(1)
          DO I = 2,MIN(MAXNPART,NITEMS-1)
          	CALL READF(FRICR(I))
          ENDDO
          SETFRICR = .TRUE.
        CASE('KCONF')
        	CALL READF(KCONF)
        CASE('KSP')
        	CALL READF(KSP)
      	CASE('KSTER')
        	CALL READF(KSTER)
        CASE('KT')
          CALL READF(KT)
        CASE('KTETH')
        	CALL READF(KTETH)
        CASE('NOBROWN')
          DOBROWN = .FALSE.
        CASE('NOCONF')
        	DOMRAD = HUGE(1D0)
        	DOMLEN = HUGE(1D0)
      	CASE('NOPRINT')
      		DOPRINT = .FALSE.
        CASE('NPART')
          CALL READI(NPART)
        CASE('NTRACK')
        	CALL READI(NTRACK)
        CASE('NTRIALS')
        	CALL READI(NTRIALS)
        CASE('OUTFILE')
          CALL READA(OUTFILE)
        CASE('RANDINITPOS')
        	RANDINITPOS = .TRUE.
      	CASE('RANDSPRING')
      		RANDSPRING = .TRUE.
        CASE('REFBOUND')
        	PERBOUND = .FALSE.
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('RXNFILE')
           CALL READA(RXNFILE)
        CASE('RXNPAIR')           
           CALL READI(PC1)
           CALL READI(PC2)
           IF (PC2.LT.0) THEN ! will expand to have pc1 react with all particles
              REACTALLCT = REACTALLCT+1
              REACTWITHALL(REACTALLCT) = PC1
           ELSE
              NRXN = NRXN+1
              RXNPAIR(NRXN,1) = PC1
              RXNPAIR(NRXN,2) = PC2
           ENDIF
        CASE('SNAPSHOTFILE')
          CALL READA (SNAPSHOTFILE)
        CASE('SNAPSHOTS')
          DUMPSNAPSHOTS = .TRUE.
          IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
          IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
          IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STERRAD')
          CALL READF(STERRAD(1))
          STERRAD = STERRAD(1)
          DO I = 2,MIN(MAXNPART,NITEMS-1)
             CALL READF(STERRAD(I))
          ENDDO
        CASE('TETHINDS')
        	DO I = 1,NITEMS-1
        		NTETH = NTETH + 1
        		IF(NTETH.GT.NTETHMAX) THEN
        			PRINT*, 'ERROR: exceeded maximum allowed tethers', NTETHMAX
        			STOP 1
      			END IF
      			CALL READI(TETHINDS(NTETH))
    			END DO
        CASE('VEL')
           CALL READF(VEL)
           IF (NITEMS.GT.2) CALL READO(VELBIDIR)
        CASE('VERBOSE')
          CALL READO(VERBOSE)
        CASE('VISCOSITY')
          CALL READF(VISCOSITY)
        CASE('WALKINDS')
        	DO I = 1,NITEMS-1
        		NWALK = NWALK + 1
        		IF(NWALK.GT.NWALKMAX) THEN
        			PRINT*, 'ERROR: exceeded maximum allowed walkers', NWALKMAX
        			STOP 1
      			END IF
      			CALL READI(WALKINDS(NWALK))
    			END DO
        CASE DEFAULT
          print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO


  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------

  IF (NPART.LT.0.OR.NPART.GT.MAXNPART) THEN
     PRINT*, 'ERROR IN NPART VALUE', NPART, MAXNPART
     STOP 1
  ENDIF
  IF (ANY(STERRAD.LT.0).OR.ANY(CONTACTRAD.LT.0)) THEN
     PRINT*, 'ERROR IN STERRAD OR CONTACTRAD VALUE', STERRAD,CONTACTRAD
     STOP 1
  ENDIF

  ! -------------
  ! expand reactions with all particles
  ! -------------
  DO RC = 1,REACTALLCT
     PC1 = REACTWITHALL(RC)
     DO PC2 = 1,NPART
        IF (PC2.NE.PC1) THEN
           NRXN = NRXN + 1
           RXNPAIR(NRXN,1) = PC1
           RXNPAIR(NRXN,2) = PC2
        ENDIF
     ENDDO
  ENDDO
  
  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(RXNFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument
     ! and additionally the millisecond time
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF
	
	IF(DOPRINT) THEN
		print*, 'Initiating Mersenne twister random number generator with seed:', SEED
		CALL SGRND(SEED)

		print*, '------------Parameter values : -------------------'
		print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
		print*, 'Output file: ', TRIM(OUTFILE)
		IF (DUMPSNAPSHOTS) THEN
		   PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
		ENDIF

  print*, 'Number of particles:', NPART
  print*, 'Start in random position?:', RANDINITPOS
		print*, 'STERRAD:', STERRAD(1:NPART)
		print*, 'CONTACTRAD:', CONTACTRAD(1:NPART)
  PRINT*, 'Friction coefficient:', FRIC(1:NPART)
  print*, 'KSTER, KCONF:', KSTER, KCONF
		IF(NWALK.GT.0) THEN
		  PRINT*, 'Number of walking particle(s): ',NWALK
		  PRINT*, 'Velocity: ',VEL, velbidir
		END IF
		IF(NTRIALS.GT.1) THEN
			PRINT*, 'Number of trials: ',NTRIALS
END IF
PRINT*, 'Reactive pairs:'
DO RC = 1,NRXN
   PRINT*, RXNPAIR(RC,:)
ENDDO
		print*, '----------------------------------------------------'
	END IF
	
END SUBROUTINE READKEY

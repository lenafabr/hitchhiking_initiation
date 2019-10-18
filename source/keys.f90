MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE

  ! Particle properties
  INTEGER :: NPART
  INTEGER, PARAMETER :: MAXNPART=100 !maximum number of input values that can be stored
  DOUBLE PRECISION :: STERRAD(MAXNPART), CONTACTRAD(MAXNPART)
  DOUBLE PRECISION :: FRIC(MAXNPART), FRICR(MAXNPART)
  DOUBLE PRECISION :: VEL !velocity in the z direction
  LOGICAL :: VELBIDIR
  ! FRICTION Coefficients explicitly set
  LOGICAL :: SETFRIC=.FALSE.,SETFRICR=.FALSE.
  LOGICAL :: RANDINITPOS=.FALSE., RANDSPRING=.FALSE.
  LOGICAL :: NOCONF = .FALSE.
  
  ! WALKING PARTICLES
  INTEGER, PARAMETER :: NWALKMAX = 1000
  INTEGER :: WALKINDS(NWALKMAX),NWALK !index of walking particles
  
  ! TETHERED PARTICLES
  INTEGER, PARAMETER :: NTETHMAX = 1000
  INTEGER :: TETHINDS(NTETHMAX),NTETH !index of tethered particles

  ! rxn pairs
  INTEGER, PARAMETER :: NRXNMAX = 1000
  INTEGER :: NRXN
  INTEGER :: RXNPAIR(NRXNMAX,2)
  CHARACTER*100 :: RXNFILE
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY
  LOGICAL :: DOPRINT
  INTEGER :: NTRIALS
  
  ! -----------------
  ! brownian dynamics
  ! -----------------
  DOUBLE PRECISION :: VISCOSITY
  INTEGER :: BDSTEPS
  INTEGER :: BDPRINTEVERY
  LOGICAL :: DOBROWN
  DOUBLE PRECISION :: DELT,KT

  ! ------------------------
  ! environment properties
  ! ------------------------
  DOUBLE PRECISION :: DOMLEN !half length of cylindrical domain
  DOUBLE PRECISION :: DOMRAD !radius of cylindrical domain
  DOUBLE PRECISION :: KSP !spring constant for walker
  DOUBLE PRECISION :: KTETH !spring constant for tether
  DOUBLE PRECISION :: KSTER !magnitude of steric repulsion
  DOUBLE PRECISION :: ATTLEN !distance from center at which spring is attached
  DOUBLE PRECISION :: KCONF !confinement strength
	INTEGER :: NTRACK !number of microtubule tracks
 LOGICAL :: PERBOUND !periodic z boundary? (opposed to reflecting)

 INTEGER :: CHECKTRACKPART ! for this particle, save first time to hit track
 DOUBLE PRECISION :: TRACKDIST ! contact distance from particle to track
 LOGICAL :: STOPHITTRACK

END MODULE KEYS

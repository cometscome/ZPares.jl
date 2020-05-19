C
C  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
C
C
C  This version of MUMPS is provided to you free of charge. It is public
C  domain, based on public domain software developed during the Esprit IV
C  European project PARASOL (1996-1999). Since this first public domain
C  version in 1999, research and developments have been supported by the
C  following institutions: CERFACS, CNRS, ENS Lyon, INPT(ENSEEIHT)-IRIT,
C  INRIA, and University of Bordeaux.
C
C  The MUMPS team at the moment of releasing this version includes
C  Patrick Amestoy, Maurice Bremond, Alfredo Buttari, Abdou Guermouche,
C  Guillaume Joslin, Jean-Yves L'Excellent, Francois-Henry Rouet, Bora
C  Ucar and Clement Weisbecker.
C
C  We are also grateful to Emmanuel Agullo, Caroline Bousquet, Indranil
C  Chowdhury, Philippe Combes, Christophe Daniel, Iain Duff, Vincent Espirat,
C  Aurelia Fevre, Jacko Koster, Stephane Pralet, Chiara Puglisi, Gregoire
C  Richard, Tzvetomila Slavova, Miroslav Tuma and Christophe Voemel who
C  have been contributing to this project.
C
C  Up-to-date copies of the MUMPS package can be obtained
C  from the Web pages:
C  http://mumps.enseeiht.fr/  or  http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED. ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1] and [2]) the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C   [2] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
C   S. Pralet, Hybrid scheduling for the parallel solution of linear
C   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
C
C******************************************************************
C
C  This file contains dummy MPI/BLACS/ScaLAPACK libraries to allow
C  linking/running MUMPS on a platform where MPI is not installed.
C
C******************************************************************
C
C MPI
C
C******************************************************************
C       SUBROUTINE MPI_BSEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM,
C      &            IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR
C       INTEGER BUF(*)
C       WRITE(*,*) 'Error. MPI_BSEND should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_BSEND
C C***********************************************************************
C       SUBROUTINE MPI_BUFFER_ATTACH(BUF, COUNT,  IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, IERR
C       INTEGER BUF(*)
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_BUFFER_ATTACH
C C***********************************************************************
C       SUBROUTINE MPI_BUFFER_DETACH(BUF, COUNT,  IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, IERR
C       INTEGER BUF(*)
C            IERR = 0
C       RETURN
C       END SUBROUTINE MPI_BUFFER_DETACH
C       SUBROUTINE MPI_GATHER( SENDBUF, COUNT, 
C      &         DATATYPE, RECVBUF, RECCOUNT, RECTYPE,
C      &         ROOT, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER COUNT, DATATYPE, RECCOUNT, RECTYPE, ROOT, COMM, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       IF ( RECCOUNT .NE. COUNT ) THEN
C         WRITE(*,*) 'ERROR in MPI_GATHER, RECCOUNT != COUNT'
C         STOP
C       ELSE
C         CALL MUMPS_COPY( COUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C         IF ( IERR .NE. 0 ) THEN
C           WRITE(*,*) 'ERROR in MPI_GATHER, DATATYPE=',DATATYPE
C           STOP
C         END IF
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GATHER
C C***********************************************************************
C       SUBROUTINE MPI_GATHERV( SENDBUF, COUNT, 
C      &         DATATYPE, RECVBUF, RECCOUNT, DISPLS, RECTYPE,
C      &         ROOT, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER COUNT, DATATYPE, RECTYPE, ROOT, COMM, IERR
C       INTEGER RECCOUNT(1)
C       INTEGER SENDBUF(*), RECVBUF(*)
C       INTEGER DISPLS(*)
C C
C C     Note that DISPLS is ignored in this version. One may
C C     want to copy in reception buffer with a shift DISPLS(1).
C C     This requires passing the offset DISPLS(1) to
C C     "MUMPS_COPY_DATATYPE" routines.
C C
C       IF ( RECCOUNT(1) .NE. COUNT ) THEN
C         WRITE(*,*) 'ERROR in MPI_GATHERV, RECCOUNT(1) != COUNT'
C         STOP
C       ELSE
C         CALL MUMPS_COPY( COUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C         IF ( IERR .NE. 0 ) THEN
C           WRITE(*,*) 'ERROR in MPI_GATHERV, DATATYPE=',DATATYPE
C           STOP
C         END IF
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GATHERV
C C***********************************************************************
C       SUBROUTINE MPI_ALLREDUCE( SENDBUF, RECVBUF, COUNT, DATATYPE,
C      &                          OPERATION, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER COUNT, DATATYPE, OPERATION, COMM, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       CALL MUMPS_COPY( COUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C       IF ( IERR .NE. 0 ) THEN
C         WRITE(*,*) 'ERROR in MPI_ALLREDUCE, DATATYPE=',DATATYPE
C         STOP
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_ALLREDUCE
C C***********************************************************************
C       SUBROUTINE MPI_REDUCE( SENDBUF, RECVBUF, COUNT, DATATYPE, OP,
C      &           ROOT, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER COUNT, DATATYPE, OP, ROOT, COMM, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       CALL MUMPS_COPY( COUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C       IF ( IERR .NE. 0 ) THEN
C         WRITE(*,*) 'ERROR in MPI_REDUCE, DATATYPE=',DATATYPE
C         STOP
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_REDUCE
C C***********************************************************************
C       SUBROUTINE MPI_REDUCE_SCATTER( SENDBUF, RECVBUF, RCVCOUNT, 
C      &           DATATYPE, OP, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER RCVCOUNT, DATATYPE, OP, ROOT, COMM, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       CALL MUMPS_COPY( RCVCOUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C       IF ( IERR .NE. 0 ) THEN
C         WRITE(*,*) 'ERROR in MPI_REDUCE_SCATTER, DATATYPE=',DATATYPE
C         STOP
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_REDUCE_SCATTER
C C***********************************************************************
C       SUBROUTINE MPI_ABORT( COMM, IERRCODE, IERR )
C       IMPLICIT NONE
C       INTEGER COMM, IERRCODE, IERR
C       WRITE(*,*) "** MPI_ABORT called"
C       STOP
C       END SUBROUTINE MPI_ABORT
C C***********************************************************************
C       SUBROUTINE MPI_ALLTOALL( SENDBUF, SENDCNT, SENDTYPE,
C      &                         RECVBUF, RECVCNT, RECVTYPE, COMM, IERR )
C       IMPLICIT NONE
C       INTEGER SENDCNT, SENDTYPE, RECVCNT, RECVTYPE, COMM, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       IF ( RECVCNT .NE. SENDCNT ) THEN
C         WRITE(*,*) 'ERROR in MPI_ALLTOALL, RECVCOUNT != SENDCOUNT'
C         STOP
C       ELSE IF ( RECVTYPE .NE. SENDTYPE ) THEN
C         WRITE(*,*) 'ERROR in MPI_ALLTOALL, RECVTYPE != SENDTYPE'
C         STOP
C       ELSE
C         CALL MUMPS_COPY( SENDCNT, SENDBUF, RECVBUF, SENDTYPE, IERR )
C         IF ( IERR .NE. 0 ) THEN
C           WRITE(*,*) 'ERROR in MPI_ALLTOALL, SENDTYPE=',SENDTYPE
C           STOP
C         END IF
C       END IF
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_ALLTOALL
C C***********************************************************************
C       SUBROUTINE MPI_ATTR_PUT( COMM, KEY, VAL, IERR )
C       IMPLICIT NONE
C       INTEGER COMM, KEY, VAL, IERR
C       RETURN
C       END SUBROUTINE MPI_ATTR_PUT
C C***********************************************************************
C       SUBROUTINE MPI_BARRIER( COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_BARRIER
C C***********************************************************************
C       SUBROUTINE MPI_GET_PROCESSOR_NAME( NAME, RESULTLEN, IERROR)
C       CHARACTER (LEN=*) NAME
C       INTEGER RESULTLEN,IERROR
C       RESULTLEN = 1
C       IERROR = 0
C       NAME = 'X'
C       RETURN
C       END SUBROUTINE MPI_GET_PROCESSOR_NAME
C C***********************************************************************
C       SUBROUTINE MPI_BCAST( BUFFER, COUNT, DATATYPE, ROOT, COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, ROOT, COMM, IERR
C       INTEGER BUFFER( * )
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_BCAST
C C***********************************************************************
C       SUBROUTINE MPI_CANCEL( IREQ, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IREQ, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_CANCEL
C C***********************************************************************
C       SUBROUTINE MPI_COMM_CREATE( COMM, GROUP, COMM2, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, GROUP, COMM2, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_CREATE
C C***********************************************************************
C       SUBROUTINE MPI_COMM_DUP( COMM, COMM2, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, COMM2, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_DUP
C C***********************************************************************
C       SUBROUTINE MPI_COMM_FREE( COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_FREE
C C***********************************************************************
C       SUBROUTINE MPI_COMM_GROUP( COMM, GROUP, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, GROUP, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_GROUP
C C***********************************************************************
C       SUBROUTINE MPI_COMM_RANK( COMM, RANK, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, RANK, IERR
C       RANK = 0
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_RANK
C C***********************************************************************
C       SUBROUTINE MPI_COMM_SIZE( COMM, SIZE, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, SIZE, IERR
C       SIZE = 1
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_SIZE
C C***********************************************************************
C       SUBROUTINE MPI_COMM_SPLIT( COMM, COLOR, KEY, COMM2, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COMM, COLOR, KEY, COMM2, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_COMM_SPLIT
C C***********************************************************************
C c     SUBROUTINE MPI_ERRHANDLER_SET( COMM, ERRHANDLER, IERR )
C c     IMPLICIT NONE
C c     INCLUDE 'mpif.h'
C c     INTEGER COMM, ERRHANDLER, IERR
C c     IERR = 0
C c     RETURN
C c     END SUBROUTINE MPI_ERRHANDLER_SET
C C***********************************************************************
C       SUBROUTINE MPI_FINALIZE( IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_FINALIZE
C C***********************************************************************
C       SUBROUTINE MPI_GET_COUNT( STATUS, DATATYPE, COUNT, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER DATATYPE, COUNT, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       WRITE(*,*) 'Error. MPI_GET_COUNT should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GET_COUNT
C C***********************************************************************
C       SUBROUTINE MPI_GROUP_FREE( GROUP, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER GROUP, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GROUP_FREE
C C***********************************************************************
C       SUBROUTINE MPI_GROUP_RANGE_EXCL( GROUP, N, RANGES, GROUP2, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER GROUP, N, GROUP2, IERR
C       INTEGER RANGES(*)
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GROUP_RANGE_EXCL
C C***********************************************************************
C       SUBROUTINE MPI_GROUP_SIZE( GROUP, SIZE, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER GROUP, SIZE, IERR
C       SIZE = 1 ! Or should it be zero ?
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_GROUP_SIZE
C C***********************************************************************
C       SUBROUTINE MPI_INIT(IERR)
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_INIT
C C***********************************************************************
C       SUBROUTINE MPI_INITIALIZED( FLAG, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       LOGICAL FLAG
C       INTEGER IERR
C       FLAG = .TRUE.
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_INITIALIZED
C C***********************************************************************
C       SUBROUTINE MPI_IPROBE( SOURCE, TAG, COMM, FLAG, STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER SOURCE, TAG, COMM, IERR
C       INTEGER STATUS(MPI_STATUS_SIZE)
C       LOGICAL FLAG
C       FLAG = .FALSE.
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_IPROBE
C C***********************************************************************
C       SUBROUTINE MPI_IRECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM,
C      &           IREQ, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, IREQ, IERR
C       INTEGER BUF(*)
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_IRECV
C C***********************************************************************
C       SUBROUTINE MPI_ISEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM,
C      &           IREQ, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR, IREQ
C       INTEGER BUF(*)
C       WRITE(*,*) 'Error. MPI_ISEND should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_ISEND
C C***********************************************************************
C       SUBROUTINE MPI_TYPE_COMMIT( NEWTYP, IERR_MPI )
C       IMPLICIT NONE
C       INTEGER NEWTYP, IERR_MPI
C       RETURN
C       END SUBROUTINE MPI_TYPE_COMMIT
C C***********************************************************************
C       SUBROUTINE MPI_TYPE_FREE( NEWTYP, IERR_MPI )
C       IMPLICIT NONE
C       INTEGER NEWTYP, IERR_MPI
C       RETURN
C       END SUBROUTINE MPI_TYPE_FREE
C C***********************************************************************
C       SUBROUTINE MPI_TYPE_CONTIGUOUS( LENGTH, DATATYPE, NEWTYPE,
C      &                                IERR_MPI )
C       IMPLICIT NONE
C       INTEGER LENGTH, DATATYPE, NEWTYPE, IERR_MPI
C       RETURN
C       END SUBROUTINE MPI_TYPE_CONTIGUOUS
C C***********************************************************************
C       SUBROUTINE MPI_OP_CREATE( FUNC, COMMUTE, OP, IERR )
C       IMPLICIT NONE
C       EXTERNAL FUNC
C       LOGICAL COMMUTE
C       INTEGER OP, IERR
C       OP = 0
C       RETURN
C       END SUBROUTINE MPI_OP_CREATE
C C***********************************************************************
C       SUBROUTINE MPI_OP_FREE( OP, IERR )
C       IMPLICIT NONE
C       INTEGER OP, IERR
C       RETURN
C       END SUBROUTINE MPI_OP_FREE
C C***********************************************************************
C       SUBROUTINE MPI_PACK( INBUF, INCOUNT, DATATYPE, OUTBUF, OUTCOUNT,
C      &           POSITION, COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER INCOUNT, DATATYPE, OUTCOUNT, POSITION, COMM, IERR
C       INTEGER INBUF(*), OUTBUF(*)
C       WRITE(*,*) 'Error. MPI_PACKED should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_PACK
C C***********************************************************************
C       SUBROUTINE MPI_PACK_SIZE( INCOUNT, DATATYPE, COMM, SIZE, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER INCOUNT, DATATYPE, COMM, SIZE, IERR
C       WRITE(*,*) 'Error. MPI_PACK_SIZE should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_PACK_SIZE
C C***********************************************************************
C       SUBROUTINE MPI_PROBE( SOURCE, TAG, COMM, STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER SOURCE, TAG, COMM, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       WRITE(*,*) 'Error. MPI_PROBE should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_PROBE
C C***********************************************************************
C       SUBROUTINE MPI_RECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM,
C      &           STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, IERR
C       INTEGER BUF(*), STATUS(MPI_STATUS_SIZE)
C       WRITE(*,*) 'Error. MPI_RECV should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_RECV
C C***********************************************************************
C       SUBROUTINE MPI_REQUEST_FREE( IREQ, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IREQ, IERR
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_REQUEST_FREE
C C***********************************************************************
C       SUBROUTINE MPI_SEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR
C       INTEGER BUF(*)
C       WRITE(*,*) 'Error. MPI_SEND should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_SEND
C C***********************************************************************
C       SUBROUTINE MPI_SSEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERR)
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR
C       INTEGER BUF(*)
C       WRITE(*,*) 'Error. MPI_SSEND should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_SSEND
C C***********************************************************************
C       SUBROUTINE MPI_TEST( IREQ, FLAG, STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IREQ, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       LOGICAL FLAG
C       FLAG = .FALSE.
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_TEST
C C***********************************************************************
C       SUBROUTINE MPI_UNPACK( INBUF, INSIZE, POSITION, OUTBUF, OUTCOUNT,
C      &           DATATYPE, COMM, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER INSIZE, POSITION, OUTCOUNT, DATATYPE, COMM, IERR
C       INTEGER INBUF(*), OUTBUF(*)
C       WRITE(*,*) 'Error. MPI_UNPACK should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_UNPACK
C C***********************************************************************
C       SUBROUTINE MPI_WAIT( IREQ, STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER IREQ, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       WRITE(*,*) 'Error. MPI_WAIT should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_WAIT
C C***********************************************************************
C       SUBROUTINE MPI_WAITALL( COUNT, ARRAY_OF_REQUESTS, STATUS, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       INTEGER ARRAY_OF_REQUESTS( COUNT )
C       WRITE(*,*) 'Error. MPI_WAITALL should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_WAITALL
C C***********************************************************************
C       SUBROUTINE MPI_WAITANY( COUNT, ARRAY_OF_REQUESTS, INDEX, STATUS,
C      &           IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, INDEX, IERR
C       INTEGER STATUS( MPI_STATUS_SIZE )
C       INTEGER ARRAY_OF_REQUESTS( COUNT )
C       WRITE(*,*) 'Error. MPI_WAITANY should not be called.'
C       STOP
C       IERR = 0
C       RETURN
C       END SUBROUTINE MPI_WAITANY
C C***********************************************************************
C       DOUBLE PRECISION FUNCTION MPI_WTIME( )
C C     elapsed time
C       DOUBLE PRECISION VAL
C C     write(*,*) 'Entering MPI_WTIME'
C       CALL MUMPS_ELAPSE( VAL )
C       MPI_WTIME = VAL
C C     write(*,*) 'Exiting MPI_WTIME'
C       RETURN
C       END FUNCTION MPI_WTIME


C***********************************************************************
C
C  Utilities to copy data
C
C***********************************************************************

C       SUBROUTINE MUMPS_COPY( COUNT, SENDBUF, RECVBUF, DATATYPE, IERR )
C       IMPLICIT NONE
C       INCLUDE 'mpif.h'
C       INTEGER COUNT, DATATYPE, IERR
C       INTEGER SENDBUF(*), RECVBUF(*)
C       IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
C         CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, COUNT )
C       ELSEIF ( DATATYPE .EQ. MPI_LOGICAL ) THEN
C         CALL MUMPS_COPY_LOGICAL( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
C         CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION .OR.
C      &          DATATYPE .EQ. MPI_REAL8 ) THEN
C         CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
C         CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
C         CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
C         CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
C       ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
C         CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, COUNT )
C       ELSE
C         IERR=1
C         RETURN
C       END IF
C       IERR=0
C       RETURN
C       END SUBROUTINE MUMPS_COPY

C       SUBROUTINE MUMPS_COPY_INTEGER( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       INTEGER S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END SUBROUTINE MUMPS_COPY_INTEGER
C       SUBROUTINE MUMPS_COPY_LOGICAL( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       LOGICAL S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END
C       SUBROUTINE MUMPS_COPY_2INTEGER( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       INTEGER S(N+N),R(N+N)
C       INTEGER I
C       DO I = 1, N+N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END SUBROUTINE MUMPS_COPY_2INTEGER
C       SUBROUTINE MUMPS_COPY_REAL( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       REAL S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END
C       SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       DOUBLE PRECISION S(N+N),R(N+N)
C       INTEGER I
C       DO I = 1, N+N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION
C       SUBROUTINE MUMPS_COPY_DOUBLE_PRECISION( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       DOUBLE PRECISION S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END
C       SUBROUTINE MUMPS_COPY_COMPLEX( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C       COMPLEX S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END SUBROUTINE MUMPS_COPY_COMPLEX
C       SUBROUTINE MUMPS_COPY_DOUBLE_COMPLEX( S, R, N )
C       IMPLICIT NONE
C       INTEGER N
C C     DOUBLE COMPLEX S(N),R(N)
C       COMPLEX(kind=kind(0.0D0)) :: S(N),R(N)
C       INTEGER I
C       DO I = 1, N
C         R(I) = S(I)
C       END DO
C       RETURN
C       END


C***********************************************************************
C
C     BLACS
C
C***********************************************************************
      SUBROUTINE blacs_gridinit( CNTXT, C, NPROW, NPCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL
      CHARACTER C
        WRITE(*,*) 'Error. BLACS_GRIDINIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridinit
C***********************************************************************
      SUBROUTINE blacs_gridinfo( CNTXT, NPROW, NPCOL, MYROW, MYCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL, MYROW, MYCOL
        WRITE(*,*) 'Error. BLACS_GRIDINFO should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridinfo
C***********************************************************************
      SUBROUTINE blacs_gridexit( CNTXT )
      IMPLICIT NONE
      INTEGER CNTXT
        WRITE(*,*) 'Error. BLACS_GRIDEXIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE blacs_gridexit


C***********************************************************************
C
C     ScaLAPACK
C
C***********************************************************************
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC,
     &           ICTXT, LLD, INFO )
      IMPLICIT NONE
      INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
        WRITE(*,*) 'Error. DESCINIT should not be called.'
        STOP
      RETURN
      END SUBROUTINE DESCINIT
C***********************************************************************
      INTEGER FUNCTION numroc( N, NB, IPROC, ISRCPROC, NPROCS ) 
      INTEGER N, NB, IPROC, ISRCPROC, NPROCS
C     Can be called
      IF ( NPROCS .ne. 1 ) THEN
        WRITE(*,*) 'Error. Last parameter from NUMROC should be 1'
        STOP
      ENDIF
      IF ( IPROC .ne. 0 ) THEN
        WRITE(*,*) 'Error. IPROC should be 0 in NUMROC.'
        STOP
      ENDIF
      NUMROC = N
      RETURN
      END FUNCTION numroc
C***********************************************************************
      SUBROUTINE pcpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcpotrf
C***********************************************************************
      SUBROUTINE pcgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgetrf
C***********************************************************************
      SUBROUTINE pctrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pctrtrs
C***********************************************************************
      SUBROUTINE pzpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
C     DOUBLE COMPLEX     A( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * )
        WRITE(*,*) 'Error. PZPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzpotrf
C***********************************************************************
      SUBROUTINE pzgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
C     DOUBLE COMPLEX     A( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * )
        WRITE(*,*) 'Error. PZGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgetrf
C***********************************************************************
      SUBROUTINE pztrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
C     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pztrtrs
C***********************************************************************
      SUBROUTINE pspotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pspotrf
C***********************************************************************
      SUBROUTINE psgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgetrf
C***********************************************************************
      SUBROUTINE pstrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pstrtrs
C***********************************************************************
      SUBROUTINE pdpotrf( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdpotrf
C***********************************************************************
      SUBROUTINE pdgetrf( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgetrf
C***********************************************************************
      SUBROUTINE pdtrtrs( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdtrtrs
C***********************************************************************
      SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW,
     &                    MYCOL, LRINDX, LCINDX, RSRC, CSRC )
      IMPLICIT NONE
      INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL,
     &                   MYROW, NPCOL, NPROW, RSRC
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END SUBROUTINE INFOG2L
C***********************************************************************
      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
        INDXG2P = 0
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END FUNCTION INDXG2P
C***********************************************************************
      SUBROUTINE pcscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      COMPLEX            ALPHA
      COMPLEX            X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PCSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcscal
C***********************************************************************
      SUBROUTINE pzscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
C     DOUBLE COMPLEX     ALPHA
C     DOUBLE COMPLEX     X( * )
      COMPLEX(kind=kind(0.0D0)) :: ALPHA, X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PZSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzscal
C***********************************************************************
      SUBROUTINE pdscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PDSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdscal
C***********************************************************************
      SUBROUTINE psscal(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      REAL               ALPHA
      REAL               X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PSSCAL should not be called.'
        STOP
      RETURN
      END SUBROUTINE psscal
C***********************************************************************
      SUBROUTINE pzdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
C     DOUBLE COMPLEX X(*), Y(*)
      COMPLEX(kind=kind(0.0D0)) :: X(*), Y(*)
      DOUBLE PRECISION DOT
        DOT = 0.0d0
        WRITE(*,*) 'Error. PZDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzdot
C***********************************************************************
      SUBROUTINE pcdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      COMPLEX X(*), Y(*)
      REAL DOT
        DOT = 0.0e0
        WRITE(*,*) 'Error. PCDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcdot
C***********************************************************************
      SUBROUTINE pddot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE PRECISION X(*), Y(*), DOT
        DOT = 0.0d0
        WRITE(*,*) 'Error. PDDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pddot
C***********************************************************************
      SUBROUTINE psdot
     &    ( N, DOT, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )
      IMPLICIT NONE
      INTEGER N, IX, JX, IY, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      REAL X(*), Y(*), DOT
        DOT = 0.0e0
        WRITE(*,*) 'Error. PSDOT should not be called.'
        STOP
      RETURN
      END SUBROUTINE psdot
C***********************************************************************
      SUBROUTINE zgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
C     DOUBLE COMPLEX A(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE zgebs2d
C***********************************************************************
      SUBROUTINE cgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE cgebs2d
C***********************************************************************
      SUBROUTINE sgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE sgebs2d
C***********************************************************************
      SUBROUTINE dgebs2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBS2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE dgebs2d
C***********************************************************************
      SUBROUTINE zgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
C     DOUBLE COMPLEX A(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE zgebr2d
C***********************************************************************
      SUBROUTINE cgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE cgebr2d
C***********************************************************************
      SUBROUTINE sgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE sgebr2d
C***********************************************************************
      SUBROUTINE dgebr2d( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBR2D should not be called.'
        STOP
      RETURN
      END SUBROUTINE dgebr2d
C***********************************************************************
      SUBROUTINE pcgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgetrs
C***********************************************************************
      SUBROUTINE pzgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
c     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgetrs
C***********************************************************************
      SUBROUTINE psgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgetrs
C***********************************************************************
      SUBROUTINE pdgetrs( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     &                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgetrs
C***********************************************************************
      SUBROUTINE pcpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      COMPLEX         A( * ), B( * )
        WRITE(*,*) 'Error. PCPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcpotrs
C***********************************************************************
      SUBROUTINE pzpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
c     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzpotrs
C***********************************************************************
      SUBROUTINE pspotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      REAL            A( * ), B( * )
        WRITE(*,*) 'Error. PSPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pspotrs
C***********************************************************************
      SUBROUTINE pdpotrs( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     &           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      DOUBLE          PRECISION A( * ), B( * )
        WRITE(*,*) 'Error. PDPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdpotrs
C***********************************************************************
      SUBROUTINE pscnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL NORM2
      COMPLEX X( * )
        WRITE(*,*) 'Error. PCNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pscnrm2
C***********************************************************************
      SUBROUTINE pdznrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2
C     DOUBLE COMPLEX X( * )
      COMPLEX(kind=kind(0.0D0)) :: X( * )
        WRITE(*,*) 'Error. PZNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdznrm2
C***********************************************************************
      SUBROUTINE psnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL    NORM2, X( * )
        WRITE(*,*) 'Error. PSNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE psnrm2
C***********************************************************************
      SUBROUTINE pdnrm2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2, X( * )
        WRITE(*,*) 'Error. PDNRM2 should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdnrm2
C***********************************************************************
      REAL FUNCTION pclange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      COMPLEX      A( * ), WORK( * )
      PCLANGE = 0.0e0
        WRITE(*,*) 'Error. PCLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pclange
C***********************************************************************
      DOUBLE PRECISION FUNCTION pzlange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PZLANGE = 0.0d0
        WRITE(*,*) 'Error. PZLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pzlange
C***********************************************************************
      REAL FUNCTION pslange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PSLANGE = 0.0e0
        WRITE(*,*) 'Error. PSLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pslange
C***********************************************************************
      DOUBLE PRECISION FUNCTION pdlange( NORM, M, N, A, IA,  JA,
     &                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      DOUBLE       PRECISION A( * ), WORK( * )
      PDLANGE = 0.0d0
        WRITE(*,*) 'Error. PDLANGE should not be called.'
        STOP
      RETURN
      END FUNCTION pdlange
C***********************************************************************
      SUBROUTINE pcgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      COMPLEX         A( * ), WORK( * )
        WRITE(*,*) 'Error. PCGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgecon
C***********************************************************************
      SUBROUTINE pzgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
C     DOUBLE COMPLEX  A( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A( * ), WORK( * )
        WRITE(*,*) 'Error. PZGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgecon
C***********************************************************************
      SUBROUTINE psgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      REAL            A( * ), WORK( * )
        WRITE(*,*) 'Error. PSGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgecon
C***********************************************************************
      SUBROUTINE pdgecon( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     &           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE          PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      DOUBLE          PRECISION A( * ), WORK( * )
        WRITE(*,*) 'Error. PDGECON should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgecon
C***********************************************************************
      SUBROUTINE pcgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
      COMPLEX    A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcgeqpf
C***********************************************************************
      SUBROUTINE pzgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
C     DOUBLE COMPLEX A( * ), TAU( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzgeqpf
C***********************************************************************
      SUBROUTINE psgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      REAL       A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE psgeqpf
C***********************************************************************
      SUBROUTINE pdgeqpf( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     &           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      DOUBLE PRECISION A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDGEQPF should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdgeqpf
C***********************************************************************
      SUBROUTINE pcaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      COMPLEX A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PCAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcaxpy
C***********************************************************************
      SUBROUTINE pzaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
C     DOUBLE COMPLEX A(*),X(*),Y(*)
      COMPLEX(kind=kind(0.0D0)) :: A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PZAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzaxpy
C***********************************************************************
      SUBROUTINE psaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      REAL A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PSAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE psaxpy
C***********************************************************************
      SUBROUTINE pdaxpy(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     &           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE PRECISION A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PDAXPY should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdaxpy
C***********************************************************************
      SUBROUTINE pctrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      COMPLEX            ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pctrsm 
C***********************************************************************
      SUBROUTINE pztrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
C     DOUBLE COMPLEX     ALPHA
      COMPLEX(kind=kind(0.0D0)) ::     ALPHA
      INTEGER            DESCA( * ), DESCB( * )
C     DOUBLE COMPLEX     A( * ), B( * )
      COMPLEX(kind=kind(0.0D0)) ::     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pztrsm 
C***********************************************************************
      SUBROUTINE pstrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      REAL               ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pstrsm 
C***********************************************************************
      SUBROUTINE pdtrsm ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      DOUBLE PRECISION   ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRSM should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdtrsm 
C***********************************************************************
      SUBROUTINE pcunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      COMPLEX   A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCUNMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pcunmqr
C***********************************************************************
      SUBROUTINE pzunmqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
C     DOUBLE COMPLEX A(  *  ), C( * ), TAU( * ), WORK( * )
      COMPLEX(kind=kind(0.0D0)) :: A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZUNMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pzunmqr
C***********************************************************************
      SUBROUTINE psormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      REAL      A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSORMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE psormqr
C***********************************************************************
      SUBROUTINE pdormqr( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     &                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     &                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER         IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER         DESCA( * ), DESCC( * )
      DOUBLE PRECISION  A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDORMQR should not be called.'
        STOP
      RETURN
      END SUBROUTINE pdormqr
C***********************************************************************
      SUBROUTINE chk1mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                    DESCAPOS0, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
      INTEGER            DESCA( * )
        WRITE(*,*) 'Error. CHK1MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE chk1mat
C***********************************************************************
      SUBROUTINE pchk2mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                     DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB,
     &                     DESCB, DESCBPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA,
     &                   MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0,
     &                   NEXTRA
      INTEGER            DESCA( * ), DESCB( * ), EX( NEXTRA ),
     &                   EXPOS( NEXTRA )
        WRITE(*,*) 'Error. PCHK2MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pchk2mat
C***********************************************************************
      SUBROUTINE pxerbla( CONTXT, SRNAME, INFO )
      IMPLICIT NONE
      INTEGER CONTXT, INFO
      CHARACTER SRNAME
        WRITE(*,*) 'Error. PXERBLA should not be called.'
        STOP
      RETURN
      END SUBROUTINE pxerbla
C***********************************************************************
      SUBROUTINE descset( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     &                    LLD )
      IMPLICIT NONE
      INTEGER            ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. DESCSET should not be called.'
        STOP
      RETURN
      END SUBROUTINE descset


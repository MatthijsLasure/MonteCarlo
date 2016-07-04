! File unit numbers used throughout the code.

MODULE MCconstants

    INTEGER, PARAMETER :: IOout  = 7    ! Was 501, but all references should be removed in the code
	INTEGER, PARAMETER :: IOerr  = 1    ! Was 500, but all references should be removed in the code
	INTEGER, PARAMETER :: IOwork = 8    ! For work on files that are closed again in the same block
                                    	! and before calling any other function that uses files
	INTEGER, PARAMETER :: IOdump = 2

	INTEGER, PARAMETER :: IOstartrange = 10

END MODULE MCconstants

/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    mtimesx_RealTimesReal.c
 * Programmer:  James Tursa
 * Version:     1.41
 * Date:        February 23, 2011
 * Copyright:   (c) 2009, 2010, 2011 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the distribution
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * mtimesx_RealTimesReal.c is a support file for the mtimesx.c mex routine.
 *
 * Change Log:
 * 2009/Sep/27 --> 1.00, Initial Release
 * 2009/Dec/03 --> 1.01, Fixed scalar * sparse for scalar with inf or NaN
 * 2009/Dec/05 --> 1.02, Fixed bug, added line scalarmultiply = 0;
 * 2009/Dec/08 --> 1.10, Added singleton expansion capability
 * 2009/Dec/10 --> 1.11, Slight code simplification for singleton expansion
 * 2010/Feb/23 --> 1.20, Fixed bug for dgemv and sgemv calls
 * 2010/Aug/02 --> 1.30, Added (nD scalar) * (nD array) capability
 *                       Replaced buggy mxRealloc with custom code
 * 2010/Oct/04 --> 1.40, Added OpenMP support for custom code
 *                       Expanded sparse * single and sparse * nD support
 *                       Fixed (nD complex scalar)C * (nD array) bug.
 * 2011/Feb/23 --> 1.41, Fixed typos in _syrk and _syr2k BLAS prototypes.
 *
 ****************************************************************************/

/*---------------------------------------------------------------------------------
 * Complex type for function return
 *--------------------------------------------------------------------------------- */

struct RealKindComplex {RealKind r; RealKind i;};

/*---------------------------------------------------------------------------------
 * Function Prototypes
 *--------------------------------------------------------------------------------- */

int AllRealZero(RealKind *x, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqP1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqP1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqM1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqM1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxM1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxM1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxM1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxM1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP0TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxP0TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxP0TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxPxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind ai, RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxPxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n);
void RealKindEqPxPxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealKindEqPxPxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p);
void RealTimesScalar(RealKind *Cpr, RealKind *Cpi, RealKind *Bpr, RealKind *Bpi, char transb,
                     mwSize m2, mwSize n2, RealKind ar, RealKind ai, mwSize n, mwSize p,int);
void RealTimesScalarX(RealKind *Cpr, RealKind *Cpi, RealKind *Bpr, RealKind *Bpi, char transb,
                     mwSize m2, mwSize n2, RealKind ar, RealKind ai, mwSize n, mwSize p, int);
void xFILLPOS(RealKind *Cpr, mwSignedIndex n);
void xFILLNEG(RealKind *Cpr, mwSignedIndex n);
RealKind xDOT(mwSignedIndex *, RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
void     xGER(mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
void     xSYR(char *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *);
void    xSYR2(char *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
void    xGEMV(char *, mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *);
void    xGEMM(char *, char *, mwSignedIndex *, mwSignedIndex *, mwSignedIndex *,
              RealKind *, RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *,
              RealKind *, RealKind *, mwSignedIndex *);
void    xSYRK(char *, char *, mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *,
              mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *);
void   xSYR2K(char *, char *, mwSignedIndex *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *,
              RealKind *, mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *);
void    xAXPY(mwSignedIndex *, RealKind *, RealKind *, mwSignedIndex *, RealKind *, mwSignedIndex *);
struct RealKindComplex RealKindDotProduct(mwSignedIndex, RealKind *, RealKind *, RealKind,
                                                         RealKind *, RealKind *, RealKind, int);
struct RealKindComplex RealKindDotProductX(mwSignedIndex, RealKind *, RealKind *, RealKind,
                                                          RealKind *, RealKind *, RealKind, int);
void RealKindOuterProduct(mwSignedIndex, mwSignedIndex, RealKind *, RealKind *, char,
                          RealKind *, RealKind *, char, RealKind *, RealKind *, int);
void RealKindOuterProductX(mwSignedIndex, mwSignedIndex, RealKind *, RealKind *, char,
                          RealKind *, RealKind *, char, RealKind *, RealKind *, int);
mxArray *RealScalarTimesReal(mxArray *, char, mwSize, mwSize, mxArray *, char, mwSize, mwSize);

/*-------------------------------------------------------------------------------------
 * Function for multiplying two MATLAB variables
 *------------------------------------------------------------------------------------- */

mxArray *RealTimesReal(mxArray *A, char transa, mxArray *B, char transb)
{
    mwSignedIndex inc = 1;
    RealKind Zero = zero;
    RealKind One = one;
    RealKind Minusone = minusone;
    char uplo = 'L';  /* Arbitrary choice. Pick the symmetric case as lower triangular */

    mwSize m1, n1, m2, n2, Andim, Bndim, Cndim, Ap, Bp, Cp, ip, p, Asize, Bsize, Csize, ndim;
    mwSize Ablock, Bblock;
    mwSize *Adims, *Bdims, *Cdims, *Adimz, *Bdimz, *Cindx;
    register mwSignedIndex i, j;
    mwSignedIndex m, n, k, l, lda, ldb, ldc;
    mxArray *C, *result, *rhs[4];
    RealKind *Apr, *Api, *Bpr, *Bpi, *Cpr, *Cpi, *Apr0, *Api0, *Bpr0, *Bpi0;
    RealKind *apr, *api, *bpr, *bpi, *cpr, *cpi;
    RealKind ai, bi, sr, si, aibi;
    RealKind Apr11, Apr12, Apr13, Apr14, 
             Apr21, Apr22, Apr23, Apr24, 
             Apr31, Apr32, Apr33, Apr34,
             Apr41, Apr42, Apr43, Apr44;
    RealKind Api11, Api12, Api13, Api14, 
             Api21, Api22, Api23, Api24, 
             Api31, Api32, Api33, Api34,
             Api41, Api42, Api43, Api44;
    RealKind Bpr11, Bpr12, Bpr13, Bpr14, 
             Bpr21, Bpr22, Bpr23, Bpr24, 
             Bpr31, Bpr32, Bpr33, Bpr34,
             Bpr41, Bpr42, Bpr43, Bpr44;
    RealKind Bpi11, Bpi12, Bpi13, Bpi14, 
             Bpi21, Bpi22, Bpi23, Bpi24, 
             Bpi31, Bpi32, Bpi33, Bpi34,
             Bpi41, Bpi42, Bpi43, Bpi44;
    char ptransa, ptransb;
    char transstring[2] = "_";
    struct RealKindComplex z;
    int scalarmultiply;
	int scalar_method, dot_method, outer_method;
	int singleton_expansion;
	int destroyA = 0, destroyB = 0;

/*--------------------------------------------------------------------------------
 * Get sizes. Note that in the multi-dimensional case, mxGetN returns the product
 * of all of the dimension sizes 2 through end.
 *-------------------------------------------------------------------------------- */

	debug_message = debug;
	threads_used = 0;
    m1 = mxGetM(A);
    n1 = mxGetN(A);
    m2 = mxGetM(B);
    n2 = mxGetN(B);

/*--------------------------------------------------------------------------------
 * Get pointers to the data areas of the operands.
 *-------------------------------------------------------------------------------- */

    Apr0 = Apr = mxGetData(A);
    Api0 = Api = mxGetImagData(A);
    Bpr0 = Bpr = mxGetData(B);
    Bpi0 = Bpi = mxGetImagData(B);

/*--------------------------------------------------------------------------------
 * Simplify transa & transb if appropriate.
 *-------------------------------------------------------------------------------- */

	if( !Api ) {
		if( transa == 'C' ) transa = 'T';
		if( transa == 'G' ) transa = 'N';
	}
	if( !Bpi ) {
		if( transb == 'C' ) transb = 'T';
		if( transb == 'G' ) transb = 'N';
	}

/*--------------------------------------------------------------------------------
 * Scalar expansion cases (custom sparse array code for these cases only).
 * If there is a inf or NaN in the scalar and the other variable is sparse,
 * then don't do the custom code because the sparse zeros will not remain
 * zero. So in that case just fall through and call mtimesx_sparse.
 * (1 x 1) * (K x N) or (M x K) * (1 x 1)
 *-------------------------------------------------------------------------------- */

    scalarmultiply = 0;
    if( m1 == 1 && n1 == 1 ) {
		if( debug ) {
			mexPrintf("MTIMESX: (1 x 1) * (array)\n");
		}
		scalarmultiply = 1;
        if( mxIsSparse(B) ) {
            if( mxIsInf(*Apr) || mxIsNaN(*Apr) ) {
                scalarmultiply = 0;
            } else if( (Api != NULL) && (mxIsInf(*Api) || mxIsNaN(*Api)) ) {
                scalarmultiply = 0;
            }
        }
    } else if( m2 == 1 && n2 == 1 ) {
		if( debug ) {
			mexPrintf("MTIMESX: (array) * (1 x 1)\n");
		}
        scalarmultiply = 1;
        if( mxIsSparse(A) ) {
            if( mxIsInf(*Bpr) || mxIsNaN(*Bpr) ) {
                scalarmultiply = 0;
            } else if( (Bpi != NULL) && (mxIsInf(*Bpi) || mxIsNaN(*Bpi)) ) {
                scalarmultiply = 0;
            }
        }
    }
    if( scalarmultiply ) {
        return RealScalarTimesReal(A, transa, m1, n1, B, transb, m2, n2);
    }

/*--------------------------------------------------------------------------------
 * Multi-dimensional sparse matrix results are not supported in MATLAB, so we are
 * forced to convert the sparse matrix to a full matrix for these cases. Just hope
 * the memory isn't blown.
 *-------------------------------------------------------------------------------- */

	if( mxIsSparse(A) && mxGetNumberOfDimensions(B) > 2 ) {
		k = mexCallMATLAB(1, rhs, 1, &A, "full");
		A = rhs[0];
	    m1 = mxGetM(A);
		n1 = mxGetN(A);
	    Apr0 = Apr = mxGetData(A);
		Api0 = Api = mxGetImagData(A);
		destroyA = 1;
	}
	if( mxIsSparse(B) && mxGetNumberOfDimensions(A) > 2 ) {
		k = mexCallMATLAB(1, rhs, 1, &B, "full");
		B = rhs[0];
	    m2 = mxGetM(B);
		n2 = mxGetN(B);
	    Bpr0 = Bpr = mxGetData(B);
		Bpi0 = Bpi = mxGetImagData(B);
		destroyB = 1;
	}

/*--------------------------------------------------------------------------------
 * Generic sparse matrix or vector operations are not directly supported.
 * So just call an m-file to do the work. Won't save any time, but at least
 * the function will be supported.
 *-------------------------------------------------------------------------------- */

    if( mxIsSparse(A) || mxIsSparse(B) ) {
		if( debug ) {
			mexPrintf("MTIMESX: Unsupported sparse operation ... calling MATLAB intrinsic function mtimes\n");
		}
        rhs[0] = A;
        transstring[0] = transa;
        rhs[1] = mxCreateString(transstring);
        rhs[2] = B;
        transstring[0] = transb;
        rhs[3] = mxCreateString(transstring);
        mexCallMATLAB(1, &result, 4, rhs, "mtimesx_sparse");
        mxDestroyArray(rhs[3]);
        mxDestroyArray(rhs[1]);
        return result;
    }

/*-------------------------------------------------------------------------------
 * Rename array sizes for convenience. Also makes sure that the integer arguments
 * to the BLAS routines are of type mwSignedIndex.
 *------------------------------------------------------------------------------- */

    Andim = mxGetNumberOfDimensions(A);
    Adims = (mwSize *) mxGetDimensions(A);
    Bndim = mxGetNumberOfDimensions(B);
    Bdims = (mwSize *) mxGetDimensions(B);
    m1 = Adims[0];
    n1 = Adims[1];
    m2 = Bdims[0];
    n2 = Bdims[1];

    if( transa == 'N' || transa == 'G' ) {
        m = m1;
        k = n1;
    } else {
        m = n1;
        k = m1;
    }
    if( transb == 'N' || transb == 'G' ) {
        l = m2;
        n = n2;
    } else {
        l = n2;
        n = m2;
    }
    lda = m1;
    ldb = m2;
    ldc = m;

/*-------------------------------------------------------------------------------
 * Check for conforming sizes. Allow nD scalar multiply.
 *------------------------------------------------------------------------------- */

    if( (k != l) && (m*k != 1) && (l*n != 1) ) {
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
        mexErrMsgTxt("Inner matrix dimensions must agree.");
    }

    ndim = (Andim <= Bndim) ? Andim : Bndim;
    for( Cp=2; Cp<ndim; Cp++ ) {
        if( Adims[Cp] != Bdims[Cp] && Adims[Cp] != 1 && Bdims[Cp] != 1 ) {
			if( destroyA ) mxDestroyArray(A);
			if( destroyB ) mxDestroyArray(B);
            mexErrMsgTxt("Dimensions 3:end must agree or be 1 in ND case.");
        }
    }

/*-------------------------------------------------------------------------------
 * Check for nD scalar multiply.
 *------------------------------------------------------------------------------- */

    if( m*k == 1 ) {
        scalarmultiply = 1;
    } else if( l*n == 1 ) {
        scalarmultiply = 2;
    } else {
        scalarmultiply = 0;
    }

/*-------------------------------------------------------------------------------
 * Construct the dimensions of the result. Also use the p variable to keep track
 * of the total number of individual matrix multiples that are involved. The
 * first two dimensions are simply the result of a single matrix multiply, with
 * accouting for the transa and transb pre-operations. The remaining dimensions
 * are copied from A or B, whichever happens to be non-singleton.
 *------------------------------------------------------------------------------- */

    Cndim = (Andim > Bndim) ? Andim : Bndim;
    Cindx = mxMalloc( Cndim * sizeof(*Cindx) );
    Cdims = mxMalloc( Cndim * sizeof(*Cdims) );
    if( scalarmultiply == 1 ) {
        Cdims[0] = l;
        Cdims[1] = n;
    } else if( scalarmultiply == 2 ) {
        Cdims[0] = m;
        Cdims[1] = k;
    } else {
        Cdims[0] = m;
        Cdims[1] = n;
    }
    Adimz = mxMalloc( Cndim * sizeof(*Adimz) );
    Adimz[0] = Adims[0];
    Adimz[1] = Adims[1];
    Bdimz = mxMalloc( Cndim * sizeof(*Bdimz) );
    Bdimz[0] = Bdims[0];
    Bdimz[1] = Bdims[1];
    p = 1;
    for( Cp=2; Cp<Cndim; Cp++ ) {
        Adimz[Cp] = (Cp < Andim) ? Adims[Cp] : 1;
        Bdimz[Cp] = (Cp < Bndim) ? Bdims[Cp] : 1;
        Cdims[Cp] = (Adimz[Cp] > Bdimz[Cp]) ? Adimz[Cp] : Bdimz[Cp];
        p *= Cdims[Cp];
    }
    for( Cp=0; Cp<Cndim; Cp++ ) {
        Cindx[Cp] = 0;
    }

/*------------------------------------------------------------------------------
 * Create output array
 *------------------------------------------------------------------------------ */

    if( mxGetNumberOfElements(A) == 0 || mxGetNumberOfElements(B) == 0 ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
        mxFree(Cindx);
        mxFree(Cdims);
        mxFree(Adimz);
        mxFree(Bdimz);
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
        return result;
    }
    if( mxIsComplex(A) || mxIsComplex(B) ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxCOMPLEX);
    } else {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
    }
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);

/*----------------------------------------------------------------------------
 * See if we can do a simple reshape to do the nD multiply all at once
 *---------------------------------------------------------------------------- */

    if( Andim == 2 && Bndim > 2 && !scalarmultiply && (k == 1 || n == 1 || transb == 'N' || transb == 'G') ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Reshaping nD multiply as a single multiply\n");
		}
        if( transb == 'T' ) transb = 'N';
        if( transb == 'C' ) transb = 'G';
        m2  = k;
        n2  = n * p;
        n   = n2;
        p   = 1;
        ldb = m2;
    }

    if( Bndim == 2 && Andim > 2 && !scalarmultiply && m == 1 ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Reshaping and reordering nD multiply as a single multiply\n");
		}
		apr = Apr; Apr = Bpr; Bpr = apr;
		api = Api; Api = Bpi; Bpi = api;
		ptransa = transa; transa = transb; transb = ptransa;
		if( transa == 'N' ) {
			transa = 'T';
		} else if( transa == 'T' ) {
			transa = 'N';
		} else if( transa == 'G' ) {
			transa = 'C';
		} else {
			transa = 'G';
		}
		if( transb == 'C' ) {
			transb = 'G';
		} else if( transb == 'T' ) {
			transb = 'N';
		}
		m1  = m2;
		n1  = n2;
        m2  = k;
        n2  = m * p;
		if( transa == 'N' || transa == 'G' ) {
			m = m1;
		} else {
			m = n1;
		}
		l   = k;
        n   = n2;
        p   = 1;
		lda = m1;
        ldb = k;
		ldc = m;
    }

/*------------------------------------------------------------------------------
 * Set up conjugate factors
 *------------------------------------------------------------------------------ */

    ai = ( transa == 'C' || transa == 'G' ) ? -one : one;
    bi = ( transb == 'C' || transb == 'G' ) ? -one : one;
	aibi = - ai * bi;

/*----------------------------------------------------------------------------
 * Individual matrix block sizes
 *---------------------------------------------------------------------------- */

    Asize = m1 * n1;
    Bsize = m2 * n2;
    Csize = Cdims[0] * Cdims[1];

#ifdef _OPENMP

/*----------------------------------------------------------------------------
 * Check to see if singleton expansion is present
 *---------------------------------------------------------------------------- */

	singleton_expansion = 0;
	for( i=2; i<Andim; i++ ) {
		if( Adims[i] == 1 ) {
			singleton_expansion = 1;
			break;
		}
	}
	if( !singleton_expansion ) {
		for( i=2; i<Bndim; i++ ) {
			if( Bdims[i] == 1 ) {
				singleton_expansion = 1;
				break;
			}
		}
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (1x1)*(KxN) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		scalarmultiply== 1 && p >= OMP_SPECIAL_SMALL && !singleton_expansion ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Ablock = Asize;
		Bblock = Bsize;
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind sr_, si_;
			RealKind *Apr_, *Bpr_, *Cpr_, *Api_, *Bpi_, *Cpi_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			Api_ = (Api?Api+offset*Asize:NULL);
			Bpi_ = (Bpi?Bpi+offset*Bsize:NULL);
			Cpi_ = (Cpi?Cpi+offset*Csize:NULL);
			for( ip_=0; ip_<p_; ip_++ ) {
		        sr_ = *Apr_;
				si_ = Api_ ? (transa=='N'||transa=='T'?*Api_:-*Api_) : zero;
				RealTimesScalar(Cpr_, Cpi_, Bpr_, Bpi_, transb, m2, n2, sr_, si_, Bblock, 1, METHOD_LOOPS);
				Apr_ += Asize;
				Bpr_ += Bsize;
				Cpr_ += Csize;
				if( Api_ ) {
					Api_ += Asize;
				}
				if( Bpi_ ) {
					Bpi_ += Bsize;
				}
				if( Cpi_ ) {
					Cpi_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (MxK)*(1x1) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		scalarmultiply== 2 && p >= OMP_SPECIAL_SMALL && !singleton_expansion ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Ablock = Asize;
		Bblock = Bsize;
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind sr_, si_;
			RealKind *Apr_, *Bpr_, *Cpr_, *Api_, *Bpi_, *Cpi_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			Api_ = (Api?Api+offset*Asize:NULL);
			Bpi_ = (Bpi?Bpi+offset*Bsize:NULL);
			Cpi_ = (Cpi?Cpi+offset*Csize:NULL);
			for( ip_=0; ip_<p_; ip_++ ) {
		        sr_ = *Bpr_;
				si_ = Bpi_ ? (transb=='N'||transb=='T'?*Bpi_:-*Bpi_) : zero;
				RealTimesScalar(Cpr_, Cpi_, Apr_, Api_, transa, m1, n1, sr_, si_, Ablock, 1, METHOD_LOOPS);
				Apr_ += Asize;
				Bpr_ += Bsize;
				Cpr_ += Csize;
				if( Api_ ) {
					Api_ += Asize;
				}
				if( Bpi_ ) {
					Bpi_ += Bsize;
				}
				if( Cpi_ ) {
					Cpi_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (2x2)*(2x2) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 2 && k == 2 && n == 2 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[2];
						Cpr_[1] = Apr_[2] * Bpr_[0] + Apr_[3] * Bpr_[2];
						Cpr_[2] = Apr_[0] * Bpr_[1] + Apr_[1] * Bpr_[3];
						Cpr_[3] = Apr_[2] * Bpr_[1] + Apr_[3] * Bpr_[3];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1];
						Cpr_[1] = Apr_[2] * Bpr_[0] + Apr_[3] * Bpr_[1];
						Cpr_[2] = Apr_[0] * Bpr_[2] + Apr_[1] * Bpr_[3];
						Cpr_[3] = Apr_[2] * Bpr_[2] + Apr_[3] * Bpr_[3];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			} else {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[2] * Bpr_[2];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[3] * Bpr_[2];
						Cpr_[2] = Apr_[0] * Bpr_[1] + Apr_[2] * Bpr_[3];
						Cpr_[3] = Apr_[1] * Bpr_[1] + Apr_[3] * Bpr_[3];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[2] * Bpr_[1];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[3] * Bpr_[1];
						Cpr_[2] = Apr_[0] * Bpr_[2] + Apr_[2] * Bpr_[3];
						Cpr_[3] = Apr_[1] * Bpr_[2] + Apr_[3] * Bpr_[3];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (2x2)*(2x1) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 2 && k == 2 && l == 2 && n == 1 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1];
					Cpr_[1] = Apr_[2] * Bpr_[0] + Apr_[3] * Bpr_[1];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[2] * Bpr_[1];
					Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[3] * Bpr_[1];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMp inline (1x2)*(2x2) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 1 && k == 2 && n == 2 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transb == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[2];
					Cpr_[1] = Apr_[0] * Bpr_[1] + Apr_[1] * Bpr_[3];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1];
					Cpr_[1] = Apr_[0] * Bpr_[2] + Apr_[1] * Bpr_[3];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (3x3)*(3x3) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 3 && k == 3 && n == 3 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[3] + Apr_[2] * Bpr_[6];
						Cpr_[1] = Apr_[3] * Bpr_[0] + Apr_[4] * Bpr_[3] + Apr_[5] * Bpr_[6];
						Cpr_[2] = Apr_[6] * Bpr_[0] + Apr_[7] * Bpr_[3] + Apr_[8] * Bpr_[6];
						Cpr_[3] = Apr_[0] * Bpr_[1] + Apr_[1] * Bpr_[4] + Apr_[2] * Bpr_[7];
						Cpr_[4] = Apr_[3] * Bpr_[1] + Apr_[4] * Bpr_[4] + Apr_[5] * Bpr_[7];
						Cpr_[5] = Apr_[6] * Bpr_[1] + Apr_[7] * Bpr_[4] + Apr_[8] * Bpr_[7];
						Cpr_[6] = Apr_[0] * Bpr_[2] + Apr_[1] * Bpr_[5] + Apr_[2] * Bpr_[8];
						Cpr_[7] = Apr_[3] * Bpr_[2] + Apr_[4] * Bpr_[5] + Apr_[5] * Bpr_[8];
						Cpr_[8] = Apr_[6] * Bpr_[2] + Apr_[7] * Bpr_[5] + Apr_[8] * Bpr_[8];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1] + Apr_[2] * Bpr_[2];
						Cpr_[1] = Apr_[3] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[5] * Bpr_[2];
						Cpr_[2] = Apr_[6] * Bpr_[0] + Apr_[7] * Bpr_[1] + Apr_[8] * Bpr_[2];
						Cpr_[3] = Apr_[0] * Bpr_[3] + Apr_[1] * Bpr_[4] + Apr_[2] * Bpr_[5];
						Cpr_[4] = Apr_[3] * Bpr_[3] + Apr_[4] * Bpr_[4] + Apr_[5] * Bpr_[5];
						Cpr_[5] = Apr_[6] * Bpr_[3] + Apr_[7] * Bpr_[4] + Apr_[8] * Bpr_[5];
						Cpr_[6] = Apr_[0] * Bpr_[6] + Apr_[1] * Bpr_[7] + Apr_[2] * Bpr_[8];
						Cpr_[7] = Apr_[3] * Bpr_[6] + Apr_[4] * Bpr_[7] + Apr_[5] * Bpr_[8];
						Cpr_[8] = Apr_[6] * Bpr_[6] + Apr_[7] * Bpr_[7] + Apr_[8] * Bpr_[8];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			} else {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[3] * Bpr_[3] + Apr_[6] * Bpr_[6];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[4] * Bpr_[3] + Apr_[7] * Bpr_[6];
						Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[5] * Bpr_[3] + Apr_[8] * Bpr_[6];
						Cpr_[3] = Apr_[0] * Bpr_[1] + Apr_[3] * Bpr_[4] + Apr_[6] * Bpr_[7];
						Cpr_[4] = Apr_[1] * Bpr_[1] + Apr_[4] * Bpr_[4] + Apr_[7] * Bpr_[7];
						Cpr_[5] = Apr_[2] * Bpr_[1] + Apr_[5] * Bpr_[4] + Apr_[8] * Bpr_[7];
						Cpr_[6] = Apr_[0] * Bpr_[2] + Apr_[3] * Bpr_[5] + Apr_[6] * Bpr_[8];
						Cpr_[7] = Apr_[1] * Bpr_[2] + Apr_[4] * Bpr_[5] + Apr_[7] * Bpr_[8];
						Cpr_[8] = Apr_[2] * Bpr_[2] + Apr_[5] * Bpr_[5] + Apr_[8] * Bpr_[8];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[3] * Bpr_[1] + Apr_[6] * Bpr_[2];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[7] * Bpr_[2];
						Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[5] * Bpr_[1] + Apr_[8] * Bpr_[2];
						Cpr_[3] = Apr_[0] * Bpr_[3] + Apr_[3] * Bpr_[4] + Apr_[6] * Bpr_[5];
						Cpr_[4] = Apr_[1] * Bpr_[3] + Apr_[4] * Bpr_[4] + Apr_[7] * Bpr_[5];
						Cpr_[5] = Apr_[2] * Bpr_[3] + Apr_[5] * Bpr_[4] + Apr_[8] * Bpr_[5];
						Cpr_[6] = Apr_[0] * Bpr_[6] + Apr_[3] * Bpr_[7] + Apr_[6] * Bpr_[8];
						Cpr_[7] = Apr_[1] * Bpr_[6] + Apr_[4] * Bpr_[7] + Apr_[7] * Bpr_[8];
						Cpr_[8] = Apr_[2] * Bpr_[6] + Apr_[5] * Bpr_[7] + Apr_[8] * Bpr_[8];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (3x3)*(3x1) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 3 && k == 3 && l == 3 && n == 1 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1] + Apr_[2] * Bpr_[2];
					Cpr_[1] = Apr_[3] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[5] * Bpr_[2];
					Cpr_[2] = Apr_[6] * Bpr_[0] + Apr_[7] * Bpr_[1] + Apr_[8] * Bpr_[2];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[3] * Bpr_[1] + Apr_[6] * Bpr_[2];
					Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[7] * Bpr_[2];
					Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[5] * Bpr_[1] + Apr_[8] * Bpr_[2];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMp inline (1x3)*(3x3) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 1 && k == 3 && n == 3 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transb == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[3] + Apr_[2] * Bpr_[6];
					Cpr_[1] = Apr_[0] * Bpr_[1] + Apr_[1] * Bpr_[4] + Apr_[2] * Bpr_[7];
					Cpr_[2] = Apr_[0] * Bpr_[2] + Apr_[1] * Bpr_[5] + Apr_[2] * Bpr_[8];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[1] + Apr_[2] * Bpr_[2];
					Cpr_[1] = Apr_[0] * Bpr_[3] + Apr_[1] * Bpr_[4] + Apr_[2] * Bpr_[5];
					Cpr_[2] = Apr_[0] * Bpr_[6] + Apr_[1] * Bpr_[7] + Apr_[2] * Bpr_[8];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (4x4)*(4x4) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 4 && k == 4 && n == 4 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[ 0] * Bpr_[0] + Apr_[ 1] * Bpr_[4] + Apr_[ 2] * Bpr_[ 8] + Apr_[ 3] * Bpr_[12];
						Cpr_[1] = Apr_[ 4] * Bpr_[0] + Apr_[ 5] * Bpr_[4] + Apr_[ 6] * Bpr_[ 8] + Apr_[ 7] * Bpr_[12];
						Cpr_[2] = Apr_[ 8] * Bpr_[0] + Apr_[ 9] * Bpr_[4] + Apr_[10] * Bpr_[ 8] + Apr_[11] * Bpr_[12];
						Cpr_[3] = Apr_[12] * Bpr_[0] + Apr_[13] * Bpr_[4] + Apr_[14] * Bpr_[ 8] + Apr_[15] * Bpr_[12];
						Cpr_[4] = Apr_[ 0] * Bpr_[1] + Apr_[ 1] * Bpr_[5] + Apr_[ 2] * Bpr_[ 9] + Apr_[ 3] * Bpr_[13];
						Cpr_[5] = Apr_[ 4] * Bpr_[1] + Apr_[ 5] * Bpr_[5] + Apr_[ 6] * Bpr_[ 9] + Apr_[ 7] * Bpr_[13];
						Cpr_[6] = Apr_[ 8] * Bpr_[1] + Apr_[ 9] * Bpr_[5] + Apr_[10] * Bpr_[ 9] + Apr_[11] * Bpr_[13];
						Cpr_[7] = Apr_[12] * Bpr_[1] + Apr_[13] * Bpr_[5] + Apr_[14] * Bpr_[ 9] + Apr_[15] * Bpr_[13];
						Cpr_[8] = Apr_[ 0] * Bpr_[2] + Apr_[ 1] * Bpr_[6] + Apr_[ 2] * Bpr_[10] + Apr_[ 3] * Bpr_[14];
						Cpr_[9] = Apr_[ 4] * Bpr_[2] + Apr_[ 5] * Bpr_[6] + Apr_[ 6] * Bpr_[10] + Apr_[ 7] * Bpr_[14];
						Cpr_[10]= Apr_[ 8] * Bpr_[2] + Apr_[ 9] * Bpr_[6] + Apr_[10] * Bpr_[10] + Apr_[11] * Bpr_[14];
						Cpr_[11]= Apr_[12] * Bpr_[2] + Apr_[13] * Bpr_[6] + Apr_[14] * Bpr_[10] + Apr_[15] * Bpr_[14];
						Cpr_[12]= Apr_[ 0] * Bpr_[3] + Apr_[ 1] * Bpr_[7] + Apr_[ 2] * Bpr_[11] + Apr_[ 3] * Bpr_[15];
						Cpr_[13]= Apr_[ 4] * Bpr_[3] + Apr_[ 5] * Bpr_[7] + Apr_[ 6] * Bpr_[11] + Apr_[ 7] * Bpr_[15];
						Cpr_[14]= Apr_[ 8] * Bpr_[3] + Apr_[ 9] * Bpr_[7] + Apr_[10] * Bpr_[11] + Apr_[11] * Bpr_[15];
						Cpr_[15]= Apr_[12] * Bpr_[3] + Apr_[13] * Bpr_[7] + Apr_[14] * Bpr_[11] + Apr_[15] * Bpr_[15];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[ 0] * Bpr_[ 0] + Apr_[ 1] * Bpr_[ 1] + Apr_[ 2] * Bpr_[ 2] + Apr_[ 3] * Bpr_[ 3];
						Cpr_[1] = Apr_[ 4] * Bpr_[ 0] + Apr_[ 5] * Bpr_[ 1] + Apr_[ 6] * Bpr_[ 2] + Apr_[ 7] * Bpr_[ 3];
						Cpr_[2] = Apr_[ 8] * Bpr_[ 0] + Apr_[ 9] * Bpr_[ 1] + Apr_[10] * Bpr_[ 2] + Apr_[11] * Bpr_[ 3];
						Cpr_[3] = Apr_[12] * Bpr_[ 0] + Apr_[13] * Bpr_[ 1] + Apr_[14] * Bpr_[ 2] + Apr_[15] * Bpr_[ 3];
						Cpr_[4] = Apr_[ 0] * Bpr_[ 4] + Apr_[ 1] * Bpr_[ 5] + Apr_[ 2] * Bpr_[ 6] + Apr_[ 3] * Bpr_[ 7];
						Cpr_[5] = Apr_[ 4] * Bpr_[ 4] + Apr_[ 5] * Bpr_[ 5] + Apr_[ 6] * Bpr_[ 6] + Apr_[ 7] * Bpr_[ 7];
						Cpr_[6] = Apr_[ 8] * Bpr_[ 4] + Apr_[ 9] * Bpr_[ 5] + Apr_[10] * Bpr_[ 6] + Apr_[11] * Bpr_[ 7];
						Cpr_[7] = Apr_[12] * Bpr_[ 4] + Apr_[13] * Bpr_[ 5] + Apr_[14] * Bpr_[ 6] + Apr_[15] * Bpr_[ 7];
						Cpr_[8] = Apr_[ 0] * Bpr_[ 8] + Apr_[ 1] * Bpr_[ 9] + Apr_[ 2] * Bpr_[10] + Apr_[ 3] * Bpr_[11];
						Cpr_[9] = Apr_[ 4] * Bpr_[ 8] + Apr_[ 5] * Bpr_[ 9] + Apr_[ 6] * Bpr_[10] + Apr_[ 7] * Bpr_[11];
						Cpr_[10]= Apr_[ 8] * Bpr_[ 8] + Apr_[ 9] * Bpr_[ 9] + Apr_[10] * Bpr_[10] + Apr_[11] * Bpr_[11];
						Cpr_[11]= Apr_[12] * Bpr_[ 8] + Apr_[13] * Bpr_[ 9] + Apr_[14] * Bpr_[10] + Apr_[15] * Bpr_[11];
						Cpr_[12]= Apr_[ 0] * Bpr_[12] + Apr_[ 1] * Bpr_[13] + Apr_[ 2] * Bpr_[14] + Apr_[ 3] * Bpr_[15];
						Cpr_[13]= Apr_[ 4] * Bpr_[12] + Apr_[ 5] * Bpr_[13] + Apr_[ 6] * Bpr_[14] + Apr_[ 7] * Bpr_[15];
						Cpr_[14]= Apr_[ 8] * Bpr_[12] + Apr_[ 9] * Bpr_[13] + Apr_[10] * Bpr_[14] + Apr_[11] * Bpr_[15];
						Cpr_[15]= Apr_[12] * Bpr_[12] + Apr_[13] * Bpr_[13] + Apr_[14] * Bpr_[14] + Apr_[15] * Bpr_[15];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			} else {
				if( transb == 'T' ) {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[4] * Bpr_[4] + Apr_[8] * Bpr_[8] + Apr_[12]* Bpr_[12];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[5] * Bpr_[4] + Apr_[9] * Bpr_[8] + Apr_[13]* Bpr_[12];
						Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[6] * Bpr_[4] + Apr_[10]* Bpr_[8] + Apr_[14]* Bpr_[12];
						Cpr_[3] = Apr_[3] * Bpr_[0] + Apr_[7] * Bpr_[4] + Apr_[11]* Bpr_[8] + Apr_[15]* Bpr_[12];
						Cpr_[4] = Apr_[0] * Bpr_[1] + Apr_[4] * Bpr_[5] + Apr_[8] * Bpr_[9] + Apr_[12]* Bpr_[13];
						Cpr_[5] = Apr_[1] * Bpr_[1] + Apr_[5] * Bpr_[5] + Apr_[9] * Bpr_[9] + Apr_[13]* Bpr_[13];
						Cpr_[6] = Apr_[2] * Bpr_[1] + Apr_[6] * Bpr_[5] + Apr_[10]* Bpr_[9] + Apr_[14]* Bpr_[13];
						Cpr_[7] = Apr_[3] * Bpr_[1] + Apr_[7] * Bpr_[5] + Apr_[11]* Bpr_[9] + Apr_[15]* Bpr_[13];
						Cpr_[8] = Apr_[0] * Bpr_[2] + Apr_[4] * Bpr_[6] + Apr_[8] * Bpr_[10]+ Apr_[12]* Bpr_[14];
						Cpr_[9] = Apr_[1] * Bpr_[2] + Apr_[5] * Bpr_[6] + Apr_[9] * Bpr_[10]+ Apr_[13]* Bpr_[14];
						Cpr_[10]= Apr_[2] * Bpr_[2] + Apr_[6] * Bpr_[6] + Apr_[10]* Bpr_[10]+ Apr_[14]* Bpr_[14];
						Cpr_[11]= Apr_[3] * Bpr_[2] + Apr_[7] * Bpr_[6] + Apr_[11]* Bpr_[10]+ Apr_[15]* Bpr_[14];
						Cpr_[12]= Apr_[0] * Bpr_[3] + Apr_[4] * Bpr_[7] + Apr_[8] * Bpr_[11]+ Apr_[12]* Bpr_[15];
						Cpr_[13]= Apr_[1] * Bpr_[3] + Apr_[5] * Bpr_[7] + Apr_[9] * Bpr_[11]+ Apr_[13]* Bpr_[15];
						Cpr_[14]= Apr_[2] * Bpr_[3] + Apr_[6] * Bpr_[7] + Apr_[10]* Bpr_[11]+ Apr_[14]* Bpr_[15];
						Cpr_[15]= Apr_[3] * Bpr_[3] + Apr_[7] * Bpr_[7] + Apr_[11]* Bpr_[11]+ Apr_[15]* Bpr_[15];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				} else {
					for( ip_=0; ip_<p_; ip_++ ) {
						Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[8] * Bpr_[2] + Apr_[12]* Bpr_[3];
						Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[5] * Bpr_[1] + Apr_[9] * Bpr_[2] + Apr_[13]* Bpr_[3];
						Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[6] * Bpr_[1] + Apr_[10]* Bpr_[2] + Apr_[14]* Bpr_[3];
						Cpr_[3] = Apr_[3] * Bpr_[0] + Apr_[7] * Bpr_[1] + Apr_[11]* Bpr_[2] + Apr_[15]* Bpr_[3];
						Cpr_[4] = Apr_[0] * Bpr_[4] + Apr_[4] * Bpr_[5] + Apr_[8] * Bpr_[6] + Apr_[12]* Bpr_[7];
						Cpr_[5] = Apr_[1] * Bpr_[4] + Apr_[5] * Bpr_[5] + Apr_[9] * Bpr_[6] + Apr_[13]* Bpr_[7];
						Cpr_[6] = Apr_[2] * Bpr_[4] + Apr_[6] * Bpr_[5] + Apr_[10]* Bpr_[6] + Apr_[14]* Bpr_[7];
						Cpr_[7] = Apr_[3] * Bpr_[4] + Apr_[7] * Bpr_[5] + Apr_[11]* Bpr_[6] + Apr_[15]* Bpr_[7];
						Cpr_[8] = Apr_[0] * Bpr_[8] + Apr_[4] * Bpr_[9] + Apr_[8] * Bpr_[10]+ Apr_[12]* Bpr_[11];
						Cpr_[9] = Apr_[1] * Bpr_[8] + Apr_[5] * Bpr_[9] + Apr_[9] * Bpr_[10]+ Apr_[13]* Bpr_[11];
						Cpr_[10]= Apr_[2] * Bpr_[8] + Apr_[6] * Bpr_[9] + Apr_[10]* Bpr_[10]+ Apr_[14]* Bpr_[11];
						Cpr_[11]= Apr_[3] * Bpr_[8] + Apr_[7] * Bpr_[9] + Apr_[11]* Bpr_[10]+ Apr_[15]* Bpr_[11];
						Cpr_[12]= Apr_[0] * Bpr_[12]+ Apr_[4] * Bpr_[13]+ Apr_[8] * Bpr_[14]+ Apr_[12]* Bpr_[15];
						Cpr_[13]= Apr_[1] * Bpr_[12]+ Apr_[5] * Bpr_[13]+ Apr_[9] * Bpr_[14]+ Apr_[13]* Bpr_[15];
						Cpr_[14]= Apr_[2] * Bpr_[12]+ Apr_[6] * Bpr_[13]+ Apr_[10]* Bpr_[14]+ Apr_[14]* Bpr_[15];
						Cpr_[15]= Apr_[3] * Bpr_[12]+ Apr_[7] * Bpr_[13]+ Apr_[11]* Bpr_[14]+ Apr_[15]* Bpr_[15];
						Apr_ += Asize;
						Bpr_ += Bsize;
						Cpr_ += Csize;
					}
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMP inline (4x4)*(4x1) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 4 && k == 4 && l == 4 && n == 1 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transa == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[ 0] * Bpr_[0] + Apr_[ 1] * Bpr_[1] + Apr_[ 2] * Bpr_[2] + Apr_[ 3] * Bpr_[3];
					Cpr_[1] = Apr_[ 4] * Bpr_[0] + Apr_[ 5] * Bpr_[1] + Apr_[ 6] * Bpr_[2] + Apr_[ 7] * Bpr_[3];
					Cpr_[2] = Apr_[ 8] * Bpr_[0] + Apr_[ 9] * Bpr_[1] + Apr_[10] * Bpr_[2] + Apr_[11] * Bpr_[3];
					Cpr_[3] = Apr_[12] * Bpr_[0] + Apr_[13] * Bpr_[1] + Apr_[14] * Bpr_[2] + Apr_[15] * Bpr_[3];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[4] * Bpr_[1] + Apr_[ 8] * Bpr_[2] + Apr_[12] * Bpr_[3];
					Cpr_[1] = Apr_[1] * Bpr_[0] + Apr_[5] * Bpr_[1] + Apr_[ 9] * Bpr_[2] + Apr_[13] * Bpr_[3];
					Cpr_[2] = Apr_[2] * Bpr_[0] + Apr_[6] * Bpr_[1] + Apr_[10] * Bpr_[2] + Apr_[14] * Bpr_[3];
					Cpr_[3] = Apr_[3] * Bpr_[0] + Apr_[7] * Bpr_[1] + Apr_[11] * Bpr_[2] + Apr_[15] * Bpr_[3];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

/*----------------------------------------------------------------------------
 * Check to see if we can do a fast OpenMp inline (1x4)*(4x4) nD multiply
 *---------------------------------------------------------------------------- */

    if( (mtimesx_mode == MTIMESX_LOOPS_OMP || mtimesx_mode == MTIMESX_SPEED_OMP) && max_threads > 1 &&
		m == 1 && k == 4 && n == 4 && p >= OMP_SPECIAL_SMALL && !singleton_expansion && !Cpi ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS (unrolled into inline multiplies)\n");
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}
		Asize *= (Andim > 2);
		Bsize *= (Bndim > 2);
		omp_set_dynamic(1);
#pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Bpr_, *Cpr_;
			mwSize ip_, p_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = p / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				p_ = p - offset;
			} else {
				p_ = blocksize;
			}
			Apr_ = Apr + offset * Asize;
			Bpr_ = Bpr + offset * Bsize;
			Cpr_ = Cpr + offset * Csize;
			if( transb == 'T' ) {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[0] + Apr_[1] * Bpr_[4] + Apr_[2] * Bpr_[ 8] + Apr_[3] * Bpr_[12];
					Cpr_[1] = Apr_[0] * Bpr_[1] + Apr_[1] * Bpr_[5] + Apr_[2] * Bpr_[ 9] + Apr_[3] * Bpr_[13];
					Cpr_[2] = Apr_[0] * Bpr_[2] + Apr_[1] * Bpr_[6] + Apr_[2] * Bpr_[10] + Apr_[3] * Bpr_[14];
					Cpr_[3] = Apr_[0] * Bpr_[3] + Apr_[1] * Bpr_[7] + Apr_[2] * Bpr_[11] + Apr_[3] * Bpr_[15];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			} else {
				for( ip_=0; ip_<p_; ip_++ ) {
					Cpr_[0] = Apr_[0] * Bpr_[ 0] + Apr_[1] * Bpr_[ 1] + Apr_[2] * Bpr_[ 2] + Apr_[3] * Bpr_[ 3];
					Cpr_[1] = Apr_[0] * Bpr_[ 4] + Apr_[1] * Bpr_[ 5] + Apr_[2] * Bpr_[ 6] + Apr_[3] * Bpr_[ 7];
					Cpr_[2] = Apr_[0] * Bpr_[ 8] + Apr_[1] * Bpr_[ 9] + Apr_[2] * Bpr_[10] + Apr_[3] * Bpr_[11];
					Cpr_[3] = Apr_[0] * Bpr_[12] + Apr_[1] * Bpr_[13] + Apr_[2] * Bpr_[14] + Apr_[3] * Bpr_[15];
					Apr_ += Asize;
					Bpr_ += Bsize;
					Cpr_ += Csize;
				}
			}
		}
	    mxFree(Cindx);
	    mxFree(Cdims);
	    mxFree(Adimz);
	    mxFree(Bdimz);
/*		
        if( AllRealZero(Cpi, m*n*p) ) {
            mxFree(Cpi);
            mxSetImagData(C, NULL);
        }
*/
		if( destroyA ) mxDestroyArray(A);
		if( destroyB ) mxDestroyArray(B);
		return result;
	}

#endif

/*----------------------------------------------------------------------------
 * Get dot product method to use.
 *---------------------------------------------------------------------------- */

	if( mtimesx_mode == MTIMESX_BLAS || mtimesx_mode == MTIMESX_MATLAB ) {
		dot_method = METHOD_BLAS;
	} else if( mtimesx_mode == MTIMESX_LOOPS ) {
	    dot_method = METHOD_LOOPS;
	} else if( mtimesx_mode == MTIMESX_LOOPS_OMP ) {
	    dot_method = METHOD_LOOPS_OMP;
	} else {
#if !defined(_MSC_VER) || _MSC_VER < 1500  /* Version 2008, 9.0 */
		if( (Apr != Bpr) && !Api && !Bpi ) {
			dot_method = METHOD_BLAS;
		} else if( mtimesx_mode == MTIMESX_SPEED_OMP && max_threads > 1 ) {
#else
		if( mtimesx_mode == MTIMESX_SPEED_OMP && max_threads > 1 ) {
#endif
		    if( Apr != Bpr ) {
			    if( Api && Bpi ) {
				    dot_method = METHOD_LOOPS_OMP;
			    } else {
				    dot_method = METHOD_LOOPS;
			    }
		    } else {
			    if( Api && (ai * bi == -one) ) {
				    dot_method = METHOD_BLAS;
			    } else {
				    dot_method = METHOD_LOOPS_OMP;
			    }
		    }
		} else {
#ifdef __LCC__
			if( (Apr == Bpr && (!Api || (ai * bi == -one))) || omp_get_num_procs() > 2 ) {
				dot_method = METHOD_BLAS;
			} else {
				dot_method = METHOD_LOOPS;
			}
#else
			if( Apr == Bpr && (!Api || (ai * bi == -one)) ) {
				dot_method = METHOD_BLAS;
			} else {
				dot_method = METHOD_LOOPS;
			}
#endif
		}
	}

/*----------------------------------------------------------------------------
 * Get outer product method to use.
 *---------------------------------------------------------------------------- */

	switch( mtimesx_mode )
	{
	case MTIMESX_BLAS:
		outer_method = METHOD_BLAS;
		break;
	case MTIMESX_LOOPS:
		outer_method = METHOD_LOOPS;
		break;
	case MTIMESX_LOOPS_OMP:
		outer_method = METHOD_LOOPS_OMP;
		break;
	case MTIMESX_SPEED_OMP:
		if( max_threads > 1 ) {
		    outer_method = METHOD_LOOPS_OMP;
			break;
		}
	case MTIMESX_MATLAB:
	case MTIMESX_SPEED:
#ifdef __LCC__
		if( Api && Bpi && omp_get_num_procs() <= 2 ) {
#else
		if( (Apr == Bpr) || (Api && Bpi) ) {
#endif
			outer_method = METHOD_LOOPS;
		} else {
			outer_method = METHOD_BLAS;
		}
		break;
	}

/*----------------------------------------------------------------------------
 * Get scalar product method to use.
 *---------------------------------------------------------------------------- */

	switch( mtimesx_mode )
	{
	case MTIMESX_BLAS:
		scalar_method = METHOD_BLAS;
		break;
	case MTIMESX_LOOPS:
		scalar_method = METHOD_LOOPS;
		break;
	case MTIMESX_LOOPS_OMP:
		scalar_method = METHOD_LOOPS_OMP;
		break;
	case MTIMESX_SPEED_OMP:
		if( max_threads > 1 ) {
		    scalar_method = METHOD_LOOPS_OMP;
			break;
		}
	case MTIMESX_MATLAB:
	case MTIMESX_SPEED:
		if( ai != zero && Bpi ) {
			scalar_method = METHOD_LOOPS;
		} else {
			scalar_method = METHOD_BLAS;
		}
		break;
	}

/*----------------------------------------------------------------------------
 * Outer Loop to process all of the individual matrix multiplies
 *---------------------------------------------------------------------------- */

	if( debug ) {
		mexPrintf("MTIMESX: Performing %d individual multiplies\n",p);
	}

    for( ip=0; ip<p; ip++ ) {
        ptransa = transa;  /* Restore the original transa and transb, because */
        ptransb = transb;  /* they might have been changed in previous iteration */
		if( debug_message ) {
			mexPrintf("MTIMESX: (%d x %d) * (%d x %d)\n",m,k,l,n);
		}

/*----------------------------------------------------------------------------
 * Scalar product (1 x 1) * (K x N)
 *---------------------------------------------------------------------------- */

    if( scalarmultiply == 1 ) {
        sr = *Apr;
		si = Api ? (transa=='N'||transa=='T'?*Api:-*Api) : zero;
        RealTimesScalar(Cpr, Cpi, Bpr, Bpi, transb, m2, n2, sr, si, Bsize, 1, scalar_method);

/*----------------------------------------------------------------------------
 * Scalar product (M x K) * (1 x 1)
 *---------------------------------------------------------------------------- */

    } else if( scalarmultiply == 2 ) {
        sr = *Bpr;
		si = Bpi ? (transb=='N'||transb=='T'?*Bpi:-*Bpi) : zero;
        RealTimesScalar(Cpr, Cpi, Apr, Api, transa, m1, n1, sr, si, Asize, 1, scalar_method);

/*---------------------------------------------------------------------------------
 * Small matrix times small matrix (M x K) * (K x N) use inline code. Only use this
 * method if running in the 'SPEED' mode and M, K, N are all <= 4.
 *--------------------------------------------------------------------------------- */

    } else if( mtimesx_mode != MTIMESX_BLAS && mtimesx_mode != MTIMESX_MATLAB && m <= 4 && k <= 4 && n <= 4 ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: LOOPS (unrolled into inline switch statements)\n");
		}
		/* Form B elements, taking size and transb into account */
		switch( k ) {

		case 1: /* (m x 1)*(1 x n) */
			switch( n ) {
     		case 1:
				mexErrMsgTxt("Internal Error (m x 1)*(1 x 1), contact author.");
				break;

    		case 2: /* (m x 1)*(1 x 2) */
				Bpr11 = Bpr[0]; Bpr12 = Bpr[1];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; Bpi12 = Bpi[1];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1];
					}
				}
				break;
    		case 3: /* (m x 1)*(1 x 3) */
				Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2];
					}
				}
				break;
    		case 4: /* (m x 1)*(1 x 4) */
				Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2]; Bpr14 = Bpr[3];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2]; Bpi14 = Bpi[3];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2]; Bpi14 =-Bpi[3];
					}
				}
				break;
			}
			break;

		case 2: /* (m x 2)*(2 x n) */
			switch( n ) {
    		case 1: /* (m x 2)*(2 x 1) */
				Bpr11 = Bpr[0]; 
				Bpr21 = Bpr[1];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; 
						Bpi21 = Bpi[1];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0];
						Bpi21 =-Bpi[1];
					}
				}
				break;
    		case 2: /* (m x 2)*(2 x 2) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[2];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[3];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1];
   				    Bpr21 = Bpr[2]; Bpr22 = Bpr[3];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[2];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[3];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[2];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[3];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1];
   				        Bpi21 = Bpi[2]; Bpi22 = Bpi[3];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1];
   				        Bpi21 =-Bpi[2]; Bpi22 =-Bpi[3];
					}
				}
				break;
    		case 3: /* (m x 2)*(2 x 3) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[2]; Bpr13 = Bpr[4];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[3]; Bpr23 = Bpr[5];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2];
   				    Bpr21 = Bpr[3]; Bpr22 = Bpr[4]; Bpr23 = Bpr[5];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[2]; Bpi13 = Bpi[4];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[3]; Bpi23 = Bpi[5];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[2]; Bpi13 =-Bpi[4];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[3]; Bpi23 =-Bpi[5];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2];
   				        Bpi21 = Bpi[3]; Bpi22 = Bpi[4]; Bpi23 = Bpi[5];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2];
   				        Bpi21 =-Bpi[3]; Bpi22 =-Bpi[4]; Bpi23 =-Bpi[5];
					}
				}
				break;
    		case 4: /* (m x 2)*(2 x 4) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[2]; Bpr13 = Bpr[4]; Bpr14 = Bpr[6];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[3]; Bpr23 = Bpr[5]; Bpr24 = Bpr[7];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2]; Bpr14 = Bpr[3];
   				    Bpr21 = Bpr[4]; Bpr22 = Bpr[5]; Bpr23 = Bpr[6]; Bpr24 = Bpr[7];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[2]; Bpi13 = Bpi[4]; Bpi14 = Bpi[6];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[3]; Bpi23 = Bpi[5]; Bpi24 = Bpi[7];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[2]; Bpi13 =-Bpi[4]; Bpi14 =-Bpi[6];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[3]; Bpi23 =-Bpi[5]; Bpi24 =-Bpi[7];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2]; Bpi14 = Bpi[3];
   				        Bpi21 = Bpi[4]; Bpi22 = Bpi[5]; Bpi23 = Bpi[6]; Bpi24 = Bpi[7];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2]; Bpi14 =-Bpi[3];
   				        Bpi21 =-Bpi[4]; Bpi22 =-Bpi[5]; Bpi23 =-Bpi[6]; Bpi24 =-Bpi[7];
					}
				}
				break;
			}
			break;

		case 3: /* (m x 3)*(3 x n) */
			switch( n ) {
    		case 1: /* (m x 3)*(3 x 1) */
				Bpr11 = Bpr[0]; 
				Bpr21 = Bpr[1];
				Bpr31 = Bpr[2];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; 
						Bpi21 = Bpi[1];
						Bpi31 = Bpi[2];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0];
						Bpi21 =-Bpi[1];
						Bpi31 =-Bpi[2];
					}
				}
				break;
    		case 2: /* (m x 3)*(3 x 2) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[3];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[4];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[5];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1];
   				    Bpr21 = Bpr[2]; Bpr22 = Bpr[3];
   				    Bpr31 = Bpr[4]; Bpr32 = Bpr[5];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[3];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[4];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[5];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[3];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[4];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[5];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1];
   				        Bpi21 = Bpi[2]; Bpi22 = Bpi[3];
   				        Bpi31 = Bpi[4]; Bpi32 = Bpi[5];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1];
   				        Bpi21 =-Bpi[2]; Bpi22 =-Bpi[3];
   				        Bpi31 =-Bpi[4]; Bpi32 =-Bpi[5];
					}
				}
				break;
    		case 3: /* (m x 3)*(3 x 3) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[3]; Bpr13 = Bpr[6];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[4]; Bpr23 = Bpr[7];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[5]; Bpr33 = Bpr[8];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2];
   				    Bpr21 = Bpr[3]; Bpr22 = Bpr[4]; Bpr23 = Bpr[5];
   				    Bpr31 = Bpr[6]; Bpr32 = Bpr[7]; Bpr33 = Bpr[8];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[3]; Bpi13 = Bpi[6];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[4]; Bpi23 = Bpi[7];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[5]; Bpi33 = Bpi[8];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[3]; Bpi13 =-Bpi[6];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[4]; Bpi23 =-Bpi[7];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[5]; Bpi33 =-Bpi[8];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2];
   				        Bpi21 = Bpi[3]; Bpi22 = Bpi[4]; Bpi23 = Bpi[5];
   				        Bpi31 = Bpi[6]; Bpi32 = Bpi[7]; Bpi33 = Bpi[8];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2];
   				        Bpi21 =-Bpi[3]; Bpi22 =-Bpi[4]; Bpi23 =-Bpi[5];
   				        Bpi31 =-Bpi[6]; Bpi32 =-Bpi[7]; Bpi33 =-Bpi[8];
					}
				}
				break;
    		case 4: /* (m x 3)*(3 x 4) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[3]; Bpr13 = Bpr[6]; Bpr14 = Bpr[9];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[4]; Bpr23 = Bpr[7]; Bpr24 = Bpr[10];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[5]; Bpr33 = Bpr[8]; Bpr34 = Bpr[11];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2]; Bpr14 = Bpr[3];
   				    Bpr21 = Bpr[4]; Bpr22 = Bpr[5]; Bpr23 = Bpr[6]; Bpr24 = Bpr[7];
   				    Bpr31 = Bpr[8]; Bpr32 = Bpr[9]; Bpr33 = Bpr[10];Bpr34 = Bpr[11];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[3]; Bpi13 = Bpi[6]; Bpi14 = Bpi[9];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[4]; Bpi23 = Bpi[7]; Bpi24 = Bpi[10];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[5]; Bpi33 = Bpi[8]; Bpi34 = Bpi[11];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[3]; Bpi13 =-Bpi[6]; Bpi14 =-Bpi[9];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[4]; Bpi23 =-Bpi[7]; Bpi24 =-Bpi[10];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[5]; Bpi33 =-Bpi[8]; Bpi34 =-Bpi[11];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2]; Bpi14 = Bpi[3];
   				        Bpi21 = Bpi[4]; Bpi22 = Bpi[5]; Bpi23 = Bpi[6]; Bpi24 = Bpi[7];
   				        Bpi31 = Bpi[8]; Bpi32 = Bpi[9]; Bpi33 = Bpi[10];Bpi34 = Bpi[11];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2]; Bpi14 =-Bpi[3];
   				        Bpi21 =-Bpi[4]; Bpi22 =-Bpi[5]; Bpi23 =-Bpi[6]; Bpi24 =-Bpi[7];
   				        Bpi31 =-Bpi[8]; Bpi32 =-Bpi[9]; Bpi33 =-Bpi[10];Bpi34 =-Bpi[11];
					}
				}
				break;
			}
			break;

		case 4: /* (m x 4)*(4 x n) */
			switch( n ) {
    		case 1: /* (m x 4)*(4 x 1) */
				Bpr11 = Bpr[0]; 
				Bpr21 = Bpr[1];
				Bpr31 = Bpr[2];
				Bpr41 = Bpr[3];
				if( Bpi ) {
					if( transb == 'N' || transb == 'T' ) {
						Bpi11 = Bpi[0]; 
						Bpi21 = Bpi[1];
						Bpi31 = Bpi[2];
						Bpi41 = Bpi[3];
					} else { /* transb == 'G' || transb == 'C' */
						Bpi11 =-Bpi[0];
						Bpi21 =-Bpi[1];
						Bpi31 =-Bpi[2];
						Bpi41 =-Bpi[3];
					}
				}
				break;
    		case 2: /* (m x 4)*(4 x 2) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[4];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[5];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[6];
   				    Bpr41 = Bpr[3]; Bpr42 = Bpr[7];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1];
   				    Bpr21 = Bpr[2]; Bpr22 = Bpr[3];
   				    Bpr31 = Bpr[4]; Bpr32 = Bpr[5];
   				    Bpr41 = Bpr[6]; Bpr42 = Bpr[7];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[4];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[5];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[6];
   				        Bpi41 = Bpi[3]; Bpi42 = Bpi[7];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[4];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[5];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[6];
   				        Bpi41 =-Bpi[3]; Bpi42 =-Bpi[7];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1];
   				        Bpi21 = Bpi[2]; Bpi22 = Bpi[3];
   				        Bpi31 = Bpi[4]; Bpi32 = Bpi[5];
   				        Bpi41 = Bpi[6]; Bpi42 = Bpi[7];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1];
   				        Bpi21 =-Bpi[2]; Bpi22 =-Bpi[3];
   				        Bpi31 =-Bpi[4]; Bpi32 =-Bpi[5];
   				        Bpi41 =-Bpi[6]; Bpi42 =-Bpi[7];
					}
				}
				break;
    		case 3: /* (m x 4)*(4 x 3) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[4]; Bpr13 = Bpr[8];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[5]; Bpr23 = Bpr[9];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[6]; Bpr33 = Bpr[10];
   				    Bpr41 = Bpr[3]; Bpr42 = Bpr[7]; Bpr43 = Bpr[11];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2];
   				    Bpr21 = Bpr[3]; Bpr22 = Bpr[4]; Bpr23 = Bpr[5];
   				    Bpr31 = Bpr[6]; Bpr32 = Bpr[7]; Bpr33 = Bpr[8];
   				    Bpr41 = Bpr[9]; Bpr42 = Bpr[10];Bpr43 = Bpr[11];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[4]; Bpi13 = Bpi[8];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[5]; Bpi23 = Bpi[9];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[6]; Bpi33 = Bpi[10];
   				        Bpi41 = Bpi[3]; Bpi42 = Bpi[7]; Bpi43 = Bpi[11];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[4]; Bpi13 =-Bpi[8];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[5]; Bpi23 =-Bpi[9];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[6]; Bpi33 =-Bpi[10];
   				        Bpi41 =-Bpi[3]; Bpi42 =-Bpi[7]; Bpi43 =-Bpi[11];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2];
   				        Bpi21 = Bpi[3]; Bpi22 = Bpi[4]; Bpi23 = Bpi[5];
   				        Bpi31 = Bpi[6]; Bpi32 = Bpi[7]; Bpi33 = Bpi[8];
   				        Bpi41 = Bpi[9]; Bpi42 = Bpi[10];Bpi43 = Bpi[11];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2];
   				        Bpi21 =-Bpi[3]; Bpi22 =-Bpi[4]; Bpi23 =-Bpi[5];
   				        Bpi31 =-Bpi[6]; Bpi32 =-Bpi[7]; Bpi33 =-Bpi[8];
   				        Bpi41 =-Bpi[9]; Bpi42 =-Bpi[10];Bpi43 =-Bpi[11];
					}
				}
				break;
    		case 4: /* (m x 4)*(4 x 4) */
				if( transb == 'N' || transb == 'G' ) {
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[4]; Bpr13 = Bpr[8]; Bpr14 = Bpr[12];
   				    Bpr21 = Bpr[1]; Bpr22 = Bpr[5]; Bpr23 = Bpr[9]; Bpr24 = Bpr[13];
   				    Bpr31 = Bpr[2]; Bpr32 = Bpr[6]; Bpr33 = Bpr[10];Bpr34 = Bpr[14];
   				    Bpr41 = Bpr[3]; Bpr42 = Bpr[7]; Bpr43 = Bpr[11];Bpr44 = Bpr[15];
				} else { /* transa == 'T' || transa == 'C' */
   				    Bpr11 = Bpr[0]; Bpr12 = Bpr[1]; Bpr13 = Bpr[2]; Bpr14 = Bpr[3];
   				    Bpr21 = Bpr[4]; Bpr22 = Bpr[5]; Bpr23 = Bpr[6]; Bpr24 = Bpr[7];
   				    Bpr31 = Bpr[8]; Bpr32 = Bpr[9]; Bpr33 = Bpr[10];Bpr34 = Bpr[11];
   				    Bpr41 = Bpr[12];Bpr42 = Bpr[13];Bpr43 = Bpr[14];Bpr44 = Bpr[15];
				}
				if( Bpi ) {
					if( transb == 'N' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[4]; Bpi13 = Bpi[8]; Bpi14 = Bpi[12];
   				        Bpi21 = Bpi[1]; Bpi22 = Bpi[5]; Bpi23 = Bpi[9]; Bpi24 = Bpi[13];
   				        Bpi31 = Bpi[2]; Bpi32 = Bpi[6]; Bpi33 = Bpi[10];Bpi34 = Bpi[14];
   				        Bpi41 = Bpi[3]; Bpi42 = Bpi[7]; Bpi43 = Bpi[11];Bpi44 = Bpi[15];
					} else if( transb == 'G' ) {
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[4]; Bpi13 =-Bpi[8]; Bpi14 =-Bpi[12];
   				        Bpi21 =-Bpi[1]; Bpi22 =-Bpi[5]; Bpi23 =-Bpi[9]; Bpi24 =-Bpi[13];
   				        Bpi31 =-Bpi[2]; Bpi32 =-Bpi[6]; Bpi33 =-Bpi[10];Bpi34 =-Bpi[14];
   				        Bpi41 =-Bpi[3]; Bpi42 =-Bpi[7]; Bpi43 =-Bpi[11];Bpi44 =-Bpi[15];
					} else if( transb == 'T' ) {
   				        Bpi11 = Bpi[0]; Bpi12 = Bpi[1]; Bpi13 = Bpi[2]; Bpi14 = Bpi[3];
   				        Bpi21 = Bpi[4]; Bpi22 = Bpi[5]; Bpi23 = Bpi[6]; Bpi24 = Bpi[7];
   				        Bpi31 = Bpi[8]; Bpi32 = Bpi[9]; Bpi33 = Bpi[10];Bpi34 = Bpi[11];
   				        Bpi41 = Bpi[12];Bpi42 = Bpi[13];Bpi43 = Bpi[14];Bpi44 = Bpi[15];
					} else {/* transb == 'C' */
   				        Bpi11 =-Bpi[0]; Bpi12 =-Bpi[1]; Bpi13 =-Bpi[2]; Bpi14 =-Bpi[3];
   				        Bpi21 =-Bpi[4]; Bpi22 =-Bpi[5]; Bpi23 =-Bpi[6]; Bpi24 =-Bpi[7];
   				        Bpi31 =-Bpi[8]; Bpi32 =-Bpi[9]; Bpi33 =-Bpi[10];Bpi34 =-Bpi[11];
   				        Bpi41 =-Bpi[12];Bpi42 =-Bpi[13];Bpi43 =-Bpi[14];Bpi44 =-Bpi[15];
					}
				}
				break;
			}
			break;
		}
		/* Form A elements and do the multiply */
		switch( m ) {
		case 1: /* (1 x k)*(k x n) */
			switch( k ) {
     		case 1: /* (1 x 1)*(1 x n) */
				mexErrMsgTxt("Internal Error (1 x 1)*(1 x n), contact author.");
				break;

			case 2: /* (1 x 2)*(2 x n) */
				Apr11 = Apr[0]; Apr12 = Apr[1];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0]; Api12 = Api[1];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0]; Api12 =-Api[1];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24
							   - Api11 * Bpi14 - Api12 * Bpi24;
						Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24
							   + Api11 * Bpr14 + Api12 * Bpr24;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23
							   - Api11 * Bpi13 - Api12 * Bpi23;
						Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23
							   + Api11 * Bpr13 + Api12 * Bpr23;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22
							   - Api11 * Bpi12 - Api12 * Bpi22;
						Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22
							   + Api11 * Bpr12 + Api12 * Bpr22;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21
							   - Api11 * Bpi11 - Api12 * Bpi21;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21
							   + Api11 * Bpr11 + Api12 * Bpr21;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Api11 * Bpr14 + Api12 * Bpr24;
					    case 3:
						    Cpi[2] = Api11 * Bpr13 + Api12 * Bpr23;
					    case 2:
						    Cpi[1] = Api11 * Bpr12 + Api12 * Bpr22;
					    case 1:
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24;
					    case 3:
						    Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23;
					    case 2:
						    Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22;
					    case 1:
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21;
					    }
					}
				}
				break;

			case 3: /* (1 x 3)*(3 x n) */
				Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34;
						Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33;
						Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32;
						Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
					    case 3:
						    Cpi[2] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
					    case 2:
						    Cpi[1] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
					    case 1:
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34;
					    case 3:
						    Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33;
					    case 2:
						    Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32;
					    case 1:
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31;
					    }
					}
				}
				break;

			case 4: /* (1 x 4)*(4 x n) */
				Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2]; Apr14 = Apr[3];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2]; Api14 = Api[3];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2]; Api14 =-Api[3];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34 - Api14 * Bpi44;
						Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43 
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33 - Api14 * Bpi43;
						Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32 - Api14 * Bpi42;
						Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31 - Api14 * Bpi41;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[3] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44;
					case 3:
						Cpr[2] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43;
					case 2:
						Cpr[1] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42;
					case 1:
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					    case 3:
						    Cpi[2] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					    case 2:
						    Cpi[1] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					    case 1:
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[3] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44;
					    case 3:
						    Cpi[2] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43;
					    case 2:
						    Cpi[1] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42;
					    case 1:
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41;
					    }
					}
				}
				break;
			}
			break;

		case 2: /* (2 x k)*(k x n) */
			switch( k ) {
     		case 1: /* (2 x 1)*(1 x n) */
				Apr11 = Apr[0];
				Apr21 = Apr[1];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0];
						Api21 = Api[1];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0];
						Api21 =-Api[1];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14
							   - Api21 * Bpi14;
						Cpr[6] = Apr11 * Bpr14
							   - Api11 * Bpi14;
						Cpi[7] = Apr21 * Bpi14
							   + Api21 * Bpr14;
						Cpi[6] = Apr11 * Bpi14
							   + Api11 * Bpr14;
					case 3:
						Cpr[5] = Apr21 * Bpr13
							   - Api21 * Bpi13;
						Cpr[4] = Apr11 * Bpr13
							   - Api11 * Bpi13;
						Cpi[5] = Apr21 * Bpi13
							   + Api21 * Bpr13;
						Cpi[4] = Apr11 * Bpi13
							   + Api11 * Bpr13;
					case 2:
						Cpr[3] = Apr21 * Bpr12
							   - Api21 * Bpi12;
						Cpr[2] = Apr11 * Bpr12
							   - Api11 * Bpi12;
						Cpi[3] = Apr21 * Bpi12
							   + Api21 * Bpr12;
						Cpi[2] = Apr11 * Bpi12
							   + Api11 * Bpr12;
					case 1:
						Cpr[1] = Apr21 * Bpr11
							   - Api21 * Bpi11;
						Cpr[0] = Apr11 * Bpr11
							   - Api11 * Bpi11;
						Cpi[1] = Apr21 * Bpi11
							   + Api21 * Bpr11;
						Cpi[0] = Apr11 * Bpi11
							   + Api11 * Bpr11;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14;
						Cpr[6] = Apr11 * Bpr14;
					case 3:
						Cpr[5] = Apr21 * Bpr13;
						Cpr[4] = Apr11 * Bpr13;
					case 2:
						Cpr[3] = Apr21 * Bpr12;
						Cpr[2] = Apr11 * Bpr12;
					case 1:
						Cpr[1] = Apr21 * Bpr11;
						Cpr[0] = Apr11 * Bpr11;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[7] = Api21 * Bpr14;
						    Cpi[6] = Api11 * Bpr14;
					    case 3:
						    Cpi[5] = Api21 * Bpr13;
						    Cpi[4] = Api11 * Bpr13;
					    case 2:
						    Cpi[3] = Api21 * Bpr12;
						    Cpi[2] = Api11 * Bpr12;
					    case 1:
						    Cpi[1] = Api21 * Bpr11;
						    Cpi[0] = Api11 * Bpr11;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[7] = Apr21 * Bpi14;
						    Cpi[6] = Apr11 * Bpi14;
					    case 3:
						    Cpi[5] = Apr21 * Bpi13;
						    Cpi[4] = Apr11 * Bpi13;
					    case 2:
						    Cpi[3] = Apr21 * Bpi12;
						    Cpi[2] = Apr11 * Bpi12;
					    case 1:
						    Cpi[1] = Apr21 * Bpi11;
						    Cpi[0] = Apr11 * Bpi11;
					    }
					}
				}
				break;

			case 2: /* (2 x 2)*(2 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[2];
  				    Apr21 = Apr[1]; Apr22 = Apr[3];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1];
  				    Apr21 = Apr[2]; Apr22 = Apr[3];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[2];
                        Api21 = Api[1]; Api22 = Api[3];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1];
                        Api21 = Api[2]; Api22 = Api[3];
                    } else if( transa == 'C' ) {
                        Api11 = -Api[0]; Api12 = -Api[1];
                        Api21 = -Api[2]; Api22 = -Api[3];
                    } else {/* transa == 'G' */
                        Api11 = -Api[0]; Api12 = -Api[2];
                        Api21 = -Api[1]; Api22 = -Api[3];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24
							   - Api21 * Bpi14 - Api22 * Bpi24;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24
							   - Api11 * Bpi14 - Api12 * Bpi24;
						Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24
							   + Api21 * Bpr14 + Api22 * Bpr24;
						Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24
							   + Api11 * Bpr14 + Api12 * Bpr24;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23
							   - Api21 * Bpi13 - Api22 * Bpi23;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23
							   - Api11 * Bpi13 - Api12 * Bpi23;
						Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23
							   + Api21 * Bpr13 + Api22 * Bpr23;
						Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23
							   + Api11 * Bpr13 + Api12 * Bpr23;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22
							   - Api21 * Bpi12 - Api22 * Bpi22;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22
							   - Api11 * Bpi12 - Api12 * Bpi22;
						Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22
							   + Api21 * Bpr12 + Api22 * Bpr22;
						Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22
							   + Api11 * Bpr12 + Api12 * Bpr22;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21
							   - Api21 * Bpi11 - Api22 * Bpi21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21
							   - Api11 * Bpi11 - Api12 * Bpi21;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21
							   + Api21 * Bpr11 + Api22 * Bpr21;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21
							   + Api11 * Bpr11 + Api12 * Bpr21;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[7] = Api21 * Bpr14 + Api22 * Bpr24;
							Cpi[6] = Api11 * Bpr14 + Api12 * Bpr24;
						case 3:
							Cpi[5] = Api21 * Bpr13 + Api22 * Bpr23;
							Cpi[4] = Api11 * Bpr13 + Api12 * Bpr23;
						case 2:
							Cpi[3] = Api21 * Bpr12 + Api22 * Bpr22;
							Cpi[2] = Api11 * Bpr12 + Api12 * Bpr22;
						case 1:
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24;
							Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24;
						case 3:
							Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23;
							Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23;
						case 2:
							Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22;
							Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22;
						case 1:
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21;
					    }
					}
				}
				break;

			case 3: /* (2 x 3)*(3 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[2]; Apr13 = Apr[4];
  				    Apr21 = Apr[1]; Apr22 = Apr[3]; Apr23 = Apr[5];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2];
  				    Apr21 = Apr[3]; Apr22 = Apr[4]; Apr23 = Apr[5];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[2]; Api13 = Api[4];
                        Api21 = Api[1]; Api22 = Api[3]; Api23 = Api[5];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2];
                        Api21 = Api[3]; Api22 = Api[4]; Api23 = Api[5];
                    } else if( transa == 'C' ) {
                        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2];
                        Api21 =-Api[3]; Api22 =-Api[4]; Api23 =-Api[5];
                    } else {/* transa == 'G' */
                        Api11 =-Api[0]; Api12 =-Api[2]; Api13 =-Api[4];
                        Api21 =-Api[1]; Api22 =-Api[3]; Api23 =-Api[5];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34;
						Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
						Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33;
						Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
						Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32;
						Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
						Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[7] = Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
							Cpi[6] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
						case 3:
							Cpi[5] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
							Cpi[4] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
						case 2:
							Cpi[3] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
							Cpi[2] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
						case 1:
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34;
							Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34;
						case 3:
							Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33;
							Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33;
						case 2:
							Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32;
							Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32;
						case 1:
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31;
					    }
					}
				}
				break;

			case 4: /* (2 x 4)*(4 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[2]; Apr13 = Apr[4]; Apr14 = Apr[6];
				    Apr21 = Apr[1]; Apr22 = Apr[3]; Apr23 = Apr[5]; Apr24 = Apr[7];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2]; Apr14 = Apr[3];
				    Apr21 = Apr[4]; Apr22 = Apr[5]; Apr23 = Apr[6]; Apr24 = Apr[7];
				}
				if( Api ) {
				    if( transa == 'N' ) {
				        Api11 = Api[0]; Api12 = Api[2]; Api13 = Api[4]; Api14 = Api[6];
				        Api21 = Api[1]; Api22 = Api[3]; Api23 = Api[5]; Api24 = Api[7];
				    } else if( transa == 'T' ) {
				        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2]; Api14 = Api[3];
				        Api21 = Api[4]; Api22 = Api[5]; Api23 = Api[6]; Api24 = Api[7];
					} else if( transa == 'G' ) {
				        Api11 =-Api[0]; Api12 =-Api[2]; Api13 =-Api[4]; Api14 =-Api[6];
				        Api21 =-Api[1]; Api22 =-Api[3]; Api23 =-Api[5]; Api24 =-Api[7];
				    } else { /* transa == 'C' */
				        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2]; Api14 =-Api[3];
				        Api21 =-Api[4]; Api22 =-Api[5]; Api23 =-Api[6]; Api24 =-Api[7];
				    }
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34 - Api24 * Bpi44;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34 - Api14 * Bpi44;
						Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33 - Api24 * Bpi43;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33 - Api14 * Bpi43;
						Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32 - Api24 * Bpi42;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32 - Api14 * Bpi42;
						Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31 - Api24 * Bpi41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31 - Api14 * Bpi41;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[7] = Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44;
						Cpr[6] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44;
					case 3:
						Cpr[5] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43;
						Cpr[4] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43;
					case 2:
						Cpr[3] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42;
						Cpr[2] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42;
					case 1:
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[7] = Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						    Cpi[6] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					    case 3:
						    Cpi[5] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						    Cpi[4] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					    case 2:
						    Cpi[3] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						    Cpi[2] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					    case 1:
						    Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[7] = Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44;
						    Cpi[6] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44;
					    case 3:
						    Cpi[5] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43;
						    Cpi[4] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43;
					    case 2:
						    Cpi[3] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42;
						    Cpi[2] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42;
					    case 1:
						    Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41;
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41;
					    }
					}
				}
				break;
			}
			break;

		case 3: /* (3 x k)*(k x n) */
			switch( k ) {
     		case 1: /* (3 x 1)*(1 x n) */
				Apr11 = Apr[0];
				Apr21 = Apr[1];
				Apr31 = Apr[2];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0];
						Api21 = Api[1];
						Api31 = Api[2];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0];
						Api21 =-Api[1];
						Api31 =-Api[2];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14
							   - Api31 * Bpi14;
						Cpr[10]= Apr21 * Bpr14
							   - Api21 * Bpi14;
						Cpr[9] = Apr11 * Bpr14
							   - Api11 * Bpi14;
						Cpi[11]= Apr31 * Bpi14
							   + Api31 * Bpr14;
						Cpi[10]= Apr21 * Bpi14
							   + Api21 * Bpr14;
						Cpi[9] = Apr11 * Bpi14
							   + Api11 * Bpr14;
					case 3:
						Cpr[8] = Apr31 * Bpr13
							   - Api31 * Bpi13;
						Cpr[7] = Apr21 * Bpr13
							   - Api21 * Bpi13;
						Cpr[6] = Apr11 * Bpr13
							   - Api11 * Bpi13;
						Cpi[8] = Apr31 * Bpi13
							   + Api31 * Bpr13;
						Cpi[7] = Apr21 * Bpi13
							   + Api21 * Bpr13;
						Cpi[6] = Apr11 * Bpi13
							   + Api11 * Bpr13;
					case 2:
						Cpr[5] = Apr31 * Bpr12
							   - Api31 * Bpi12;
						Cpr[4] = Apr21 * Bpr12
							   - Api21 * Bpi12;
						Cpr[3] = Apr11 * Bpr12
							   - Api11 * Bpi12;
						Cpi[5] = Apr31 * Bpi12
							   + Api31 * Bpr12;
						Cpi[4] = Apr21 * Bpi12
							   + Api21 * Bpr12;
						Cpi[3] = Apr11 * Bpi12
							   + Api11 * Bpr12;
					case 1:
						Cpr[2] = Apr31 * Bpr11
							   - Api31 * Bpi11;
						Cpr[1] = Apr21 * Bpr11
							   - Api21 * Bpi11;
						Cpr[0] = Apr11 * Bpr11
							   - Api11 * Bpi11;
						Cpi[2] = Apr31 * Bpi11
							   + Api31 * Bpr11;
						Cpi[1] = Apr21 * Bpi11
							   + Api21 * Bpr11;
						Cpi[0] = Apr11 * Bpi11
							   + Api11 * Bpr11;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14;
						Cpr[10]= Apr21 * Bpr14;
						Cpr[9] = Apr11 * Bpr14;
					case 3:
						Cpr[8] = Apr31 * Bpr13;
						Cpr[7] = Apr21 * Bpr13;
						Cpr[6] = Apr11 * Bpr13;
					case 2:
						Cpr[5] = Apr31 * Bpr12;
						Cpr[4] = Apr21 * Bpr12;
						Cpr[3] = Apr11 * Bpr12;
					case 1:
						Cpr[2] = Apr31 * Bpr11;
						Cpr[1] = Apr21 * Bpr11;
						Cpr[0] = Apr11 * Bpr11;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[11]= Api31 * Bpr14;
						    Cpi[10]= Api21 * Bpr14;
						    Cpi[9] = Api11 * Bpr14;
					    case 3:
						    Cpi[8] = Api31 * Bpr13;
						    Cpi[7] = Api21 * Bpr13;
						    Cpi[6] = Api11 * Bpr13;
					    case 2:
						    Cpi[5] = Api31 * Bpr12;
						    Cpi[4] = Api21 * Bpr12;
						    Cpi[3] = Api11 * Bpr12;
					    case 1:
						    Cpi[2] = Api31 * Bpr11;
						    Cpi[1] = Api21 * Bpr11;
						    Cpi[0] = Api11 * Bpr11;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[11]= Apr31 * Bpi14;
						    Cpi[10]= Apr21 * Bpi14;
						    Cpi[9] = Apr11 * Bpi14;
					    case 3:
						    Cpi[8] = Apr31 * Bpi13;
						    Cpi[7] = Apr21 * Bpi13;
						    Cpi[6] = Apr11 * Bpi13;
					    case 2:
						    Cpi[5] = Apr31 * Bpi12;
						    Cpi[4] = Apr21 * Bpi12;
						    Cpi[3] = Apr11 * Bpi12;
					    case 1:
						    Cpi[2] = Apr31 * Bpi11;
						    Cpi[1] = Apr21 * Bpi11;
						    Cpi[0] = Apr11 * Bpi11;
					    }
					}
				}
				break;

			case 2: /* (3 x 2)*(2 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[3];
  				    Apr21 = Apr[1]; Apr22 = Apr[4];
  				    Apr31 = Apr[2]; Apr32 = Apr[5];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1];
  				    Apr21 = Apr[2]; Apr22 = Apr[3];
  				    Apr31 = Apr[4]; Apr32 = Apr[5];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[3];
                        Api21 = Api[1]; Api22 = Api[4];
                        Api31 = Api[2]; Api32 = Api[5];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1];
                        Api21 = Api[2]; Api22 = Api[3];
                        Api31 = Api[4]; Api32 = Api[5];
                    } else if( transa == 'C' ) {
                        Api11 = -Api[0]; Api12 = -Api[1];
                        Api21 = -Api[2]; Api22 = -Api[3];
                        Api31 = -Api[4]; Api32 = -Api[5];
                    } else {/* transa == 'G' */
                        Api11 = -Api[0]; Api12 = -Api[3];
                        Api21 = -Api[1]; Api22 = -Api[4];
                        Api31 = -Api[2]; Api32 = -Api[5];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24
							   - Api31 * Bpi14 - Api32 * Bpi24;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24
							   - Api21 * Bpi14 - Api22 * Bpi24;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24
							   - Api11 * Bpi14 - Api12 * Bpi24;
						Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24
							   + Api31 * Bpr14 + Api32 * Bpr24;
						Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24
							   + Api21 * Bpr14 + Api22 * Bpr24;
						Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24
							   + Api11 * Bpr14 + Api12 * Bpr24;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23
							   - Api31 * Bpi13 - Api32 * Bpi23;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23
							   - Api21 * Bpi13 - Api22 * Bpi23;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23
							   - Api11 * Bpi13 - Api12 * Bpi23;
						Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23
							   + Api31 * Bpr13 + Api32 * Bpr23;
						Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23
							   + Api21 * Bpr13 + Api22 * Bpr23;
						Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23
							   + Api11 * Bpr13 + Api12 * Bpr23;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22
							   - Api31 * Bpi12 - Api32 * Bpi22;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22
							   - Api21 * Bpi12 - Api22 * Bpi22;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22
							   - Api11 * Bpi12 - Api12 * Bpi22;
						Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22
							   + Api31 * Bpr12 + Api32 * Bpr22;
						Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22
							   + Api21 * Bpr12 + Api22 * Bpr22;
						Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22
							   + Api11 * Bpr12 + Api12 * Bpr22;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21
							   - Api31 * Bpi11 - Api32 * Bpi21;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21
							   - Api21 * Bpi11 - Api22 * Bpi21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21
							   - Api11 * Bpi11 - Api12 * Bpi21;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21
							   + Api31 * Bpr11 + Api32 * Bpr21;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21
							   + Api21 * Bpr11 + Api22 * Bpr21;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21
							   + Api11 * Bpr11 + Api12 * Bpr21;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[11]= Api31 * Bpr14 + Api32 * Bpr24;
							Cpi[10]= Api21 * Bpr14 + Api22 * Bpr24;
							Cpi[9] = Api11 * Bpr14 + Api12 * Bpr24;
						case 3:
							Cpi[8] = Api31 * Bpr13 + Api32 * Bpr23;
							Cpi[7] = Api21 * Bpr13 + Api22 * Bpr23;
							Cpi[6] = Api11 * Bpr13 + Api12 * Bpr23;
						case 2:
							Cpi[5] = Api31 * Bpr12 + Api32 * Bpr22;
							Cpi[4] = Api21 * Bpr12 + Api22 * Bpr22;
							Cpi[3] = Api11 * Bpr12 + Api12 * Bpr22;
						case 1:
							Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21;
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24;
							Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24;
							Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24;
						case 3:
							Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23;
							Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23;
							Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23;
						case 2:
							Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22;
							Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22;
							Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22;
						case 1:
							Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21;
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21;
					    }
					}
				}
				break;

			case 3: /* (3 x 3)*(3 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[3]; Apr13 = Apr[6];
  				    Apr21 = Apr[1]; Apr22 = Apr[4]; Apr23 = Apr[7];
  				    Apr31 = Apr[2]; Apr32 = Apr[5]; Apr33 = Apr[8];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2];
  				    Apr21 = Apr[3]; Apr22 = Apr[4]; Apr23 = Apr[5];
  				    Apr31 = Apr[6]; Apr32 = Apr[7]; Apr33 = Apr[8];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[3]; Api13 = Api[6];
                        Api21 = Api[1]; Api22 = Api[4]; Api23 = Api[7];
                        Api31 = Api[2]; Api32 = Api[5]; Api33 = Api[8];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2];
                        Api21 = Api[3]; Api22 = Api[4]; Api23 = Api[5];
                        Api31 = Api[6]; Api32 = Api[7]; Api33 = Api[8];
                    } else if( transa == 'C' ) {
                        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2];
                        Api21 =-Api[3]; Api22 =-Api[4]; Api23 =-Api[5];
                        Api31 =-Api[6]; Api32 =-Api[7]; Api33 =-Api[8];
                    } else {/* transa == 'G' */
                        Api11 =-Api[0]; Api12 =-Api[3]; Api13 =-Api[6];
                        Api21 =-Api[1]; Api22 =-Api[4]; Api23 =-Api[7];
                        Api31 =-Api[2]; Api32 =-Api[5]; Api33 =-Api[8];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34
							   - Api31 * Bpi14 - Api32 * Bpi24 - Api33 * Bpi34;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34;
						Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34
							   + Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34;
						Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
						Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33
							   - Api31 * Bpi13 - Api32 * Bpi23 - Api33 * Bpi33;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33;
						Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33
							   + Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33;
						Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
						Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32
							   - Api31 * Bpi12 - Api32 * Bpi22 - Api33 * Bpi32;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32;
						Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32
							   + Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32;
						Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
						Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31
							   - Api31 * Bpi11 - Api32 * Bpi21 - Api33 * Bpi31;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31
							   + Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[11]= Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34;
							Cpi[10]= Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
							Cpi[9] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
						case 3:
							Cpi[8] = Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33;
							Cpi[7] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
							Cpi[6] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
						case 2:
							Cpi[5] = Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32;
							Cpi[4] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
							Cpi[3] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
						case 1:
							Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31;
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34;
							Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34;
							Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34;
						case 3:
							Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33;
							Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33;
							Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33;
						case 2:
							Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32;
							Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32;
							Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32;
						case 1:
							Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31;
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31;
					    }
					}
				}
				break;

			case 4: /* (3 x 4)*(4 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[3]; Apr13 = Apr[6]; Apr14 = Apr[9];
				    Apr21 = Apr[1]; Apr22 = Apr[4]; Apr23 = Apr[7]; Apr24 = Apr[10];
				    Apr31 = Apr[2]; Apr32 = Apr[5]; Apr33 = Apr[8]; Apr34 = Apr[11];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2]; Apr14 = Apr[3];
				    Apr21 = Apr[4]; Apr22 = Apr[5]; Apr23 = Apr[6]; Apr24 = Apr[7];
				    Apr31 = Apr[8]; Apr32 = Apr[9]; Apr33 = Apr[10];Apr34 = Apr[11];
				}
				if( Api ) {
				    if( transa == 'N' ) {
				        Api11 = Api[0]; Api12 = Api[3]; Api13 = Api[6]; Api14 = Api[9];
				        Api21 = Api[1]; Api22 = Api[4]; Api23 = Api[7]; Api24 = Api[10];
				        Api31 = Api[2]; Api32 = Api[5]; Api33 = Api[8]; Api34 = Api[11];
				    } else if( transa == 'T' ) {
				        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2]; Api14 = Api[3];
				        Api21 = Api[4]; Api22 = Api[5]; Api23 = Api[6]; Api24 = Api[7];
				        Api31 = Api[8]; Api32 = Api[9]; Api33 = Api[10];Api34 = Api[11];
					} else if( transa == 'G' ) {
				        Api11 =-Api[0]; Api12 =-Api[3]; Api13 =-Api[6]; Api14 =-Api[9];
				        Api21 =-Api[1]; Api22 =-Api[4]; Api23 =-Api[7]; Api24 =-Api[10];
				        Api31 =-Api[2]; Api32 =-Api[5]; Api33 =-Api[8]; Api34 =-Api[11];
				    } else { /* transa == 'C' */
				        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2]; Api14 =-Api[3];
				        Api21 =-Api[4]; Api22 =-Api[5]; Api23 =-Api[6]; Api24 =-Api[7];
				        Api31 =-Api[8]; Api32 =-Api[9]; Api33 =-Api[10];Api34 =-Api[11];
				    }
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34 + Apr34 * Bpr44
							   - Api31 * Bpi14 - Api32 * Bpi24 - Api33 * Bpi34 - Api34 * Bpi44;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34 - Api24 * Bpi44;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34 - Api14 * Bpi44;
						Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34 + Apr34 * Bpi44
							   + Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34 + Api34 * Bpr44;
						Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33 + Apr34 * Bpr43
							   - Api31 * Bpi13 - Api32 * Bpi23 - Api33 * Bpi33 - Api34 * Bpi43;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33 - Api24 * Bpi43;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33 - Api14 * Bpi43;
						Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33 + Apr34 * Bpi43
							   + Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33 + Api34 * Bpr43;
						Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32 + Apr34 * Bpr42
							   - Api31 * Bpi12 - Api32 * Bpi22 - Api33 * Bpi32 - Api34 * Bpi42;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32 - Api24 * Bpi42;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32 - Api14 * Bpi42;
						Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32 + Apr34 * Bpi42
							   + Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32 + Api34 * Bpr42;
						Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31 + Apr34 * Bpr41
							   - Api31 * Bpi11 - Api32 * Bpi21 - Api33 * Bpi31 - Api34 * Bpi41;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31 - Api24 * Bpi41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31 - Api14 * Bpi41;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31 + Apr34 * Bpi41
							   + Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31 + Api34 * Bpr41;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[11]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34 + Apr34 * Bpr44;
						Cpr[10]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44;
						Cpr[9] = Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44;
					case 3:
						Cpr[8] = Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33 + Apr34 * Bpr43;
						Cpr[7] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43;
						Cpr[6] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43;
					case 2:
						Cpr[5] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32 + Apr34 * Bpr42;
						Cpr[4] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42;
						Cpr[3] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42;
					case 1:
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31 + Apr34 * Bpr41;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[11]= Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34 + Api34 * Bpr44;
						    Cpi[10]= Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						    Cpi[9] = Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					    case 3:
						    Cpi[8] = Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33 + Api34 * Bpr43;
						    Cpi[7] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						    Cpi[6] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					    case 2:
						    Cpi[5] = Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32 + Api34 * Bpr42;
						    Cpi[4] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						    Cpi[3] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					    case 1:
						    Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31 + Api34 * Bpr41;
						    Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[11]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34 + Apr34 * Bpi44;
						    Cpi[10]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44;
						    Cpi[9] = Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44;
					    case 3:
						    Cpi[8] = Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33 + Apr34 * Bpi43;
						    Cpi[7] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43;
						    Cpi[6] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43;
					    case 2:
						    Cpi[5] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32 + Apr34 * Bpi42;
						    Cpi[4] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42;
						    Cpi[3] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42;
					    case 1:
						    Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31 + Apr34 * Bpi41;
						    Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41;
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41;
					    }
					}
				}
				break;
			}
			break;

		case 4: /* (4 x k)*(k x n) */
			switch( k ) {
     		case 1: /* (4 x 1)*(1 x n) */
				Apr11 = Apr[0];
				Apr21 = Apr[1];
				Apr31 = Apr[2];
				Apr41 = Apr[3];
				if( Api ) {
					if( transa == 'N' || transa == 'T' ) {
        				Api11 = Api[0];
						Api21 = Api[1];
						Api31 = Api[2];
						Api41 = Api[3];
					} else { /* transa == 'G' || transa == 'C' */
        				Api11 =-Api[0];
						Api21 =-Api[1];
						Api31 =-Api[2];
						Api41 =-Api[3];
					}
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14
							   - Api41 * Bpi14;
						Cpr[14]= Apr31 * Bpr14
							   - Api31 * Bpi14;
						Cpr[13]= Apr21 * Bpr14
							   - Api21 * Bpi14;
						Cpr[12]= Apr11 * Bpr14
							   - Api11 * Bpi14;
						Cpi[15]= Apr41 * Bpi14
							   + Api41 * Bpr14;
						Cpi[14]= Apr31 * Bpi14
							   + Api31 * Bpr14;
						Cpi[13]= Apr21 * Bpi14
							   + Api21 * Bpr14;
						Cpi[12]= Apr11 * Bpi14
							   + Api11 * Bpr14;
					case 3:
						Cpr[11]= Apr41 * Bpr13
							   - Api41 * Bpi13;
						Cpr[10]= Apr31 * Bpr13
							   - Api31 * Bpi13;
						Cpr[9] = Apr21 * Bpr13
							   - Api21 * Bpi13;
						Cpr[8] = Apr11 * Bpr13
							   - Api11 * Bpi13;
						Cpi[11]= Apr41 * Bpi13
							   + Api41 * Bpr13;
						Cpi[10]= Apr31 * Bpi13
							   + Api31 * Bpr13;
						Cpi[9] = Apr21 * Bpi13
							   + Api21 * Bpr13;
						Cpi[8] = Apr11 * Bpi13
							   + Api11 * Bpr13;
					case 2:
						Cpr[7] = Apr41 * Bpr12
							   - Api41 * Bpi12;
						Cpr[6] = Apr31 * Bpr12
							   - Api31 * Bpi12;
						Cpr[5] = Apr21 * Bpr12
							   - Api21 * Bpi12;
						Cpr[4] = Apr11 * Bpr12
							   - Api11 * Bpi12;
						Cpi[7] = Apr41 * Bpi12
							   + Api41 * Bpr12;
						Cpi[6] = Apr31 * Bpi12
							   + Api31 * Bpr12;
						Cpi[5] = Apr21 * Bpi12
							   + Api21 * Bpr12;
						Cpi[4] = Apr11 * Bpi12
							   + Api11 * Bpr12;
					case 1:
						Cpr[3] = Apr41 * Bpr11
							   - Api41 * Bpi11;
						Cpr[2] = Apr31 * Bpr11
							   - Api31 * Bpi11;
						Cpr[1] = Apr21 * Bpr11
							   - Api21 * Bpi11;
						Cpr[0] = Apr11 * Bpr11
							   - Api11 * Bpi11;
						Cpi[3] = Apr41 * Bpi11
							   + Api41 * Bpr11;
						Cpi[2] = Apr31 * Bpi11
							   + Api31 * Bpr11;
						Cpi[1] = Apr21 * Bpi11
							   + Api21 * Bpr11;
						Cpi[0] = Apr11 * Bpi11
							   + Api11 * Bpr11;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14;
						Cpr[14]= Apr31 * Bpr14;
						Cpr[13]= Apr21 * Bpr14;
						Cpr[12]= Apr11 * Bpr14;
					case 3:
						Cpr[11]= Apr41 * Bpr13;
						Cpr[10]= Apr31 * Bpr13;
						Cpr[9] = Apr21 * Bpr13;
						Cpr[8] = Apr11 * Bpr13;
					case 2:
						Cpr[7] = Apr41 * Bpr12;
						Cpr[6] = Apr31 * Bpr12;
						Cpr[5] = Apr21 * Bpr12;
						Cpr[4] = Apr11 * Bpr12;
					case 1:
						Cpr[3] = Apr41 * Bpr11;
						Cpr[2] = Apr31 * Bpr11;
						Cpr[1] = Apr21 * Bpr11;
						Cpr[0] = Apr11 * Bpr11;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[15]= Api41 * Bpr14;
						    Cpi[14]= Api31 * Bpr14;
						    Cpi[13]= Api21 * Bpr14;
						    Cpi[12]= Api11 * Bpr14;
					    case 3:
						    Cpi[11]= Api41 * Bpr13;
						    Cpi[10]= Api31 * Bpr13;
						    Cpi[9] = Api21 * Bpr13;
						    Cpi[8] = Api11 * Bpr13;
					    case 2:
						    Cpi[7] = Api41 * Bpr12;
						    Cpi[6] = Api31 * Bpr12;
						    Cpi[5] = Api21 * Bpr12;
						    Cpi[4] = Api11 * Bpr12;
					    case 1:
						    Cpi[3] = Api41 * Bpr11;
						    Cpi[2] = Api31 * Bpr11;
						    Cpi[1] = Api21 * Bpr11;
						    Cpi[0] = Api11 * Bpr11;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[15]= Apr41 * Bpi14;
						    Cpi[14]= Apr31 * Bpi14;
						    Cpi[13]= Apr21 * Bpi14;
						    Cpi[12]= Apr11 * Bpi14;
					    case 3:
						    Cpi[11]= Apr41 * Bpi13;
						    Cpi[10]= Apr31 * Bpi13;
						    Cpi[9] = Apr21 * Bpi13;
						    Cpi[8] = Apr11 * Bpi13;
					    case 2:
						    Cpi[7] = Apr41 * Bpi12;
						    Cpi[6] = Apr31 * Bpi12;
						    Cpi[5] = Apr21 * Bpi12;
						    Cpi[4] = Apr11 * Bpi12;
					    case 1:
						    Cpi[3] = Apr41 * Bpi11;
						    Cpi[2] = Apr31 * Bpi11;
						    Cpi[1] = Apr21 * Bpi11;
						    Cpi[0] = Apr11 * Bpi11;
					    }
					}
				}
				break;

			case 2: /* (4 x 2)*(2 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[4];
  				    Apr21 = Apr[1]; Apr22 = Apr[5];
  				    Apr31 = Apr[2]; Apr32 = Apr[6];
  				    Apr41 = Apr[3]; Apr42 = Apr[7];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1];
  				    Apr21 = Apr[2]; Apr22 = Apr[3];
  				    Apr31 = Apr[4]; Apr32 = Apr[5];
  				    Apr41 = Apr[6]; Apr42 = Apr[7];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[4];
                        Api21 = Api[1]; Api22 = Api[5];
                        Api31 = Api[2]; Api32 = Api[6];
                        Api41 = Api[3]; Api42 = Api[7];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1];
                        Api21 = Api[2]; Api22 = Api[3];
                        Api31 = Api[4]; Api32 = Api[5];
                        Api41 = Api[6]; Api42 = Api[7];
                    } else if( transa == 'C' ) {
                        Api11 = -Api[0]; Api12 = -Api[1];
                        Api21 = -Api[2]; Api22 = -Api[3];
                        Api31 = -Api[4]; Api32 = -Api[5];
                        Api41 = -Api[6]; Api42 = -Api[7];
                    } else {/* transa == 'G' */
                        Api11 = -Api[0]; Api12 = -Api[4];
                        Api21 = -Api[1]; Api22 = -Api[5];
                        Api31 = -Api[2]; Api32 = -Api[6];
                        Api41 = -Api[3]; Api42 = -Api[7];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24
							   - Api41 * Bpi14 - Api42 * Bpi24;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24
							   - Api31 * Bpi14 - Api32 * Bpi24;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24
							   - Api21 * Bpi14 - Api22 * Bpi24;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24
							   - Api11 * Bpi14 - Api12 * Bpi24;
						Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24
							   + Api41 * Bpr14 + Api42 * Bpr24;
						Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24
							   + Api31 * Bpr14 + Api32 * Bpr24;
						Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24
							   + Api21 * Bpr14 + Api22 * Bpr24;
						Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24
							   + Api11 * Bpr14 + Api12 * Bpr24;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23
							   - Api41 * Bpi13 - Api42 * Bpi23;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23
							   - Api31 * Bpi13 - Api32 * Bpi23;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23
							   - Api21 * Bpi13 - Api22 * Bpi23;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23
							   - Api11 * Bpi13 - Api12 * Bpi23;
						Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23
							   + Api41 * Bpr13 + Api42 * Bpr23;
						Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23
							   + Api31 * Bpr13 + Api32 * Bpr23;
						Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23
							   + Api21 * Bpr13 + Api22 * Bpr23;
						Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23
							   + Api11 * Bpr13 + Api12 * Bpr23;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22
							   - Api41 * Bpi12 - Api42 * Bpi22;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22
							   - Api31 * Bpi12 - Api32 * Bpi22;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22
							   - Api21 * Bpi12 - Api22 * Bpi22;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22
							   - Api11 * Bpi12 - Api12 * Bpi22;
						Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22
							   + Api41 * Bpr12 + Api42 * Bpr22;
						Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22
							   + Api31 * Bpr12 + Api32 * Bpr22;
						Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22
							   + Api21 * Bpr12 + Api22 * Bpr22;
						Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22
							   + Api11 * Bpr12 + Api12 * Bpr22;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21
							   - Api41 * Bpi11 - Api42 * Bpi21;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21
							   - Api31 * Bpi11 - Api32 * Bpi21;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21
							   - Api21 * Bpi11 - Api22 * Bpi21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21
							   - Api11 * Bpi11 - Api12 * Bpi21;
						Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21
							   + Api41 * Bpr11 + Api42 * Bpr21;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21
							   + Api31 * Bpr11 + Api32 * Bpr21;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21
							   + Api21 * Bpr11 + Api22 * Bpr21;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21
							   + Api11 * Bpr11 + Api12 * Bpr21;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[15]= Api41 * Bpr14 + Api42 * Bpr24;
							Cpi[14]= Api31 * Bpr14 + Api32 * Bpr24;
							Cpi[13]= Api21 * Bpr14 + Api22 * Bpr24;
							Cpi[12]= Api11 * Bpr14 + Api12 * Bpr24;
						case 3:
							Cpi[11]= Api41 * Bpr13 + Api42 * Bpr23;
							Cpi[10]= Api31 * Bpr13 + Api32 * Bpr23;
							Cpi[9] = Api21 * Bpr13 + Api22 * Bpr23;
							Cpi[8] = Api11 * Bpr13 + Api12 * Bpr23;
						case 2:
							Cpi[7] = Api41 * Bpr12 + Api42 * Bpr22;
							Cpi[6] = Api31 * Bpr12 + Api32 * Bpr22;
							Cpi[5] = Api21 * Bpr12 + Api22 * Bpr22;
							Cpi[4] = Api11 * Bpr12 + Api12 * Bpr22;
						case 1:
							Cpi[3] = Api41 * Bpr11 + Api42 * Bpr21;
							Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21;
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24;
							Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24;
							Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24;
							Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24;
						case 3:
							Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23;
							Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23;
							Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23;
							Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23;
						case 2:
							Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22;
							Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22;
							Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22;
							Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22;
						case 1:
							Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21;
							Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21;
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21;
					    }
					}
				}
				break;

			case 3: /* (4 x 3)*(3 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[4]; Apr13 = Apr[8];
  				    Apr21 = Apr[1]; Apr22 = Apr[5]; Apr23 = Apr[9];
  				    Apr31 = Apr[2]; Apr32 = Apr[6]; Apr33 = Apr[10];
  				    Apr41 = Apr[3]; Apr42 = Apr[7]; Apr43 = Apr[11];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2];
  				    Apr21 = Apr[3]; Apr22 = Apr[4]; Apr23 = Apr[5];
  				    Apr31 = Apr[6]; Apr32 = Apr[7]; Apr33 = Apr[8];
  				    Apr41 = Apr[9]; Apr42 = Apr[10];Apr43 = Apr[11];
				}
                if( Api ) {
                    if(        transa == 'N' ) {
                        Api11 = Api[0]; Api12 = Api[4]; Api13 = Api[8];
                        Api21 = Api[1]; Api22 = Api[5]; Api23 = Api[9];
                        Api31 = Api[2]; Api32 = Api[6]; Api33 = Api[10];
                        Api41 = Api[3]; Api42 = Api[7]; Api43 = Api[11];
                    } else if( transa == 'T' ) {
                        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2];
                        Api21 = Api[3]; Api22 = Api[4]; Api23 = Api[5];
                        Api31 = Api[6]; Api32 = Api[7]; Api33 = Api[8];
                        Api41 = Api[9]; Api42 = Api[10];Api43 = Api[11];
                    } else if( transa == 'C' ) {
                        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2];
                        Api21 =-Api[3]; Api22 =-Api[4]; Api23 =-Api[5];
                        Api31 =-Api[6]; Api32 =-Api[7]; Api33 =-Api[8];
                        Api41 =-Api[9]; Api42 =-Api[10];Api43 =-Api[11];
                    } else {/* transa == 'G' */
                        Api11 =-Api[0]; Api12 =-Api[4]; Api13 =-Api[8];
                        Api21 =-Api[1]; Api22 =-Api[5]; Api23 =-Api[9];
                        Api31 =-Api[2]; Api32 =-Api[6]; Api33 =-Api[10];
                        Api41 =-Api[3]; Api42 =-Api[7]; Api43 =-Api[11];
                    }
                }
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24 + Apr43 * Bpr34
							   - Api41 * Bpi14 - Api42 * Bpi24 - Api43 * Bpi34;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34
							   - Api31 * Bpi14 - Api32 * Bpi24 - Api33 * Bpi34;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34;
						Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24 + Apr43 * Bpi34
							   + Api41 * Bpr14 + Api42 * Bpr24 + Api43 * Bpr34;
						Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34
							   + Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34;
						Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
						Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23 + Apr43 * Bpr33
							   - Api41 * Bpi13 - Api42 * Bpi23 - Api43 * Bpi33;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33
							   - Api31 * Bpi13 - Api32 * Bpi23 - Api33 * Bpi33;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33;
						Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23 + Apr43 * Bpi33
							   + Api41 * Bpr13 + Api42 * Bpr23 + Api43 * Bpr33;
						Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33
							   + Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33;
						Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
						Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22 + Apr43 * Bpr32
							   - Api41 * Bpi12 - Api42 * Bpi22 - Api43 * Bpi32;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32
							   - Api31 * Bpi12 - Api32 * Bpi22 - Api33 * Bpi32;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32;
						Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22 + Apr43 * Bpi32
							   + Api41 * Bpr12 + Api42 * Bpr22 + Api43 * Bpr32;
						Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32
							   + Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32;
						Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
						Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21 + Apr43 * Bpr31
							   - Api41 * Bpi11 - Api42 * Bpi21 - Api43 * Bpi31;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31
							   - Api31 * Bpi11 - Api32 * Bpi21 - Api33 * Bpi31;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31;
						Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21 + Apr43 * Bpi31
							   + Api41 * Bpr11 + Api42 * Bpr21 + Api43 * Bpr31;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31
							   + Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24 + Apr43 * Bpr34;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23 + Apr43 * Bpr33;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22 + Apr43 * Bpr32;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21 + Apr43 * Bpr31;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31;
					}
					if( Api ) {
					    switch( n ) {
        				case 4:
							Cpi[15]= Api41 * Bpr14 + Api42 * Bpr24 + Api43 * Bpr34;
							Cpi[14]= Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34;
							Cpi[13]= Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34;
							Cpi[12]= Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34;
						case 3:
							Cpi[11]= Api41 * Bpr13 + Api42 * Bpr23 + Api43 * Bpr33;
							Cpi[10]= Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33;
							Cpi[9] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33;
							Cpi[8] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33;
						case 2:
							Cpi[7] = Api41 * Bpr12 + Api42 * Bpr22 + Api43 * Bpr32;
							Cpi[6] = Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32;
							Cpi[5] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32;
							Cpi[4] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32;
						case 1:
							Cpi[3] = Api41 * Bpr11 + Api42 * Bpr21 + Api43 * Bpr31;
							Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31;
							Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31;
							Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        				case 4:
							Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24 + Apr43 * Bpi34;
							Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34;
							Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34;
							Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34;
						case 3:
							Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23 + Apr43 * Bpi33;
							Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33;
							Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33;
							Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33;
						case 2:
							Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22 + Apr43 * Bpi32;
							Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32;
							Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32;
							Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32;
						case 1:
							Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21 + Apr43 * Bpi31;
							Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31;
							Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31;
							Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31;
					    }
					}
				}
				break;

			case 4: /* (4 x 4)*(4 x n) */
				if( transa == 'N' || transa == 'G' ) {
				    Apr11 = Apr[0]; Apr12 = Apr[4]; Apr13 = Apr[8]; Apr14 = Apr[12];
				    Apr21 = Apr[1]; Apr22 = Apr[5]; Apr23 = Apr[9]; Apr24 = Apr[13];
				    Apr31 = Apr[2]; Apr32 = Apr[6]; Apr33 = Apr[10];Apr34 = Apr[14];
				    Apr41 = Apr[3]; Apr42 = Apr[7]; Apr43 = Apr[11];Apr44 = Apr[15];
				} else { /* transa == 'T' || transa == 'C' */
				    Apr11 = Apr[0]; Apr12 = Apr[1]; Apr13 = Apr[2]; Apr14 = Apr[3];
				    Apr21 = Apr[4]; Apr22 = Apr[5]; Apr23 = Apr[6]; Apr24 = Apr[7];
				    Apr31 = Apr[8]; Apr32 = Apr[9]; Apr33 = Apr[10];Apr34 = Apr[11];
				    Apr41 = Apr[12];Apr42 = Apr[13];Apr43 = Apr[14];Apr44 = Apr[15];
				}
				if( Api ) {
				    if( transa == 'N' ) {
				        Api11 = Api[0]; Api12 = Api[4]; Api13 = Api[8]; Api14 = Api[12];
				        Api21 = Api[1]; Api22 = Api[5]; Api23 = Api[9]; Api24 = Api[13];
				        Api31 = Api[2]; Api32 = Api[6]; Api33 = Api[10];Api34 = Api[14];
				        Api41 = Api[3]; Api42 = Api[7]; Api43 = Api[11];Api44 = Api[15];
				    } else if( transa == 'T' ) {
				        Api11 = Api[0]; Api12 = Api[1]; Api13 = Api[2]; Api14 = Api[3];
				        Api21 = Api[4]; Api22 = Api[5]; Api23 = Api[6]; Api24 = Api[7];
				        Api31 = Api[8]; Api32 = Api[9]; Api33 = Api[10];Api34 = Api[11];
				        Api41 = Api[12];Api42 = Api[13];Api43 = Api[14];Api44 = Api[15];
					} else if( transa == 'G' ) {
				        Api11 =-Api[0]; Api12 =-Api[4]; Api13 =-Api[8]; Api14 =-Api[12];
				        Api21 =-Api[1]; Api22 =-Api[5]; Api23 =-Api[9]; Api24 =-Api[13];
				        Api31 =-Api[2]; Api32 =-Api[6]; Api33 =-Api[10];Api34 =-Api[14];
				        Api41 =-Api[3]; Api42 =-Api[7]; Api43 =-Api[11];Api44 =-Api[15];
				    } else { /* transa == 'C' */
				        Api11 =-Api[0]; Api12 =-Api[1]; Api13 =-Api[2]; Api14 =-Api[3];
				        Api21 =-Api[4]; Api22 =-Api[5]; Api23 =-Api[6]; Api24 =-Api[7];
				        Api31 =-Api[8]; Api32 =-Api[9]; Api33 =-Api[10];Api34 =-Api[11];
				        Api41 =-Api[12];Api42 =-Api[13];Api43 =-Api[14];Api44 =-Api[15];
				    }
				}
				if( Api && Bpi ) {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24 + Apr43 * Bpr34 + Apr44 * Bpr44
							   - Api41 * Bpi14 - Api42 * Bpi24 - Api43 * Bpi34 - Api44 * Bpi44;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34 + Apr34 * Bpr44
							   - Api31 * Bpi14 - Api32 * Bpi24 - Api33 * Bpi34 - Api34 * Bpi44;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44
							   - Api21 * Bpi14 - Api22 * Bpi24 - Api23 * Bpi34 - Api24 * Bpi44;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44
							   - Api11 * Bpi14 - Api12 * Bpi24 - Api13 * Bpi34 - Api14 * Bpi44;
						Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24 + Apr43 * Bpi34 + Apr44 * Bpi44
							   + Api41 * Bpr14 + Api42 * Bpr24 + Api43 * Bpr34 + Api44 * Bpr44;
						Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34 + Apr34 * Bpi44
							   + Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34 + Api34 * Bpr44;
						Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44
							   + Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44
							   + Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23 + Apr43 * Bpr33 + Apr44 * Bpr43
							   - Api41 * Bpi13 - Api42 * Bpi23 - Api43 * Bpi33 - Api44 * Bpi43;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33 + Apr34 * Bpr43
							   - Api31 * Bpi13 - Api32 * Bpi23 - Api33 * Bpi33 - Api34 * Bpi43;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43
							   - Api21 * Bpi13 - Api22 * Bpi23 - Api23 * Bpi33 - Api24 * Bpi43;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43
							   - Api11 * Bpi13 - Api12 * Bpi23 - Api13 * Bpi33 - Api14 * Bpi43;
						Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23 + Apr43 * Bpi33 + Apr44 * Bpi43
							   + Api41 * Bpr13 + Api42 * Bpr23 + Api43 * Bpr33 + Api44 * Bpr43;
						Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33 + Apr34 * Bpi43
							   + Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33 + Api34 * Bpr43;
						Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43
							   + Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43
							   + Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22 + Apr43 * Bpr32 + Apr44 * Bpr42
							   - Api41 * Bpi12 - Api42 * Bpi22 - Api43 * Bpi32 - Api44 * Bpi42;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32 + Apr34 * Bpr42
							   - Api31 * Bpi12 - Api32 * Bpi22 - Api33 * Bpi32 - Api34 * Bpi42;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42
							   - Api21 * Bpi12 - Api22 * Bpi22 - Api23 * Bpi32 - Api24 * Bpi42;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42
							   - Api11 * Bpi12 - Api12 * Bpi22 - Api13 * Bpi32 - Api14 * Bpi42;
						Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22 + Apr43 * Bpi32 + Apr44 * Bpi42
							   + Api41 * Bpr12 + Api42 * Bpr22 + Api43 * Bpr32 + Api44 * Bpr42;
						Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32 + Apr34 * Bpi42
							   + Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32 + Api34 * Bpr42;
						Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42
							   + Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42
							   + Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21 + Apr43 * Bpr31 + Apr44 * Bpr41
							   - Api41 * Bpi11 - Api42 * Bpi21 - Api43 * Bpi31 - Api44 * Bpi41;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31 + Apr34 * Bpr41
							   - Api31 * Bpi11 - Api32 * Bpi21 - Api33 * Bpi31 - Api34 * Bpi41;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41
							   - Api21 * Bpi11 - Api22 * Bpi21 - Api23 * Bpi31 - Api24 * Bpi41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41
							   - Api11 * Bpi11 - Api12 * Bpi21 - Api13 * Bpi31 - Api14 * Bpi41;
						Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21 + Apr43 * Bpi31 + Apr44 * Bpi41
							   + Api41 * Bpr11 + Api42 * Bpr21 + Api43 * Bpr31 + Api44 * Bpr41;
						Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31 + Apr34 * Bpi41
							   + Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31 + Api34 * Bpr41;
						Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41
							   + Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41
							   + Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					}
				} else {
					switch( n ) {
        			case 4:
						Cpr[15]= Apr41 * Bpr14 + Apr42 * Bpr24 + Apr43 * Bpr34 + Apr44 * Bpr44;
						Cpr[14]= Apr31 * Bpr14 + Apr32 * Bpr24 + Apr33 * Bpr34 + Apr34 * Bpr44;
						Cpr[13]= Apr21 * Bpr14 + Apr22 * Bpr24 + Apr23 * Bpr34 + Apr24 * Bpr44;
						Cpr[12]= Apr11 * Bpr14 + Apr12 * Bpr24 + Apr13 * Bpr34 + Apr14 * Bpr44;
					case 3:
						Cpr[11]= Apr41 * Bpr13 + Apr42 * Bpr23 + Apr43 * Bpr33 + Apr44 * Bpr43;
						Cpr[10]= Apr31 * Bpr13 + Apr32 * Bpr23 + Apr33 * Bpr33 + Apr34 * Bpr43;
						Cpr[9] = Apr21 * Bpr13 + Apr22 * Bpr23 + Apr23 * Bpr33 + Apr24 * Bpr43;
						Cpr[8] = Apr11 * Bpr13 + Apr12 * Bpr23 + Apr13 * Bpr33 + Apr14 * Bpr43;
					case 2:
						Cpr[7] = Apr41 * Bpr12 + Apr42 * Bpr22 + Apr43 * Bpr32 + Apr44 * Bpr42;
						Cpr[6] = Apr31 * Bpr12 + Apr32 * Bpr22 + Apr33 * Bpr32 + Apr34 * Bpr42;
						Cpr[5] = Apr21 * Bpr12 + Apr22 * Bpr22 + Apr23 * Bpr32 + Apr24 * Bpr42;
						Cpr[4] = Apr11 * Bpr12 + Apr12 * Bpr22 + Apr13 * Bpr32 + Apr14 * Bpr42;
					case 1:
						Cpr[3] = Apr41 * Bpr11 + Apr42 * Bpr21 + Apr43 * Bpr31 + Apr44 * Bpr41;
						Cpr[2] = Apr31 * Bpr11 + Apr32 * Bpr21 + Apr33 * Bpr31 + Apr34 * Bpr41;
						Cpr[1] = Apr21 * Bpr11 + Apr22 * Bpr21 + Apr23 * Bpr31 + Apr24 * Bpr41;
						Cpr[0] = Apr11 * Bpr11 + Apr12 * Bpr21 + Apr13 * Bpr31 + Apr14 * Bpr41;
					}
					if( Api ) {
					    switch( n ) {
        			    case 4:
						    Cpi[15]= Api41 * Bpr14 + Api42 * Bpr24 + Api43 * Bpr34 + Api44 * Bpr44;
						    Cpi[14]= Api31 * Bpr14 + Api32 * Bpr24 + Api33 * Bpr34 + Api34 * Bpr44;
						    Cpi[13]= Api21 * Bpr14 + Api22 * Bpr24 + Api23 * Bpr34 + Api24 * Bpr44;
						    Cpi[12]= Api11 * Bpr14 + Api12 * Bpr24 + Api13 * Bpr34 + Api14 * Bpr44;
					    case 3:
						    Cpi[11]= Api41 * Bpr13 + Api42 * Bpr23 + Api43 * Bpr33 + Api44 * Bpr43;
						    Cpi[10]= Api31 * Bpr13 + Api32 * Bpr23 + Api33 * Bpr33 + Api34 * Bpr43;
						    Cpi[9] = Api21 * Bpr13 + Api22 * Bpr23 + Api23 * Bpr33 + Api24 * Bpr43;
						    Cpi[8] = Api11 * Bpr13 + Api12 * Bpr23 + Api13 * Bpr33 + Api14 * Bpr43;
					    case 2:
						    Cpi[7] = Api41 * Bpr12 + Api42 * Bpr22 + Api43 * Bpr32 + Api44 * Bpr42;
						    Cpi[6] = Api31 * Bpr12 + Api32 * Bpr22 + Api33 * Bpr32 + Api34 * Bpr42;
						    Cpi[5] = Api21 * Bpr12 + Api22 * Bpr22 + Api23 * Bpr32 + Api24 * Bpr42;
						    Cpi[4] = Api11 * Bpr12 + Api12 * Bpr22 + Api13 * Bpr32 + Api14 * Bpr42;
					    case 1:
						    Cpi[3] = Api41 * Bpr11 + Api42 * Bpr21 + Api43 * Bpr31 + Api44 * Bpr41;
						    Cpi[2] = Api31 * Bpr11 + Api32 * Bpr21 + Api33 * Bpr31 + Api34 * Bpr41;
						    Cpi[1] = Api21 * Bpr11 + Api22 * Bpr21 + Api23 * Bpr31 + Api24 * Bpr41;
						    Cpi[0] = Api11 * Bpr11 + Api12 * Bpr21 + Api13 * Bpr31 + Api14 * Bpr41;
					    }
					}
					if( Bpi ) {
					    switch( n ) {
        			    case 4:
						    Cpi[15]= Apr41 * Bpi14 + Apr42 * Bpi24 + Apr43 * Bpi34 + Apr44 * Bpi44;
						    Cpi[14]= Apr31 * Bpi14 + Apr32 * Bpi24 + Apr33 * Bpi34 + Apr34 * Bpi44;
						    Cpi[13]= Apr21 * Bpi14 + Apr22 * Bpi24 + Apr23 * Bpi34 + Apr24 * Bpi44;
						    Cpi[12]= Apr11 * Bpi14 + Apr12 * Bpi24 + Apr13 * Bpi34 + Apr14 * Bpi44;
					    case 3:
						    Cpi[11]= Apr41 * Bpi13 + Apr42 * Bpi23 + Apr43 * Bpi33 + Apr44 * Bpi43;
						    Cpi[10]= Apr31 * Bpi13 + Apr32 * Bpi23 + Apr33 * Bpi33 + Apr34 * Bpi43;
						    Cpi[9] = Apr21 * Bpi13 + Apr22 * Bpi23 + Apr23 * Bpi33 + Apr24 * Bpi43;
						    Cpi[8] = Apr11 * Bpi13 + Apr12 * Bpi23 + Apr13 * Bpi33 + Apr14 * Bpi43;
					    case 2:
						    Cpi[7] = Apr41 * Bpi12 + Apr42 * Bpi22 + Apr43 * Bpi32 + Apr44 * Bpi42;
						    Cpi[6] = Apr31 * Bpi12 + Apr32 * Bpi22 + Apr33 * Bpi32 + Apr34 * Bpi42;
						    Cpi[5] = Apr21 * Bpi12 + Apr22 * Bpi22 + Apr23 * Bpi32 + Apr24 * Bpi42;
						    Cpi[4] = Apr11 * Bpi12 + Apr12 * Bpi22 + Apr13 * Bpi32 + Apr14 * Bpi42;
					    case 1:
						    Cpi[3] = Apr41 * Bpi11 + Apr42 * Bpi21 + Apr43 * Bpi31 + Apr44 * Bpi41;
						    Cpi[2] = Apr31 * Bpi11 + Apr32 * Bpi21 + Apr33 * Bpi31 + Apr34 * Bpi41;
						    Cpi[1] = Apr21 * Bpi11 + Apr22 * Bpi21 + Apr23 * Bpi31 + Apr24 * Bpi41;
						    Cpi[0] = Apr11 * Bpi11 + Apr12 * Bpi21 + Apr13 * Bpi31 + Apr14 * Bpi41;
					    }
					}
				}
				break;
			}
			break;
		}

/*----------------------------------------------------------------------------
 * Vector dot product (1 x K) * (K x 1)
 *---------------------------------------------------------------------------- */

    } else if( m == 1 && n == 1 ) {
        z = RealKindDotProduct(k, Apr, Api, ai, Bpr, Bpi, bi, dot_method);
        *Cpr = z.r;
        if( Cpi ) {
            *Cpi = z.i;
        }

/*----------------------------------------------------------------------------
 * Vector outer product (M x 1) * (1 x N)
 *---------------------------------------------------------------------------- */

    } else if( k == 1 ) {
        RealKindOuterProduct(m, n, Apr, Api, transa, Bpr, Bpi, transb, Cpr, Cpi, outer_method);
      
/*----------------------------------------------------------------------------
 * Matrix times vector (M x K) * (K x 1)
 *---------------------------------------------------------------------------- */

    } else if( n == 1 ) {

/*----------------------------------------------------------------------------
 * If the first matrix is not transposed, use calls to xGEMV. Also use this
 * method if running in the 'MATLAB' or 'BLAS' mode, or the number of processors
 * is greater than 2.
 *---------------------------------------------------------------------------- */

        if( transa == 'N' || transa == 'G' || mtimesx_mode == MTIMESX_BLAS ||
			mtimesx_mode == MTIMESX_MATLAB || omp_get_num_procs() > 2 ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xGEMV) "\n");
			}
            if( transa == 'G' ) ptransa = 'N';
            xGEMV(PTRANSA, M1, N1, ONE, Apr, LDA, Bpr, INCX, ZERO, Cpr, INCY);
            if( Bpi ) {
                xGEMV(PTRANSA, M1, N1, BI, Apr, LDA, Bpi, INCX, ZERO, Cpi, INCY);
                if( Api ) {                               /* (complex matrix) * (complex vector) */
                    xGEMV(PTRANSA, M1, N1, AIBI, Api, LDA, Bpi, INCX,  ONE, Cpr, INCY);
                    xGEMV(PTRANSA, M1, N1, AI, Api, LDA, Bpr, INCX,  ONE, Cpi, INCY);
                } else {                                  /* (real matrix) * (complex vector)
                     * already done */
                }
            } else {
                if( Api ) {                               /* (complex matrix) * (real vector) */
                    xGEMV(PTRANSA, M1, N1, AI, Api, LDA, Bpr, INCX, ZERO, Cpi, INCY);
                } else {                                  /* (real matrix) * (real vector)
                     * already done */
                }
            }

/* Alternate method ... doesn't match MATLAB exactly ... not up to date
 *
 *         if( transa == 'N' || transa == 'G' || matlab ) {
 *             if( transa == 'G' ) ptransa = 'N';
 *             xGEMV(PTRANSA, M1, N1, ONE, Apr, LDA, Bpr, INCX, ZERO, Cpr, INCY);
 *             if( mxIsComplex(A) ) {
 *                 xGEMV(PTRANSA, M1, N1, AI, Api, LDA, Bpr, INCX, ZERO, Cpi, INCY);
 *                 if( mxIsComplex(B) ) {                    // (complex matrix) * (complex vector)
 *                     xGEMV(PTRANSA, M1, N1, AIBI, Api, LDA, Bpi, INCX,  ONE, Cpr, INCY);
 *                     xGEMV(PTRANSA, M1, N1, BI, Apr, LDA, Bpi, INCX,  ONE, Cpi, INCY);
 *                 } else {                                  // (complex matrix) * (real vector)
 *                     // already done
 *                 }
 *             } else {
 *                 if( mxIsComplex(B) ) {                    // (real matrix) * (complex vector)
 *                     xGEMV(PTRANSA, M1, N1, BI, Apr, LDA, Bpi, INCX, ZERO, Cpi, INCY);
 *                 } else {                                  // (real matrix) * (real vector)
 *                     // already done
 *                 }
 *             } */

/*-----------------------------------------------------------------------------------------
 * Else if the first matrix is transposed, then use calls to xDOT instead (faster) because
 * the matrix can be accessed as a series of contiguous column vectors.
 *----------------------------------------------------------------------------------------- */

        } else { /* transa == 'T' || transa == 'C' */
            apr = Apr;
            api = Api;
            if( Api ) {
                for( i=0; i<m; i++ ) {                   /* (complex matrix) * (vector) */
                    z = RealKindDotProduct(k, apr, api, ai, Bpr, Bpi, bi, dot_method);
                    Cpr[i] = z.r;
                    Cpi[i] = z.i;
                    apr += k;
                    api += k;
                }
            } else {                                     /* (real matrix) * (complex vector) */
                if( Bpi ) {
                    for( i=0; i<m; i++ ) {
                        z = RealKindDotProduct(k, apr, Api, ai, Bpr, Bpi, bi, dot_method);
                        Cpr[i] = z.r;
                        Cpi[i] = z.i;
                        apr += k;
                    }
                } else {                                 /* (real matrix) * (real vector) */
                    for( i=0; i<m; i++ ) {
                        z = RealKindDotProduct(k, apr, Api, ai, Bpr, Bpi, bi, dot_method);
                        Cpr[i] = z.r;
                        apr += k;
                    }
                }
            }
        }

/*----------------------------------------------------------------------------------------
 * Vector times matrix (1 x K) * (K x N)
 *---------------------------------------------------------------------------------------- */

    } else if( m == 1 ) {

/*----------------------------------------------------------------------------------------
 * If the second matrix is transposed, then use calls to xGEMV with the arguments reversed.
 * Also use this method if running in 'MATLAB' or 'BLAS' mode, or the number of processors
 * is greater than 2.
 *---------------------------------------------------------------------------------------- */

        if( transb == 'C' || transb == 'T' || mtimesx_mode == MTIMESX_BLAS ||
			mtimesx_mode == MTIMESX_MATLAB || omp_get_num_procs() > 2 ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xGEMV) "\n");
			}
            if( transb == 'C' || transb == 'T' ) {
                ptransb = 'N';
            } else {
                ptransb = 'T';
            }
            xGEMV(PTRANSB, M2, N2, ONE, Bpr, LDB, Apr, INCX, ZERO, Cpr, INCY);
            if( Api ) {
                xGEMV(PTRANSB, M2, N2, AI, Bpr, LDB, Api, INCX, ZERO, Cpi, INCY);
                if( Bpi ) {                               /* (complex matrix) * (complex vector) */
                    xGEMV(PTRANSB, M2, N2, AIBI, Bpi, LDB, Api, INCX,  ONE, Cpr, INCY);
                    xGEMV(PTRANSB, M2, N2, BI, Bpi, LDB, Apr, INCX,  ONE, Cpi, INCY);
                } else {                                  /* (complex matrix) * (real vector)
                     * already done */
                }
            } else {
                if( Bpi ) {                               /* (real matrix) * (complex vector) */
                    xGEMV(PTRANSB, M2, N2, BI, Bpi, LDB, Apr, INCX, ZERO, Cpi, INCY);
                } else {                                  /* (real matrix) * (real vector)
                     * already done */
                }
            }

/* Alternate method ... doesn't match MATLAB exactly ... not up to date
 *
 *         if( transb == 'C' || transb == 'T' || matlab ) {
 *             if( transb == 'C' || transb == 'T' ) {
 *                 ptransb = 'N';
 *             } else {
 *                 ptransb = 'T';
 *             }
 *             xGEMV(PTRANSB, M2, N2, ONE, Bpr, LDB, Apr, INCX, ZERO, Cpr, INCY);
 *             if( mxIsComplex(B) ) {
 *                 xGEMV(PTRANSB, M2, N2, BI, Bpi, LDB, Apr, INCX, ZERO, Cpi, INCY);
 *                 if( mxIsComplex(A) ) {                    // (complex matrix) * (complex vector)
 *                     xGEMV(PTRANSB, M2, N2, AIBI, Bpi, LDB, Api, INCX,  ONE, Cpr, INCY);
 *                     xGEMV(PTRANSB, M2, N2, AI, Bpr, LDB, Api, INCX,  ONE, Cpi, INCY);
 *                 } else {                                  // (real matrix) * (complex vector)
 *                     // already done
 *                 }
 *             } else {
 *                 if( mxIsComplex(A) ) {                    // (complex matrix) * (real vector)
 *                     xGEMV(PTRANSB, M2, N2, AI, Bpr, LDB, Api, INCX, ZERO, Cpi, INCY);
 *                 } else {                                  // (real matrix) * (real vector)
 *                     // already done
 *                 }
 *             } */

/*-----------------------------------------------------------------------------------------
 * Else if the second matrix is not transposed, then use calls to dot product instead
 * (faster) because the matrix can be accessed as a series of contiguous column vectors.
 *----------------------------------------------------------------------------------------- */

        } else {
            bpr = Bpr;
            bpi = Bpi;
            if( Bpi ) {
                for( i=0; i<n; i++ ) {
                    z = RealKindDotProduct(k, Apr, Api, ai, bpr, bpi, bi, dot_method);
                    Cpr[i] = z.r;
                    Cpi[i] = z.i;
                    bpr += k;
                    bpi += k;
                }
            } else {                                     /* (complex vector) * (real matrix) */
                if( Api ) {
                    for( i=0; i<n; i++ ) {
                        z = RealKindDotProduct(k, Apr, Api, ai, bpr, Bpi, bi, dot_method);
                        Cpr[i] = z.r;
                        Cpi[i] = z.i;
                        bpr += k;
                    }
                } else {                                 /* (real vector) * (real matrix) */
                    for( i=0; i<n; i++ ) {
                        z = RealKindDotProduct(k, Apr, Api, ai, bpr, Bpi, bi, dot_method);
                        Cpr[i] = z.r;
                        bpr += k;
                    }
                }
            }
        }
        
/*---------------------------------------------------------------------------------
 * Matrix times matrix (M x K) * (K x N) with N small and first matrix transposed.
 * Use dot product (faster) because the 1st matrix can be accessed as a series of
 * contiguous column vectors. When the column size reaches about 8 then the memory
 * access efficiency of the BLAS routines increases and this custom method is no
 * longer faster. The number 8 is likely machine / implementation dependent. Only
 * use this method if running in one of the 'SPEED' or 'LOOPS' type modes.
 *--------------------------------------------------------------------------------- */

    } else if( !(mtimesx_mode == MTIMESX_BLAS || mtimesx_mode == MTIMESX_MATLAB) &&
		       n < 7 && (transa == 'T' || transa == 'C') && (transb == 'N' || transb == 'G') ) {
        bpr = Bpr;
        bpi = Bpi;
        cpr = Cpr;
        cpi = Cpi;
        for( j=0; j<n; j++ ) {
            apr = Apr;
            api = Api;
            for( i=0; i<m; i++ ) {
                z = RealKindDotProduct(k, apr, api, ai, bpr, bpi, bi, dot_method);
                *cpr++ = z.r;
                if( cpi ) *cpi++ = z.i;
                apr += k;
                if( api ) api += k;
            }
            bpr += k;
            if( bpi ) bpi += k;
        }

/*---------------------------------------------------------------------------------------------
 * Matrix times matrix (M x K) *(K x N)
 *--------------------------------------------------------------------------------------------- */

    } else {

/*/--------------------------------------------------------------------------------------------
 * If Matrix product is actually the same matrix, use calls to the symmetric routines xSYRK and
 * xSYR2K where possible. These only work on the lower or upper triangle part of the matrix, so
 * we will have to fill out the other half manually, but even so this will be faster. Some of
 * these will not match MATLAB exactly, so only run them in cases where matching is not required.
 *--------------------------------------------------------------------------------------------- */

        if( Apr == Bpr && ((transa == 'N' && transb == 'T') ||
                           (transa == 'T' && transb == 'N')) ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYRK) " and " TOKENSTRING(xSYR2K) "\n");
			}
            xSYRK(UPLO, TRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( Api ) {
                xSYRK(UPLO, TRANSA, N, K, MINUSONE, Api, LDA, ONE, Cpr, LDC);
                xSYR2K(UPLO,TRANSA, N, K, ONE, Api, LDA, Apr, LDA, ZERO, Cpi, LDC);
                xFILLPOS(Cpi, n);
            }
            xFILLPOS(Cpr, n);

        } else if( Apr == Bpr && (!(mtimesx_mode == MTIMESX_BLAS || mtimesx_mode == MTIMESX_MATLAB) ||
			                 (!Api && !Bpi)) &&
                             ((transa == 'G' && transb == 'C') ||
                              (transa == 'C' && transb == 'G')) ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYRK) " and " TOKENSTRING(xSYR2K) "\n");
			}
            if( transa == 'G')  ptransa = 'N';
            xSYRK(UPLO, PTRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( Api ) {
                xSYRK(UPLO, PTRANSA, N, K, MINUSONE, Api, LDA, ONE, Cpr, LDC);
                xSYR2K(UPLO,PTRANSA, N, K, MINUSONE, Apr, LDA, Api, LDA, ZERO, Cpi, LDC);
                xFILLPOS(Cpi, n);
            }
            xFILLPOS(Cpr, n);

        } else if( Apr == Bpr && ((transa == 'N' && transb == 'C') ||
                                  (transa == 'T' && transb == 'G' && 
							  (!(mtimesx_mode == MTIMESX_BLAS || mtimesx_mode == MTIMESX_MATLAB) ||
							  (!Api && !Bpi)))) ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYRK) " and " TOKENSTRING(xGEMM) "\n");
			}
            if( transb == 'G' ) ptransb = 'N';
            xSYRK(UPLO, TRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( Api ) {
                xSYRK(UPLO, TRANSA, N, K, ONE, Api, LDA, ONE, Cpr, LDC);
                xGEMM(TRANSA, PTRANSB, M, N, K, ONE, Api, LDA, Apr, LDA, ZERO, Cpi, LDC);
                xFILLNEG(Cpi, n);
            }
            xFILLPOS(Cpr, n);

        } else if( Apr == Bpr && ((transa == 'C' && transb == 'N') ||
                                  (transa == 'G' && transb == 'T' &&
							  (!(mtimesx_mode == MTIMESX_BLAS || mtimesx_mode == MTIMESX_MATLAB) ||
							  (!Api && !Bpi)))) ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYRK) " and " TOKENSTRING(xGEMM) "\n");
			}
            if( transa == 'G' ) ptransa = 'N';
            xSYRK(UPLO, PTRANSA, N, K, ONE, Apr, LDA, ZERO, Cpr, LDC);
            if( Api ) {
                xSYRK(UPLO, PTRANSA, N, K, ONE, Api, LDA, ONE, Cpr, LDC);
                xGEMM(PTRANSA, TRANSB, M, N, K, ONE, Apr, LDA, Api, LDA, ZERO, Cpi, LDC);
                xFILLNEG(Cpi, n);
            }
            xFILLPOS(Cpr, n);

/*-------------------------------------------------------------------------------------------
 * Else this is not a symmetric case, so just call the general matrix multiply routine xGEMM.
 *------------------------------------------------------------------------------------------- */

        } else {
			if( debug_message ) {
				mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xGEMM) "\n");
			}
            if( transa == 'G' ) ptransa = 'N';
            if( transb == 'G' ) ptransb = 'N';
            xGEMM(PTRANSA, PTRANSB, M, N, K, ONE, Apr, LDA, Bpr, LDB, ZERO, Cpr, LDC);
            if( Bpi ) {
                xGEMM(PTRANSA, PTRANSB, M, N, K, BI, Apr, LDA, Bpi, LDB, ZERO, Cpi, LDC);
                if( Api ) {                               /* (complex matrix) * (complex matrix) */
                    xGEMM(PTRANSA, PTRANSB, M, N, K, AIBI, Api, LDA, Bpi, LDB,  ONE, Cpr, LDC);
                    xGEMM(PTRANSA, PTRANSB, M, N, K, AI, Api, LDA, Bpr, LDB,  ONE, Cpi, LDC);
                } else {                                  /* (real matrix) * (complex matrix)
                     * already done */
                }
            } else {
                if( Api ) {                               /* (complex matrix) * (real matrix) */
                    xGEMM(PTRANSA, PTRANSB, M, N, K, AI, Api, LDA, Bpr, LDB, ZERO, Cpi, LDC);
                } else {                                  /* (real matrix) * (real matrix)
                     * already done */
                }
            }

/* Alternate method ... doesn't match MATLAB exactly ... not up to date
 *
 *         } else {
 *             if( transa == 'G' ) ptransa = 'N';
 *             if( transb == 'G' ) ptransb = 'N';
 *             xGEMM(PTRANSA, PTRANSB, M, N, K, ONE, Apr, LDA, Bpr, LDB, ZERO, Cpr, LDC);
 *             if( mxIsComplex(A) ) {
 *                 xGEMM(PTRANSA, PTRANSB, M, N, K, AI, Api, LDA, Bpr, LDB, ZERO, Cpi, LDC);
 *                 if( mxIsComplex(B) ) {                    // (complex matrix) * (complex matrix)
 *                     xGEMM(PTRANSA, PTRANSB, M, N, K, AIBI, Api, LDA, Bpi, LDB,  ONE, Cpr, LDC);
 *                     xGEMM(PTRANSA, PTRANSB, M, N, K, BI, Apr, LDA, Bpi, LDB,  ONE, Cpi, LDC);
 *                 } else {                                  // (complex matrix) * (real matrix)
 *                     // already done
 *                 }
 *             } else {
 *                 if( mxIsComplex(B) ) {                    // (real matrix) * (complex matrix)
 *                     xGEMM(PTRANSA, PTRANSB, M, N, K, BI, Apr, LDA, Bpi, LDB, ZERO, Cpi, LDC);
 *                 } else {                                  // (real matrix) * (real matrix)
 *                     // already done
 *                 }
 *             } */

        }
    }

/*----------------------------------------------------------------------------
 * End Outer Loop to process all of the individual matrix multiplies. Increment
 * the matrix pointers to point to the next pair of matrices to be multiplied.
 *---------------------------------------------------------------------------- */

    if( ip < p-1 ) {
        j = 2;
        while( ++Cindx[j] == Cdims[j] ) Cindx[j++] = 0;
        Ap = Bp = 0;
        Ablock = Asize;
        Bblock = Bsize;
        for( j=2; j<Cndim; j++ ) {
            if( Cindx[j] < Adimz[j] ) Ap += Cindx[j] * Ablock;
            Ablock *= Adimz[j];
            if( Cindx[j] < Bdimz[j] ) Bp += Cindx[j] * Bblock;
            Bblock *= Bdimz[j];
        }
        Apr = Apr0 + Ap;
        if( Api ) Api = Api0 + Ap;
        Bpr = Bpr0 + Bp;
        if( Bpi ) Bpi = Bpi0 + Bp;
        Cpr += Csize;
        if( Cpi ) Cpi += Csize;
    }
	debug_message = 0;
    }

    mxFree(Cindx);
    mxFree(Cdims);
    mxFree(Adimz);
    mxFree(Bdimz);

/*---------------------------------------------------------------------------------
 * If the imaginary part is all zero, then free it and set the pointer to NULL.
 *--------------------------------------------------------------------------------- */

    Cpi = mxGetImagData(C);
    if( AllRealZero(Cpi, m*n*p) ) {
        mxFree(Cpi);
        mxSetImagData(C, NULL);
    }

/*---------------------------------------------------------------------------------
 * Done.
 *--------------------------------------------------------------------------------- */

	if( destroyA ) mxDestroyArray(A);
	if( destroyB ) mxDestroyArray(B);
    return result;

}

/*--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *-------------------------------------------------------------------------------------- */

/*--------------------------------------------------------------------------------------
 * Returns 1 (true) if all of the elements are zero. Returns 0 (false) if at least one
 * of the elements is non-zero, or if the pointer to the data is NULL.
 *-------------------------------------------------------------------------------------- */

int AllRealZero(RealKind *x, mwSignedIndex n)
{
    register mwSignedIndex i;
    if( x == NULL ) return 0;
    for(i=0; i<n; i++) {
        if( x[i] != zero ) return 0;
    }
    return 1;
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 0*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 0*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 0*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 0*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = -(*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] - Bpi[i];
            Cpi[i] = Bpr[i] + Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i];
            Cpi[i] = Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + Bpi[i];
            Cpi[i] = Bpr[i] - Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr - *bpi;
                    *Cpi++ = *bpr + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr;
                    *Cpi++ = *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + *bpi;
                    *Cpi++ = *bpr - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 - 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + Bpi[i];
            Cpi[i] = Bpi[i] - Bpr[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 - 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  Bpr[i] - Bpi[i];
            Cpi[i] = -Bpi[i] - Bpr[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 - 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr;
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + *bpi;
                    *Cpi++ = *bpi - *bpr;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 - 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr;
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   *bpr - *bpi;
                    *Cpi++ = - *bpr - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + ai*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */
void RealKindEqP1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =      Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] - ai * Bpi[i];
            Cpi[i] = ai * Bpr[i] + Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + ai*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =      Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = Bpr[i] + ai * Bpi[i];
            Cpi[i] = ai * Bpr[i] - Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + ai*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =       *bpr;
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr - ai * (*bpi);
                    *Cpi++ = ai * (*bpr) + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (1 + ai*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqP1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =       *bpr;
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = *bpr + ai * (*bpi);
                    *Cpi++ = ai * (*bpr) - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - Bpi[i];
            Cpi[i] =  Bpr[i] - Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + Bpi[i];
            Cpi[i] =  Bpr[i] + Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - (*bpi);
                    *Cpi++ =   *bpr  - (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + (*bpi);
                    *Cpi++ =   *bpr  + (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 - 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1M1TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + Bpi[i];
            Cpi[i] = -Bpi[i] - Bpr[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 - 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1M1TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - Bpi[i];
            Cpi[i] =  Bpi[i] - Bpr[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 - 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1M1TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + (*bpi);
                    *Cpi++ = -(*bpr) - (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 - 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */
void RealKindEqM1M1TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - (*bpi);
                    *Cpi++ = -(*bpr) + (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 0*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P0TimesRealKindN(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] = -Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 0*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P0TimesRealKindG(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i];
            Cpi[i] =  Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 0*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P0TimesRealKindT(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ = -(*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + 0*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1P0TimesRealKindC(RealKind *Cpr, RealKind *Cpi,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr);
                    *Cpi++ =   *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + ai*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1PxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =    - Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] - ai * Bpi[i];
            Cpi[i] =  ai * Bpr[i] - Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + ai*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1PxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] =    - Bpr[i];
            Cpi[i] = ai * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = -Bpr[i] + ai * Bpi[i];
            Cpi[i] =  ai * Bpr[i] + Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + ai*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1PxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =    - (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) - ai * (*bpi);
                    *Cpi++ = ai * (*bpr) - *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (-1 + ai*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqM1PxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =    - (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = -(*bpr) + ai * (*bpi);
                    *Cpi++ = ai * (*bpr) + *bpi;
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =      Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i] - Bpi[i];
            Cpi[i] = Bpr[i] + ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =      Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i] + Bpi[i];
            Cpi[i] = Bpr[i] - ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =       *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr) - (*bpi);
                    *Cpi++ = (*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =       *bpr;
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr) + (*bpi);
                    *Cpi++ = (*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar - 1*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxM1TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =    - Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  ar * Bpr[i] + Bpi[i];
            Cpi[i] = -Bpr[i] + ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar - 1*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxM1TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] =    - Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =  ar * Bpr[i] - Bpi[i];
            Cpi[i] = -Bpr[i] - ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar - 1*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxM1TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =    - (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) + (*bpi);
                    *Cpi++ = -(*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar - 1*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxM1TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ =    - (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) - (*bpi);
                    *Cpi++ = -(*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 0*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP0TimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
            Cpi[i] = ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 0*i) * (Bpr - Bpi * i) */
void RealKindEqPxP0TimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        for(i=0; i<n; i++) {
            Cpr[i] = ar * Bpr[i];
        }
    } else {
        for(i=0; i<n; i++) {
            Cpr[i] =   ar * Bpr[i];
            Cpi[i] =  -ar * Bpi[i];
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 0*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP0TimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr);
                    *Cpi++ =  ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + 0*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxP0TimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =   ar * (*bpr);
                    *Cpi++ =  -ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + ai*i) * (Bpr + Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxPxTimesRealKindN(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
                Cpi[i] = ai * Bpr[i];
            }
        }
    } else {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i];
                Cpi[i] =  ar * Bpi[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i] - ai * Bpi[i];
                Cpi[i] =  ai * Bpr[i] + ar * Bpi[i];
            }
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + ai*i) * (Bpr - Bpi * i)
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxPxTimesRealKindG(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSignedIndex n)
{
    register mwSignedIndex i;

    if( Bpi == NULL ) {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] = ar * Bpr[i];
                Cpi[i] = ai * Bpr[i];
            }
        }
    } else {
        if( ai == zero ) {
            for(i=0; i<n; i++) {
                Cpr[i] =   ar * Bpr[i];
                Cpi[i] =  -ar * Bpi[i];
            }
        } else {
            for(i=0; i<n; i++) {
                Cpr[i] =  ar * Bpr[i] + ai * Bpi[i];
                Cpi[i] =  ai * Bpr[i] - ar * Bpi[i];
            }
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + ai*i) * (Bpr + Bpi * i)T
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxPxTimesRealKindT(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) - ai * (*bpi);
                    *Cpi++ =  ai * (*bpr) + ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * C = (ar + ai*i) * (Bpr + Bpi * i)C
 *-------------------------------------------------------------------------------------- */

void RealKindEqPxPxTimesRealKindC(RealKind *Cpr, RealKind *Cpi, RealKind ar, RealKind ai,
                                  RealKind *Bpr, RealKind *Bpi, mwSize m2, mwSize n2, mwSignedIndex p)
{
    register mwSize i, j;
    RealKind *bpr, *Br = Bpr;
    RealKind *bpi, *Bi = Bpi;
    register mwSignedIndex ip;

    if( Bpi == NULL ) {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ = ar * (*bpr);
                    *Cpi++ = ai * (*bpr);
                    bpr += m2;
                }
                Bpr++;
            }
            Br += m2 * n2;
        }
    } else {
        for( ip=0; ip<p; ip++ ) {
            Bpr = Br;
            Bpi = Bi;
            for( i=0; i<m2; i++ ) {
                bpr = Bpr;
                bpi = Bpi;
                for( j=0; j<n2; j++ ) {
                    *Cpr++ =  ar * (*bpr) + ai * (*bpi);
                    *Cpi++ =  ai * (*bpr) - ar * (*bpi);
                    bpr += m2;
                    bpi += m2;
                }
                Bpr++;
                Bpi++;
            }
            Br += m2 * n2;
            Bi += m2 * n2;
        }
    }
}

/*--------------------------------------------------------------------------------------
 * Fill the upper triangle with contents of the lower triangle
 *-------------------------------------------------------------------------------------- */

void xFILLPOS(RealKind *Cpr, mwSignedIndex n)
{
    RealKind *source, *target;
    register mwSignedIndex i, j;

    source = Cpr + 1;
    target = Cpr + n;
    for( i=1; i<n; i++ ) {
        for( j=i; j<n; j++ ) {
            *target = *source;
            target += n;
            source++;
        }
        source += i + 1;
        target = source + n - 1;
    }
}

/*--------------------------------------------------------------------------------------
 * Add the -Cpr transpose to the current array
 *-------------------------------------------------------------------------------------- */

void xFILLNEG(RealKind *Cpr, mwSignedIndex n)
{
    RealKind *source, *target;
    register mwSignedIndex i, j;

    source = Cpr;
    target = Cpr;
    for( i=0; i<n; i++ ) {
        for( j=i; j<n; j++ ) {
            *target -= *source;
            if( i != j ) {
                *source = -(*target);
            }
            target += n;
            source++;
        }
        source += i + 1;
        target = source;
    }
}

/*------------------------------------------------------------------------------------------
 * Dot Product type calculation. For conjugate cases, where ai == -1 or bi == -1, use custom
 * loops instead of doing the actual multply by ai or bi. Also, use loop unrolling to speed
 * up the calculations and improve the accuracy. For PC WinXP, the balance between speed and
 * accuracy seemed to be optimal at a blocksize of 10. If the MATLAB mode is set, then just
 * duplicate the BLAS calls that MATLAB uses (slower since it accesses the variables twice).
 *------------------------------------------------------------------------------------------ */

#ifdef _OPENMP

/*------------------------------------------------------------------------------------------
 * Interface for OpenMP enabled compiling.
 *------------------------------------------------------------------------------------------ */

struct RealKindComplex RealKindDotProduct(mwSignedIndex k,
                                          RealKind *Apr, RealKind *Api, RealKind ai,
                                          RealKind *Bpr, RealKind *Bpi, RealKind bi,
										  int dot_method)
{
    struct RealKindComplex RealKindZZ[OMP_DOT_ARRAY_SIZE];
	struct RealKindComplex z;
	struct RealKindComplex *zz;
	int i, j, n;

	if( dot_method != METHOD_LOOPS_OMP ||
	    (mtimesx_mode == MTIMESX_SPEED_OMP && k < OMP_DOT_SMALL) ) {
		if( debug_message ) {
            if( dot_method == METHOD_BLAS ) {
			    mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xDOT) "\n");
			} else {
		        mexPrintf("MTIMESX: LOOPS dot product\n");
			}
			debug_message = 0;
		}
		z = RealKindDotProductX(k, Apr, Api, ai, Bpr, Bpi, bi, dot_method);
	} else {
		if( debug_message ) {
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS dot product\n");
			debug_message = 0;
		}
		if( max_threads <= OMP_DOT_ARRAY_SIZE ) {
			zz = RealKindZZ;
		} else {
			zz = mxMalloc(max_threads * sizeof(*zz));
		}
        omp_set_dynamic(1);
        #pragma omp parallel num_threads(max_threads)
		{
			RealKind *Apr_, *Api_, *Bpr_, *Bpi_;
			mwSize k_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = k / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				k_ = k - offset;
			} else {
				k_ = blocksize;
			}
			Apr_ = Apr + offset;
			if( Api ) {
				Api_ = Api + offset;
			} else {
				Api_ = NULL;
			}
			Bpr_ = Bpr + offset;
			if( Bpi ) {
				Bpi_ = Bpi + offset;
			} else {
				Bpi_ = NULL;
			}
			zz[thread_num] = RealKindDotProductX(k_, Apr_, Api_, ai, Bpr_, Bpi_, bi, dot_method);
		}

/* Combine thread results in a pre-determined order, so result with same inputs will be the */
/* same regardless of the order in which the threads execute and complete. Use a binary     */
/* reduction scheme to maintain accuracy.                                                   */

		n = threads_used - 1;
		while( n ) {
			if( !(n & 1) ) {
				zz[n-1].r += zz[n].r;
				zz[n-1].i += zz[n].i;
				n--;
			}
			for( i=0, j=0; i<n; i+=2, j++ ) {
				zz[j].r = zz[i].r + zz[i+1].r;
				zz[j].i = zz[i].i + zz[i+1].i;
			}
			n >>= 1;
		}
		z = zz[0];
/*
    	z.r = zero;
		z.i = zero;
		for( i=0; i<threads_used; i++ ) {
			z.r += zz[i].r;
			z.i += zz[i].i;
		}

*/
		if( zz != RealKindZZ ) {
			mxFree(zz);
		}
	}
    return z;
}

struct RealKindComplex RealKindDotProductX(mwSignedIndex k,
                                          RealKind *Apr, RealKind *Api, RealKind ai,
                                          RealKind *Bpr, RealKind *Bpi, RealKind bi,
										  int dot_method)

#else

/*------------------------------------------------------------------------------------------
 * Interface for non-OpenMP enabled compiling.
 *------------------------------------------------------------------------------------------ */

struct RealKindComplex RealKindDotProduct(mwSignedIndex k,
                                          RealKind *Apr, RealKind *Api, RealKind ai,
                                          RealKind *Bpr, RealKind *Bpi, RealKind bi,
										  int dot_method)

#endif

{
    register double sr = 0.0, si = 0.0;
    register mwSignedIndex i;
    struct RealKindComplex z;
    mwSignedIndex k10, inc = 1;

    if( dot_method == METHOD_BLAS ) {
#ifndef _OPENMP
		if( debug_message ) {
		    mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xDOT) "\n");
			debug_message = 0;
		}
#endif
        sr = xDOT( &k, Apr, &inc, Bpr, &inc );
        if( Api != NULL ) {
            if( Api != Bpi || ai == bi ) {
                si = xDOT( &k, Api, &inc, Bpr, &inc ) * ai;
            }
            if( Bpi != NULL ) {
                sr -= xDOT( &k, Api, &inc, Bpi, &inc ) * ai * bi;
                if( Api != Bpi ) {
                    si += xDOT( &k, Apr, &inc, Bpi, &inc ) * bi;
                } else if( ai == bi ){
                    si += si;
                }
            }
        } else if( Bpi != NULL ) {
            si = xDOT( &k, Apr, &inc, Bpi, &inc ) * bi;
        }
        z.r = (RealKind) sr;
        z.i = (RealKind) si;
        return z;
    }

#ifndef _OPENMP
	if( debug_message ) {
        mexPrintf("MTIMESX: LOOPS dot product\n");
		debug_message = 0;
	}
#endif
    k10 = k % 10;
    if( Api != NULL ) {
        if( Bpi != NULL ) {                    /* (complex vector) dot (complex vector) */
            if( ai == one ) {
                if( bi == one ) {
                    if( Apr == Bpr ) {

                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Apr[i] - Api[i] * Api[i];
                        si += Api[i] * Apr[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Apr[i  ] - Api[i  ] * Api[i  ]
                           +  Apr[i+1] * Apr[i+1] - Api[i+1] * Api[i+1]
                           +  Apr[i+2] * Apr[i+2] - Api[i+2] * Api[i+2]
                           +  Apr[i+3] * Apr[i+3] - Api[i+3] * Api[i+3]
                           +  Apr[i+4] * Apr[i+4] - Api[i+4] * Api[i+4]
                           +  Apr[i+5] * Apr[i+5] - Api[i+5] * Api[i+5]
                           +  Apr[i+6] * Apr[i+6] - Api[i+6] * Api[i+6]
                           +  Apr[i+7] * Apr[i+7] - Api[i+7] * Api[i+7]
                           +  Apr[i+8] * Apr[i+8] - Api[i+8] * Api[i+8]
                           +  Apr[i+9] * Apr[i+9] - Api[i+9] * Api[i+9];
                        si += Api[i  ] * Apr[i  ]
                           +  Api[i+1] * Apr[i+1]
                           +  Api[i+2] * Apr[i+2]
                           +  Api[i+3] * Apr[i+3]
                           +  Api[i+4] * Apr[i+4]
                           +  Api[i+5] * Apr[i+5]
                           +  Api[i+6] * Apr[i+6]
                           +  Api[i+7] * Apr[i+7]
                           +  Api[i+8] * Apr[i+8]
                           +  Api[i+9] * Apr[i+9];
                    }
                    si += si;
                        
                    } else {
                    
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] - Api[i] * Bpi[i];
                        si += Api[i] * Bpr[i] + Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] - Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] - Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] - Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] - Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] - Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] - Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] - Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] - Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] - Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] - Api[i+9] * Bpi[i+9];
                        si += Api[i  ] * Bpr[i  ] + Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] + Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] + Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] + Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] + Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] + Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] + Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] + Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] + Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] + Apr[i+9] * Bpi[i+9];
                    }
                    
                    }
                } else {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] + Api[i] * Bpi[i];
                        si += Api[i] * Bpr[i] - Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] + Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] + Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] + Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] + Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] + Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] + Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] + Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] + Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] + Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] + Api[i+9] * Bpi[i+9];
                        si += Api[i  ] * Bpr[i  ] - Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] - Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] - Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] - Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] - Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] - Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] - Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] - Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] - Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] - Apr[i+9] * Bpi[i+9];
                    }
                }
            } else {
                if( bi == one ) {
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] + Api[i] * Bpi[i];
                        si += Apr[i] * Bpi[i] - Api[i] * Bpr[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] + Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] + Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] + Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] + Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] + Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] + Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] + Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] + Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] + Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] + Api[i+9] * Bpi[i+9];
                        si += Apr[i  ] * Bpi[i  ] - Api[i  ] * Bpr[i  ]
                           +  Apr[i+1] * Bpi[i+1] - Api[i+1] * Bpr[i+1]
                           +  Apr[i+2] * Bpi[i+2] - Api[i+2] * Bpr[i+2]
                           +  Apr[i+3] * Bpi[i+3] - Api[i+3] * Bpr[i+3]
                           +  Apr[i+4] * Bpi[i+4] - Api[i+4] * Bpr[i+4]
                           +  Apr[i+5] * Bpi[i+5] - Api[i+5] * Bpr[i+5]
                           +  Apr[i+6] * Bpi[i+6] - Api[i+6] * Bpr[i+6]
                           +  Apr[i+7] * Bpi[i+7] - Api[i+7] * Bpr[i+7]
                           +  Apr[i+8] * Bpi[i+8] - Api[i+8] * Bpr[i+8]
                           +  Apr[i+9] * Bpi[i+9] - Api[i+9] * Bpr[i+9];
                    }
                } else {
                    
                    if( Apr == Bpr ) {
                        
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Apr[i] - Api[i] * Api[i];
                        si -= Api[i] * Apr[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Apr[i  ] - Api[i  ] * Api[i  ]
                           +  Apr[i+1] * Apr[i+1] - Api[i+1] * Api[i+1]
                           +  Apr[i+2] * Apr[i+2] - Api[i+2] * Api[i+2]
                           +  Apr[i+3] * Apr[i+3] - Api[i+3] * Api[i+3]
                           +  Apr[i+4] * Apr[i+4] - Api[i+4] * Api[i+4]
                           +  Apr[i+5] * Apr[i+5] - Api[i+5] * Api[i+5]
                           +  Apr[i+6] * Apr[i+6] - Api[i+6] * Api[i+6]
                           +  Apr[i+7] * Apr[i+7] - Api[i+7] * Api[i+7]
                           +  Apr[i+8] * Apr[i+8] - Api[i+8] * Api[i+8]
                           +  Apr[i+9] * Apr[i+9] - Api[i+9] * Api[i+9];
                        si -= Api[i  ] * Apr[i  ]
                           +  Api[i+1] * Apr[i+1]
                           +  Api[i+2] * Apr[i+2]
                           +  Api[i+3] * Apr[i+3]
                           +  Api[i+4] * Apr[i+4]
                           +  Api[i+5] * Apr[i+5]
                           +  Api[i+6] * Apr[i+6]
                           +  Api[i+7] * Apr[i+7]
                           +  Api[i+8] * Apr[i+8]
                           +  Api[i+9] * Apr[i+9];
                    }
                    si += si;
                    
                    } else {
                    
                    for( i=0; i<k10; i++ ) {
                        sr += Apr[i] * Bpr[i] - Api[i] * Bpi[i];
                        si -= Api[i] * Bpr[i] + Apr[i] * Bpi[i];
                    }
                    for( i=k10; i<k; i+=10 ) {
                        sr += Apr[i  ] * Bpr[i  ] - Api[i  ] * Bpi[i  ]
                           +  Apr[i+1] * Bpr[i+1] - Api[i+1] * Bpi[i+1]
                           +  Apr[i+2] * Bpr[i+2] - Api[i+2] * Bpi[i+2]
                           +  Apr[i+3] * Bpr[i+3] - Api[i+3] * Bpi[i+3]
                           +  Apr[i+4] * Bpr[i+4] - Api[i+4] * Bpi[i+4]
                           +  Apr[i+5] * Bpr[i+5] - Api[i+5] * Bpi[i+5]
                           +  Apr[i+6] * Bpr[i+6] - Api[i+6] * Bpi[i+6]
                           +  Apr[i+7] * Bpr[i+7] - Api[i+7] * Bpi[i+7]
                           +  Apr[i+8] * Bpr[i+8] - Api[i+8] * Bpi[i+8]
                           +  Apr[i+9] * Bpr[i+9] - Api[i+9] * Bpi[i+9];
                        si -= Api[i  ] * Bpr[i  ] + Apr[i  ] * Bpi[i  ]
                           +  Api[i+1] * Bpr[i+1] + Apr[i+1] * Bpi[i+1]
                           +  Api[i+2] * Bpr[i+2] + Apr[i+2] * Bpi[i+2]
                           +  Api[i+3] * Bpr[i+3] + Apr[i+3] * Bpi[i+3]
                           +  Api[i+4] * Bpr[i+4] + Apr[i+4] * Bpi[i+4]
                           +  Api[i+5] * Bpr[i+5] + Apr[i+5] * Bpi[i+5]
                           +  Api[i+6] * Bpr[i+6] + Apr[i+6] * Bpi[i+6]
                           +  Api[i+7] * Bpr[i+7] + Apr[i+7] * Bpi[i+7]
                           +  Api[i+8] * Bpr[i+8] + Apr[i+8] * Bpi[i+8]
                           +  Api[i+9] * Bpr[i+9] + Apr[i+9] * Bpi[i+9];
                    }
                    
                    }
                }
            }
            z.i = (RealKind) si;
        } else {                                  /* (complex vector) dot (real vector) */
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
                si += Api[i] * Bpr[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
                si += Api[i  ] * Bpr[i  ]
                   +  Api[i+1] * Bpr[i+1]
                   +  Api[i+2] * Bpr[i+2]
                   +  Api[i+3] * Bpr[i+3]
                   +  Api[i+4] * Bpr[i+4]
                   +  Api[i+5] * Bpr[i+5]
                   +  Api[i+6] * Bpr[i+6]
                   +  Api[i+7] * Bpr[i+7]
                   +  Api[i+8] * Bpr[i+8]
                   +  Api[i+9] * Bpr[i+9];
            }
            z.i = (RealKind) (si * ai);
        }
    } else {
        if( Bpi != NULL ) {                    /* (real vector) dot (complex vector) */
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
                si += Apr[i] * Bpi[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
                si += Apr[i  ] * Bpi[i  ]
                   +  Apr[i+1] * Bpi[i+1]
                   +  Apr[i+2] * Bpi[i+2]
                   +  Apr[i+3] * Bpi[i+3]
                   +  Apr[i+4] * Bpi[i+4]
                   +  Apr[i+5] * Bpi[i+5]
                   +  Apr[i+6] * Bpi[i+6]
                   +  Apr[i+7] * Bpi[i+7]
                   +  Apr[i+8] * Bpi[i+8]
                   +  Apr[i+9] * Bpi[i+9];
            }
            z.i = (RealKind) (si * bi);
        } else {                                  /* (real vector) dot (real vector) */
            for( i=0; i<k10; i++ ) {
                sr += Apr[i] * Bpr[i];
            }
            for( i=k10; i<k; i+=10 ) {
                sr += Apr[i  ] * Bpr[i  ]
                   +  Apr[i+1] * Bpr[i+1]
                   +  Apr[i+2] * Bpr[i+2]
                   +  Apr[i+3] * Bpr[i+3]
                   +  Apr[i+4] * Bpr[i+4]
                   +  Apr[i+5] * Bpr[i+5]
                   +  Apr[i+6] * Bpr[i+6]
                   +  Apr[i+7] * Bpr[i+7]
                   +  Apr[i+8] * Bpr[i+8]
                   +  Apr[i+9] * Bpr[i+9];
            }
        }
    }
    z.r = (RealKind) sr;
    return z;
}

/*----------------------------------------------------------------------------------------
 * Outer Product calculation. Use custom loops for all of the special cases involving
 * conjugates and transposes to minimize the total number of operations involved.
 *---------------------------------------------------------------------------------------- */

#ifdef _OPENMP

/*------------------------------------------------------------------------------------------
 * Interface for OpenMP enabled compiling.
 *------------------------------------------------------------------------------------------ */

void RealKindOuterProduct(mwSignedIndex m, mwSignedIndex n,
                          RealKind *Apr, RealKind *Api, char transa,
                          RealKind *Bpr, RealKind *Bpi, char transb,
                          RealKind *Cpr, RealKind *Cpi, int outer_method)
{
	int t;
	RealKind ai, bi;

	t = max_threads;
	if( n < t ) t = n;
	if( outer_method != METHOD_LOOPS_OMP ||
	    (mtimesx_mode == MTIMESX_SPEED_OMP && m*n < OMP_OUTER_SMALL) ) {
		if( debug_message ) {
            if( outer_method == METHOD_BLAS ) {
				if( debug_message ) {
		            if( Apr != Bpr || 1 ) {  /* Force calls to generic xGER instead of xSYR */
		                mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xGER) "\n");
		            } else {
			            if( Api ) {
    		                ai = ( transa == 'G' || transa == 'C' ) ? -one : one;
		                    bi = ( transb == 'G' || transb == 'C' ) ? -one : one;
				            if( ai == bi ) {
		                        mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) " and " TOKENSTRING(xSYR2) "\n");
				            } else {
		                        mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) " and " TOKENSTRING(xGER) "\n");
				            }
  			            } else {
		                    mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) "\n");
			            }
					}
					debug_message = 0;
		        }
			} else {
		        mexPrintf("MTIMESX: LOOPS outer product\n");
			}
			debug_message = 0;
		}
		RealKindOuterProductX(m, n, Apr, Api, transa, Bpr, Bpi, transb, Cpr, Cpi, outer_method);
	} else {
		if( debug_message ) {
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS outer product\n");
			debug_message = 0;
		}
        omp_set_dynamic(1);
        #pragma omp parallel num_threads(t)
		{
			RealKind *Bpr_, *Bpi_, *Cpr_, *Cpi_;
			mwSize n_, blocksize, offset;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			blocksize = n / num_threads;
			offset = thread_num * blocksize;
			if( thread_num == num_threads-1 ) {
				n_ = n - offset;
			} else {
				n_ = blocksize;
			}
			Bpr_ = Bpr + offset;
			if( Bpi ) {
				Bpi_ = Bpi + offset;
			} else {
				Bpi_ = NULL;
			}
			Cpr_ = Cpr + offset * m;
			if( Cpi ) {
				Cpi_ = Cpi + offset * m;
			} else {
				Cpi_ = NULL;
			}
			RealKindOuterProductX(m, n_, Apr, Api, transa, Bpr_, Bpi_, transb, Cpr_, Cpi_, outer_method);
		}
	}
}

void RealKindOuterProductX(mwSignedIndex m, mwSignedIndex n,
                          RealKind *Apr, RealKind *Api, char transa,
                          RealKind *Bpr, RealKind *Bpi, char transb,
                          RealKind *Cpr, RealKind *Cpi, int outer_method)

#else

/*------------------------------------------------------------------------------------------
 * Interface for non-OpenMP enabled compiling.
 *------------------------------------------------------------------------------------------ */

void RealKindOuterProduct(mwSignedIndex m, mwSignedIndex n,
                          RealKind *Apr, RealKind *Api, char transa,
                          RealKind *Bpr, RealKind *Bpi, char transb,
                          RealKind *Cpr, RealKind *Cpi, int outer_method)

#endif

{
    register mwSignedIndex i, j;
    mwSignedIndex kk;
	mwSignedIndex inc = 1;
	RealKind One = one;
	RealKind ai, bi, aibi;
	char uplo = 'L';

	if( outer_method == METHOD_BLAS ) {
		ai = ( transa == 'G' || transa == 'C' ) ? -one : one;
		bi = ( transb == 'G' || transb == 'C' ) ? -one : one;
		aibi = - ai * bi;
		if( Apr != Bpr || 1 ) {  /* Force calls to generic xGER instead of xSYR */
#ifndef _OPENMP
		    if( debug_message ) {
		        mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xGER) "\n");
			    debug_message = 0;
		    }
#endif
		    xGER( M, N, ONE, Apr, INCX, Bpr, INCY, Cpr, M );
            if( Api ) {
        	    xGER( M, N, AI, Api, INCX, Bpr, INCY, Cpi, M );
                if( Bpi ) {
        		    xGER( M, N, AIBI, Api, INCX, Bpi, INCY, Cpr, M );
        		    xGER( M, N, BI, Apr, INCX, Bpi, INCY, Cpi, M );
			    }
            } else if( Bpi ) {
        	    xGER( M, N, BI, Apr, INCX, Bpi, INCY, Cpi, M );
            }
		} else {
		    xSYR( UPLO, M, ONE, Apr, INCX, Cpr, M );
			if( Api ) {
		        xSYR( UPLO, M, AIBI, Api, INCX, Cpr, M );
				if( ai == bi ) {
#ifndef _OPENMP
		            if( debug_message ) {
		                mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) " and " TOKENSTRING(xSYR2) "\n");
			            debug_message = 0;
		            }
#endif
		            xSYR2( UPLO, M, AI, Apr, INCX, Api, INCY, Cpi, M );
                    xFILLPOS(Cpi, m);
				} else {
#ifndef _OPENMP
		            if( debug_message ) {
		                mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) " and " TOKENSTRING(xGER) "\n");
			            debug_message = 0;
		            }
#endif
        		    xGER( M, N, BI, Apr, INCX, Api, INCY, Cpi, M );
                    xFILLNEG(Cpi, m);
				}
#ifndef _OPENMP
			} else {
		        if( debug_message ) {
		            mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xSYR) "\n");
			        debug_message = 0;
		        }
#endif
			}
            xFILLPOS(Cpr, m);
		}
		return;
	}

#ifndef _OPENMP
	if( debug_message ) {
        mexPrintf("MTIMESX: LOOPS outer product\n");
		debug_message = 0;
	}
#endif
	kk = 0;
    if( Api ) {
        if( Bpi ) {
            if( (transa == 'C' || transa == 'G') && (transb == 'C' || transb == 'G') ) {
                for( j=0; j<n; j++ ) { /* (ar + bi*i)(C or G) * (br + bi*i)(C or G) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j] - Api[i] * Bpi[j];
                        Cpi[kk] = - Apr[i] * Bpi[j] - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else if( (transa == 'C' || transa == 'G') && (transb == 'N' || transb == 'T') ) {
                for( j=0; j<n; j++ ) { /* (ar + bi*i)(C or G) * (br + bi*i)(N or T) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] + Api[i] * Bpi[j];
                        Cpi[kk] = Apr[i] * Bpi[j] - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else if( (transa == 'N' || transa == 'T') && (transb == 'C' || transb == 'G') ) {
                for( j=0; j<n; j++ ) { /* (ar + bi*i)(N or T) * (br + bi*i)(C or G) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] + Api[i] * Bpi[j];
                        Cpi[kk] = Api[i] * Bpr[j] - Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            } else {                   /* (ar + bi*i)(N or T) * (br + bi*i)(N or T) */
                for( j=0; j<n; j++ ) {
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j] - Api[i] * Bpi[j];
                        Cpi[kk] = Apr[i] * Bpi[j] + Api[i] * Bpr[j];
                        kk++;
                    }
                }
            }
        } else {
            if( transa == 'C' || transa == 'G' ) {
                for( j=0; j<n; j++ ) { /* (ar + bi*i)(C or G) * (br)(N or T or C or G) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j];
                        Cpi[kk] = - Api[i] * Bpr[j];
                        kk++;
                    }
                }
            } else {
                for( j=0; j<n; j++ ) { /* (ar + bi*i)(N or T) * (br)(N or T or C or G) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j];
                        Cpi[kk] = Api[i] * Bpr[j];
                        kk++;
                    }
                }
            }
        }
    } else {
        if( Bpi ) {
            if( transb == 'C' || transb == 'G' ) {
                for( j=0; j<n; j++ ) { /* (ar)(N or T or C or G) * (br + bi*i)(C or G) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] =   Apr[i] * Bpr[j];
                        Cpi[kk] = - Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            } else {
                for( j=0; j<n; j++ ) { /* (ar)(N or T or C or G) * (br + bi*i)(N or T) */
                    for( i=0; i<m; i++ ) {
                        Cpr[kk] = Apr[i] * Bpr[j];
                        Cpi[kk] = Apr[i] * Bpi[j];
                        kk++;
                    }
                }
            }
        } else {
            for( j=0; j<n; j++ ) { /* (ar)(N or T or C or G) * (br)(N or T or C or G) */
                for( i=0; i<m; i++ ) {
                    Cpr[kk++] = Apr[i] * Bpr[j];
                }
            }
        }
    }
}

/*--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *-------------------------------------------------------------------------------------- */

mxArray *RealScalarTimesReal(mxArray *A, char transa, mwSize m1, mwSize n1,
                             mxArray *B, char transb, mwSize m2, mwSize n2)
{
    mwSize nzmax, Bndim, Cndim, Cp, p;
    mwSize *Bdims, *Cdims;
    mwSignedIndex mn, n, k;
    mxArray *C, *result, *Bt = NULL;
    RealKind *Apr, *Api, *Bpr, *Bpi, *Cpr, *Cpi;
    RealKind ar, ai, br, bi, sr, si;
    char trans;
    mxComplexity complexflag;
    mwIndex *jc;
	int scalar_method;
/*----- */
    if( m2 == 1 && n2 == 1 ) {  /* Make sure scalar is A and array is B */
        if( m1 == 1 && n1 == 1 ) {  /* Check for scalar * scalar */
			if( debug_message ) {
				mexPrintf("MTIMESX: Inline (scalar * scalar) code\n");
				debug_message = 0;
			}
            Apr = mxGetData(A);
            ar = *Apr;
            Api = mxGetImagData(A);
            ai = Api ? *Api : zero;
            if( transa == 'C' || transa == 'G' ) {
                ai = -ai;
            }
            Bpr = mxGetData(B);
            br = *Bpr;
            Bpi = mxGetImagData(B);
            bi = Bpi ? *Bpi : zero;
            if( transb == 'C' || transb == 'G' ) {
                bi = -bi;
            }
            sr = (RealKind) (((double)ar) * ((double)br) - ((double)ai) * ((double)bi));
            si = (RealKind) (((double)ar) * ((double)bi) + ((double)ai) * ((double)br));
            complexflag = (si == zero) ? mxREAL : mxCOMPLEX;
            if( mxIsSparse(A) || mxIsSparse(B) ) {
                result = mxCreateSparse(1, 1, 1, complexflag);
                *mxGetIr(result) = 0;
                jc = mxGetJc(result);
                jc[0] = 0;
                jc[1] = 1;
            } else {
                result = mxCreateNumericMatrix(1, 1, MatlabReturnType, complexflag);
            }
            Cpr = mxGetData(result);
            *Cpr = (RealKind) sr;
            if( complexflag == mxCOMPLEX ) {
                Cpi = mxGetImagData(result);
                *Cpi = (RealKind) si;
            }
            return result;
        } else {
            C = A;
            A = B;
            B = C;
            trans  = transa;
            transa = transb;
            transb = trans;
        }
    }

/*--------------------------------------------------------------------------------------
 * Check for multiplying by 1 and no actual transpose or conjugate is involved. In this
 * case just return a shared data copy, which is very fast since no data copying is done.
 * Also, if there *is* a transpose but the first two dimensions are a vector and there
 * is no actual conjugate involved, we can return a shared data copy with a slight
 * modification of the dimensions ... just switch the first two.
 *-------------------------------------------------------------------------------------- */

    Apr = mxGetData(A);
    Api = mxGetImagData(A);
    ar = *Apr;
    ai = Api ? ((transa == 'C' || transa == 'G') ? -(*Api) :*Api) : zero;

    Bpr = mxGetData(B);
    Bpi = mxGetImagData(B);
    Bndim = mxGetNumberOfDimensions(B);
    Bdims = (mwSize *) mxGetDimensions(B);

    if( ar == one && ai == zero ) {
        if( transb == 'N' || (transb == 'G' && Bpi == NULL) ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: Shared data copy\n");
			}
            result = mxCreateSharedDataCopy(B);
            return result;
        } else if( (Bdims[0] == 1 || Bdims[1] == 1) && (transb == 'T' || (transb == 'C' && Bpi == NULL)) ) {
						if( debug_message ) {
				mexPrintf("MTIMESX: Shared data copy\n");
			}
            result = mxCreateSharedDataCopy(B);
            Cdims = mxMalloc( Bndim * sizeof(*Cdims) );
            for( Cp=2; Cp<Bndim; Cp++) {
                Cdims[Cp] = Bdims[Cp];
            }
            Cdims[0] = Bdims[1];
            Cdims[1] = Bdims[0];
            mxSetDimensions(result,Cdims,Bndim);
            mxFree(Cdims);
            return result;
        }
    }

/*--------------------------------------------------------------------------------------
 * For sparse matrix, do the transpose now and then do the scalar multiply. That way if
 * B is real and A is complex we are only doing the transpose work on one array.
 *-------------------------------------------------------------------------------------- */

    if( (transb == 'T' || transb == 'C') && mxIsSparse(B) ) {
		if( debug_message ) {
			mexPrintf("MTIMESX: Callback to MATLAB to transpose a sparse matrix\n");
		}
        mexCallMATLAB(1, &Bt, 1, &B, "transpose");
        B = Bt;
        transb = (transb == 'T') ? 'N' : 'G';
    }

    if( mxIsSparse(B) ) {
        jc = mxGetJc(B);
        n2 = mxGetN(B);
        n = jc[n2];
    } else {
        n = mxGetNumberOfElements(B);
    }
    if( n < 0 ) {
        mexErrMsgTxt("Number of elements too large ... overflows a signed integer");
    }

    complexflag = (n > 0 && (mxIsComplex(A) || mxIsComplex(B))) ? mxCOMPLEX : mxREAL;

/*-------------------------------------------------------------------------------
 * Construct the dimensions of the result. Also use the p variable to keep track
 * of the total number of individual matrix multiples that are involved. The
 * first two dimensions are simply the result of a single matrix multiply, with
 * accouting for the transb pre-operation. The remaining dimensions are simply
 * copied from B.
 *------------------------------------------------------------------------------- */

    Cndim = Bndim;
    Cdims = mxMalloc( Cndim * sizeof(*Cdims) );
    if( transb == 'N' || transb == 'G' ) {
        Cdims[0] = Bdims[0];
        Cdims[1] = Bdims[1];
    } else {
        Cdims[0] = Bdims[1];
        Cdims[1] = Bdims[0];
    }
    p = 1;
    for( Cp=2; Cp<Cndim; Cp++) {
        p *= (Cdims[Cp] = Bdims[Cp]);
    }

/*------------------------------------------------------------------------------
 * Create output array
 *------------------------------------------------------------------------------ */

    if( mxIsSparse(B) ) {
        result = mxCreateSparse(Cdims[0], Cdims[1], n, complexflag);
        memcpy(mxGetIr(result), mxGetIr(B), n * sizeof(mwIndex));
        memcpy(mxGetJc(result), mxGetJc(B), (n2+1) * sizeof(mwIndex));
    } else if( mxGetNumberOfElements(B) == 0 ) {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, mxREAL);
        mxFree(Cdims);
        return result;
    } else {
        result = mxCreateNumericArray(Cndim, Cdims, MatlabReturnType, complexflag);
    }
    mxFree(Cdims);
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);

    m2 = Bdims[0];
    n2 = Bdims[1];
    mn = m2 * n2;

    if( n == 0 ) {  /* If result is empty, just return right now */
        return result;
    }
    C = result;
    Cpr = mxGetData(C);
    Cpi = mxGetImagData(C);
    Bpr = mxGetData(B);
    Bpi = mxGetImagData(B);

/*------------------------------------------------------------------------------------------------
 * Check for matlab flag. If set, then use a BLAS call for this function.
 *------------------------------------------------------------------------------------------------ */

/*    if( matlab ) {
         * Future upgrade ... insert BLAS calls here and return?
 *    } */

/*------------------------------------------------------------------------------------------------
 * If the matrix is really a vector, then there is no need to do an actual transpose. We can
 * simply strip the transpose operation away and rely on the dimensions to do the transpose.
 *------------------------------------------------------------------------------------------------ */

    if( m2 ==1 || n2 == 1 ) {
        if( transb == 'T' ) transb = 'N';
        if( transb == 'C' ) transb = 'G';
    }

/*----------------------------------------------------------------------------
 * Get scalar product method to use.
 *---------------------------------------------------------------------------- */

	switch( mtimesx_mode )
	{
	case MTIMESX_BLAS:
		scalar_method = METHOD_BLAS;
		break;
	case MTIMESX_LOOPS:
		scalar_method = METHOD_LOOPS;
		break;
	case MTIMESX_LOOPS_OMP:
		scalar_method = METHOD_LOOPS_OMP;
		break;
	case MTIMESX_SPEED_OMP:
		if( max_threads > 1 ) {
		    scalar_method = METHOD_LOOPS_OMP;
			break;
		}
	case MTIMESX_MATLAB:
	case MTIMESX_SPEED:
		if( ai != zero && Bpi ) {
			scalar_method = METHOD_LOOPS;
		} else {
			scalar_method = METHOD_BLAS;
		}
		break;
	}

/*------------------------------------------------------------------------------------------------
 * Do the scalar multiply.
 *------------------------------------------------------------------------------------------------ */

    RealTimesScalar(Cpr, Cpi, Bpr, Bpi, transb, m2, n2, ar, ai, n, p, scalar_method);

/*-------------------------------------------------------------------------------------------
 * If the imaginary part is all zero, then free it and set the pointer to NULL.
 *------------------------------------------------------------------------------------------- */

    if( AllRealZero(Cpi, n) ) {
        mxFree(Cpi);
        Cpi = NULL;
        mxSetImagData(C, NULL);
    }

/*-------------------------------------------------------------------------------------------
 * Clean up sparse matrix and realloc if appropriate.
 *------------------------------------------------------------------------------------------- */

    if( mxIsSparse(C) ) {
        nzmax = mxGetNzmax(C);
        k = spclean(C);
        if( k == 0 ) k = 1;
        if( nzmax - k > REALLOCTOL ) {
			if( debug_message ) {
				mexPrintf("MTIMESX: Reallocate sparse matrix\n");
			}
            mxSetPr(C, myRealloc(Cpr, k * sizeof(RealKind)));
            mxSetIr(C, myRealloc(mxGetIr(C), k * sizeof(mwIndex)));
            if( Cpi != NULL ) {
                mxSetPi(C, myRealloc(Cpi, k * sizeof(RealKind)));
            }
            mxSetNzmax(C, k);
        }
    }
    if( Bt != NULL ) {
        mxDestroyArray(Bt);
    }

    return result;
}

/*--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *--------------------------------------------------------------------------------------
 *-------------------------------------------------------------------------------------- */

/*-------------------------------------------------------------------------------------------
 * Scalar times array slice.
 *------------------------------------------------------------------------------------------- */

/*-------------------------------------------------------------------------------------------
 * OpenMP interface.
 *------------------------------------------------------------------------------------------- */

#ifdef _OPENMP

void RealTimesScalar(RealKind *Cpr, RealKind *Cpi, RealKind *Bpr, RealKind *Bpi, char transb,
                     mwSize m2, mwSize n2, RealKind ar, RealKind ai, mwSize n, mwSize p,
					 int scalar_method)
{
	mwSize q;

	if( transb == 'N' || transb == 'G' ) {
		q = n;
	} else {
		q = m2 * n2;
	}

/*------------------------------------------------------------------------------------------------
 * For small sizes, don't bother with the OpenMP overhead unless forced.
 *------------------------------------------------------------------------------------------------ */

	if( scalar_method != METHOD_LOOPS_OMP ||
	    (mtimesx_mode == MTIMESX_SPEED_OMP && q < OMP_SCALAR_SMALL) ) {
		if( debug_message ) {
	        if( scalar_method == METHOD_BLAS && (transb == 'N' || transb == 'G') ) {
			    mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xAXPY) "\n");
			} else {
			    mexPrintf("MTIMESX: LOOPS scalar multiply\n");
			}
			debug_message = 0;
		}
	    RealTimesScalarX(Cpr, Cpi, Bpr, Bpi, transb, m2, n2, ar, ai, n, p, scalar_method);
	} else {
		if( debug_message ) {
			mexPrintf("MTIMESX: OpenMP multi-threaded LOOPS scalar multiply\n");
			debug_message = 0;
		}
        omp_set_dynamic(1);
        #pragma omp parallel num_threads(max_threads)
		{
			RealKind *Cpr_, *Cpi_, *Bpr_, *Bpi_;
			mwSize n_, p_, blocksize, offset, f;
			int thread_num = omp_get_thread_num();
			int num_threads = omp_get_num_threads();
            #pragma omp master
			{
				threads_used = num_threads;
			}
			if( transb == 'N' || transb == 'G' ) {
			    blocksize = n / num_threads;
			    offset = thread_num * blocksize;
			    if( thread_num == num_threads-1 ) {
				    n_ = n - offset;
			    } else {
				    n_ = blocksize;
			    }
				f = offset;
			} else {
			    blocksize = p / num_threads;
			    offset = thread_num * blocksize;
			    if( thread_num == num_threads-1 ) {
				    p_ = p - offset;
			    } else {
				    p_ = blocksize;
			    }
				f = offset * m2 * n2;
			}
			Cpr_ = Cpr + f;
			if( Cpi ) {
			    Cpi_ = Cpi + f;
			} else {
			    Cpi_ = NULL;
			}
			Bpr_ = Bpr + f;
			if( Bpi ) {
			    Bpi_ = Bpi + f;
			} else {
			    Bpi_ = NULL;
			}
			RealTimesScalarX(Cpr_, Cpi_, Bpr_, Bpi_, transb, m2, n2, ar, ai, n_, p_, scalar_method);
		}
	}
}

void RealTimesScalarX(RealKind *Cpr, RealKind *Cpi, RealKind *Bpr, RealKind *Bpi, char transb,
                      mwSize m2, mwSize n2, RealKind ar, RealKind ai, mwSize n, mwSize p,
					  int scalar_method)

#else

/*-------------------------------------------------------------------------------------------
 * Non-OpenMP interface.
 *------------------------------------------------------------------------------------------- */

void RealTimesScalar(RealKind *Cpr, RealKind *Cpi, RealKind *Bpr, RealKind *Bpi, char transb,
                     mwSize m2, mwSize n2, RealKind ar, RealKind ai, mwSize n, mwSize p,
					 int scalar_method)

#endif

{
	mwSignedIndex inc = 1, k = n;

/*------------------------------------------------------------------------------------------------
 * If scalar multiply mode is BLAS mode and there is no actual transpose involved then do the
 * BLAS calls now.
 *------------------------------------------------------------------------------------------------ */

	if( scalar_method == METHOD_BLAS && (transb == 'N' || transb == 'G') ) {
#ifndef _OPENMP
		if( debug_message ) {
			mexPrintf("MTIMESX: BLAS calls to " TOKENSTRING(xAXPY) "\n");
			debug_message = 0;
		}
#endif
		xAXPY( K, AR, Bpr, INCX, Cpr, INCY );
		if( Cpi ) {
			if( ai != zero ) {
				xAXPY( K, AI, Bpr, INCX, Cpi, INCY );
				if( Bpi ) {
					if( transb == 'N' ) ai = -ai;
					xAXPY( K, AI, Bpi, INCX, Cpr, INCY );
				}
			}
			if( Bpi ) {
				if( transb == 'G' ) ar = -ar;
				xAXPY( K, AR, Bpi, INCX, Cpi, INCY );
			}
		}
		return;
	}

/*------------------------------------------------------------------------------------------------
 * Some specialized cases, no need to multiply by +1 or -1, we can program that directly into the
 * calculations without a multiply. We do need to multiply by zero, however, so that the sign of
 * any -0 that might be present gets carried over into the result, and also any inf or NaN that is
 * present gets a proper result. So no special code for multiplying by zero.
 *------------------------------------------------------------------------------------------------ */

#ifndef _OPENMP
	if( debug_message ) {
		mexPrintf("MTIMESX: LOOPS scalar multiply\n");
		debug_message = 0;
	}
#endif
    if( ar == one ) {
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqP1P1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 + 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqP1P1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 + 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqP1P1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 + 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqP1P1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 + 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqP1M1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 - 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqP1M1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 - 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqP1M1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 - 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqP1M1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 - 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {  /* this case never reached ... it is the shared data copy above */
                RealKindEqP1P0TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 + 0*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqP1P0TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (1 + 0*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqP1P0TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 + 0*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqP1P0TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (1 + 0*i) * (Bpr + Bpi * i)C */
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqP1PxTimesRealKindN(Cpr, Cpi, ai, Bpr, Bpi, n); /* C = (1 + ai*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqP1PxTimesRealKindG(Cpr, Cpi, ai, Bpr, Bpi, n); /* C = (1 + ai*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqP1PxTimesRealKindT(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); /* C = (1 + ai*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqP1PxTimesRealKindC(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); /* C = (1 + ai*i) * (Bpr + Bpi * i)C */
            }
        }
    } else if( ar == -one ) {
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqM1P1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 + 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqM1P1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 + 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqM1P1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 + 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqM1P1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 + 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqM1M1TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 - 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqM1M1TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 - 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqM1M1TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 - 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqM1M1TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 - 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {
                RealKindEqM1P0TimesRealKindN(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 + 0*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqM1P0TimesRealKindG(Cpr, Cpi, Bpr, Bpi, n); /* C = (-1 + 0*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqM1P0TimesRealKindT(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 + 0*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqM1P0TimesRealKindC(Cpr, Cpi, Bpr, Bpi, m2, n2, p); /* C = (-1 + 0*i) * (Bpr + Bpi * i)C */
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqM1PxTimesRealKindN(Cpr, Cpi, ai, Bpr, Bpi, n); /* C = (-1 + ai*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqM1PxTimesRealKindG(Cpr, Cpi, ai, Bpr, Bpi, n); /* C = (-1 + ai*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqM1PxTimesRealKindT(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); /* C = (-1 + ai*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqM1PxTimesRealKindC(Cpr, Cpi, ai, Bpr, Bpi, m2, n2, p); /* C = (-1 + ai*i) * (Bpr + Bpi * i)C */
            }
        }
    } else {  /* ar != one && ar != -one */
        if( ai == one ) {
            if( transb == 'N' ) {
                RealKindEqPxP1TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar + 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqPxP1TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar + 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqPxP1TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar + 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqPxP1TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar + 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == -one ) {
            if( transb == 'N' ) {
                RealKindEqPxM1TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar - 1*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqPxM1TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar - 1*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqPxM1TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar - 1*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqPxM1TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar - 1*i) * (Bpr + Bpi * i)C */
            }
        } else if( ai == zero ) {
            if( transb == 'N' ) {
                RealKindEqPxP0TimesRealKindN(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar + 0*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqPxP0TimesRealKindG(Cpr, Cpi, ar, Bpr, Bpi, n); /* C = (ar + 0*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqPxP0TimesRealKindT(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar + 0*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqPxP0TimesRealKindC(Cpr, Cpi, ar, Bpr, Bpi, m2, n2, p); /* C = (ar + 0*i) * (Bpr + Bpi * i)C */
            }
        } else {
            if( transb == 'N' ) {
                RealKindEqPxPxTimesRealKindN(Cpr, Cpi, ar, ai, Bpr, Bpi, n); /* C = (ar + ai*i) * (Bpr + Bpi * i) */
            } else if( transb == 'G' ) {
                RealKindEqPxPxTimesRealKindG(Cpr, Cpi, ar, ai, Bpr, Bpi, n); /* C = (ar + ai*i) * (Bpr - Bpi * i) */
            } else if( transb == 'T' ) {
                RealKindEqPxPxTimesRealKindT(Cpr, Cpi, ar, ai, Bpr, Bpi, m2, n2, p); /* C = (ar + ai*i) * (Bpr + Bpi * i)T */
            } else { /* if( transb == 'C' ) { */
                RealKindEqPxPxTimesRealKindC(Cpr, Cpi, ar, ai, Bpr, Bpi, m2, n2, p); /* C = (ar + ai*i) * (Bpr + Bpi * i)C */
            }
        }
    }
}

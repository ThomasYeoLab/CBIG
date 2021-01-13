#include "fintrf.h"
C     glmnetMex.F
C     
C     Lasso and elastic-net regularized generalized linear models
      
C     [a0,ca,ia,nin,rsq,alm,nlp,jerr] = ...
C     glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,w,ka,cl,intr,maxit)
C     [a0,ca,ia,nin,dev,alm,nlp,jerr] = ...
C     glmnetMex(parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,nc,maxit,kopt)
C     
C     Extremely efficient procedures for fitting the entire lasso or
C     elastic-net regularization path for linear regression, logistic,  
C     multinomial, Poisson and Cox regression models. The algorithm uses cyclical
C     coordinate descent in a pathwise as described in the paper on the maintainer's
C     website.
C     
C     NOTES: This is a MEX-file wrapper of GLMnet.f for MATLAB. Should be called
C     only by glmnet.m. For details about input and output arguments, see
C     GLMnet.f.
C     
C     LICENSE: GPL-2
C     
C     DATE: 13 Jul 2009
C     
C     AUTHORS:
C     Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani
C     Fortran code was written by Jerome Friedman
C     R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
C     The original MATLAB wrapper was written by Hui Jiang (14 Jul 2009),
C     and was updated and is maintained by Junyang Qian (30 Aug 2013) junyangq@stanford.edu,
C     Department of Statistics, Stanford University, Stanford, California, USA.
C     
C     REFERENCES:
C     Jerome Friedman, Trevor Hastie and Rob Tibshirani. (2008). 
C     Regularization Paths for Generalized Linear Models via Coordinate Descent
C     Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010.

C     Noah Simon, Jerome Friedman, Trevor Hastie and Rob Tibshirani. (2011).
C     Regularization Paths for Cox’s Proportional Hazards Model via Coordinate Descent
C     Journal of Statistical Software, Vol. 39(5) 1-13.

C     Robert Tibshirani, Jacob Bien, Jerome Friedman, Trevor Hastie, Noah Simon, 
C     Jonathan Taylor, Ryan J. Tibshirani. (2010).
C     Strong Rules for Discarding Predictors in Lasso-type Problems
C     Journal of the Royal Statistical Society: Series B (Statistical Methodology), 74(2), 245-266.

C     Noah Simon, Jerome Friedman and Trevor Hastie (2013). 
C     A Blockwise Descent Algorithm for Group-penalized Multiresponse and Multinomial Regression 
C     (in arXiv, submitted)
C     
C     EXAMPLE:
C     task = 11;
C     parm = 1.0;
C     x = [1 1; 2 2; 3 3];
C     y = [1 3 2]';
C     jd = 0;
C     vp = [1 1]';
C     ne = 3;
C     nx = 2;
C     nlam = 100;
C     flmin = 0.0001;
C     ulam = 0;
C     thr = 1.0e-4;
C     isd = 1;
C     w = [1 1 1]';
C     ka = 2;
C     cl = [-Inf; Inf];
C     intr = 1;
C     maxit = 1E+5;
C     [lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr] = glmnetMex(task,parm,x,y,jd,vp,ne,nx,nlam,flmin,ulam,thr,
C     isd,w,ka,cl,intr,maxit)
C     
C     DEVELOPMENT:
C       14 Jul 2009: Original version of glmnet.m written.
C       30 Aug 2013: Updated glmnet.m with more options and more models
C                    (multi-response Gaussian, cox and Poisson models) supported.  

C     OLDER UPDATES:
C       26 Jan 2010: Fixed a bug in the description of y, pointed out by
C                    Peter Rijnbeek from Erasmus University.
C       09 Mar 2010: Fixed a bug of printing "ka = 2", pointed out by
C                    Ramon Casanova from Wake Forest University.
C       25 Mar 2010: Fixed a bug when p > n in multinomial fitting, pointed
C                    out by Gerald Quon from University of Toronto
C       25 Jul 2010: Check for input matrix format and size
C       27 Sep 2010: Fixed a bug of undefined "df" in multinomial fitting,
C                    pointed by Jeff Howbert from Insilicos.
C     
C-----------------------------------------------------------------------

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C-----------------------------------------------------------------------

      mwpointer plhs(*), prhs(*)
      mwpointer mxCreateDoubleMatrix, mxGetPr, mxCreateNumericArray
      integer nlhs, nrhs
      mwsize mxGetM, mxGetN, mxGetNzmax
      integer mxIsNumeric
      integer mxIsSparse
      
C-----------------------------------------------------------------------

C     Input
      real parm,flmin,thr, intr
      integer ka,no,ni,nr,nc,ne,nx,nlam,isd,maxit,kopt,isparse,nnz,jsd
      real, dimension (:), allocatable :: x,y,w,vp,ulam,cl,sr,xs,o,d,
     $     flog,a
      integer, dimension (:), allocatable :: ix,jx,jd,irs,jcs

      mwpointer pr

C     Output
      integer lmu,nlp,jerr
      real dev0
      real, dimension (:), allocatable :: a0,ca,alm,dev,rsq
      integer, dimension (:), allocatable :: ia,nin

C     Temporary      
      mwpointer temp_pr
      mwsize temp_m, temp_n, temp_nzmax, dims(3)
      integer task,i

C     For internal parameters
      real fdev, devmax, eps, big, pmin, prec, exmx
      integer mnlam, mxit
      
C     Check for proper number of arguments.
      if (nrhs .eq. 0) then
         task = -1;
      else
         temp_pr = mxGetPr(prhs(1))
         call getinteger(temp_pr,task,1)
      endif

C     Get input

      if (task .eq. -1) then
         call get_int_parms(fdev,eps,big,mnlam,devmax,pmin,exmx)
         call get_bnorm(prec,mxit)

         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putreal(fdev,temp_pr,1)

         plhs(2) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(devmax,temp_pr,1)

         plhs(3) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(3))
         call putreal(eps,temp_pr,1)

         plhs(4) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putreal(big,temp_pr,1)

         plhs(5) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putinteger(mnlam,temp_pr,1)

         plhs(6) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(pmin,temp_pr,1)

         plhs(7) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putreal(exmx,temp_pr,1)

         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putreal(prec,temp_pr,1)

         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putinteger(mxit,temp_pr,1)

         return   
      endif 
      
      if (task .eq. 0) then
         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,fdev,1)

         temp_pr = mxGetPr(prhs(3))
         call getreal(temp_pr,devmax,1)

         temp_pr = mxGetPr(prhs(4))
         call getreal(temp_pr,eps,1)

         temp_pr = mxGetPr(prhs(5))
         call getreal(temp_pr,big,1)

         temp_pr = mxGetPr(prhs(6))
         call getinteger(temp_pr,mnlam,1)

         temp_pr = mxGetPr(prhs(7))
         call getreal(temp_pr,pmin,1)

         temp_pr = mxGetPr(prhs(8))
         call getreal(temp_pr,exmx,1)

         temp_pr = mxGetPr(prhs(9))
         call getreal(temp_pr,prec,1)

         temp_pr = mxGetPr(prhs(10))
         call getinteger(temp_pr,mxit,1)

         call chg_fract_dev(fdev)
         call chg_dev_max(devmax)
         call chg_min_flmin(eps)
         call chg_big(big)
         call chg_min_lambdas(mnlam)
         call chg_min_null_prob(pmin)
         call chg_max_exp(exmx)
         call chg_bnorm(prec, mxit)

         return
      endif

c$$$  -----------------Gaussian--------------------  
c$$$  ---input---   
      
      if (task .eq. 10 .or. task .eq. 11) then
         if (task .eq. 11) then
            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            no = temp_m
            temp_n = mxGetN(prhs(3))
            ni = temp_n
            allocate(x(1:no*ni))
            call getreal(temp_pr,x,no*ni)
            
         else
            temp_m = mxGetM(prhs(4))
            no = temp_m

            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            nnz = temp_m
            allocate(xs(1:nnz))
            call getreal(temp_pr,xs,nnz)

            temp_pr = mxGetPr(prhs(19))
            allocate(irs(1:nnz))
            call getinteger(temp_pr,irs,nnz) 

            temp_pr = mxGetPr(prhs(20))
            temp_n = mxGetM(prhs(20))
            ni = temp_n - 1
            allocate(jcs(1:(ni+1)))
            call getinteger(temp_pr,jcs,(ni+1)) 
         endif

         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,parm,1)

         temp_pr = mxGetPr(prhs(4))
         allocate(y(1:no))
         call getreal(temp_pr,y,no)
         
         temp_pr = mxGetPr(prhs(5))
         temp_m = mxGetM(prhs(5))
         temp_n = mxGetN(prhs(5))
         allocate(jd(temp_m*temp_n))
         call getinteger(temp_pr,jd,temp_m*temp_n)     
         
         temp_pr = mxGetPr(prhs(6))
         allocate(vp(1:ni))
         call getreal(temp_pr,vp,ni)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,ne,1)

         temp_pr = mxGetPr(prhs(8))
         call getinteger(temp_pr,nx,1)
         
         temp_pr = mxGetPr(prhs(9))
         call getinteger(temp_pr,nlam,1)

         temp_pr = mxGetPr(prhs(10))
         call getreal(temp_pr,flmin,1)     
         
         temp_pr = mxGetPr(prhs(11))
         temp_m = mxGetM(prhs(11))
         temp_n = mxGetN(prhs(11))
         allocate(ulam(1:temp_m * temp_n))
         call getreal(temp_pr,ulam,temp_m * temp_n)
         
         temp_pr = mxGetPr(prhs(12))
         call getreal(temp_pr,thr,1)
         
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,isd,1)

         temp_pr = mxGetPr(prhs(14))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)

         temp_pr = mxGetPr(prhs(15))
         call getinteger(temp_pr,ka,1)

         temp_pr = mxGetPr(prhs(16))
         allocate(cl(1:2*ni))
         call getreal(temp_pr,cl,2*ni)

         temp_pr = mxGetPr(prhs(17))
         call getinteger(temp_pr,intr,1)

         temp_pr = mxGetPr(prhs(18))
         call getinteger(temp_pr,maxit,1)       

c$$$  ---prepare output---

         allocate(ia(1:nx))
         call zerointeger(ia,nx)
         allocate(nin(1:nlam))
         call zerointeger(nin,nlam)
         allocate(alm(1:nlam))
         call zeroreal(alm,nlam)
         allocate(a0(1:nlam))
         call zeroreal(a0,nlam)         
         allocate(ca(1:nx*nlam))
         call zeroreal(ca,nx*nlam)
         allocate(rsq(1:nlam))
         call zeroreal(rsq,nlam)


c$$$  ---computation----

         if (task .eq. 11) then    
            call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,
     $           ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     $           nlp,jerr)
         else
            call spelnet(ka,parm,no,ni,xs,jcs,irs,y,w,jd,vp,cl,ne,nx,
     $           nlam,flmin,ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,
     $           rsq,alm,nlp,jerr)
         endif

c$$$  ----output-----

         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putinteger(lmu,temp_pr,1)

         plhs(4) = mxCreateDoubleMatrix(nx,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putinteger(ia,temp_pr,nx)
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putinteger(nin,temp_pr,lmu)
         
         plhs(7) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putreal(alm,temp_pr,lmu)
         
         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putinteger(nlp,temp_pr,1)
         
         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putinteger(jerr,temp_pr,1)

         plhs(2) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(a0,temp_pr,lmu)   
         
         plhs(3) = mxCreateDoubleMatrix(nx,lmu,0)
         temp_pr = mxGetPr(plhs(3))
         call putreal(ca,temp_pr,nx*lmu)  
         
         plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(rsq,temp_pr,lmu)  

         deallocate(y)
         deallocate(jd)
         deallocate(vp)
         deallocate(ulam)
         deallocate(a0)
         deallocate(ca)
         deallocate(ia)
         deallocate(nin)
         deallocate(alm)
         deallocate(w)
         deallocate(rsq)
         deallocate(cl)
         
         if (task .eq. 11) then
            deallocate(x)
         else
            deallocate(xs)
            deallocate(irs)
            deallocate(jcs)  
         endif
         return
      endif      
c$$$  --------------end of Gaussian---------------------------

c$$$  -----------------mGaussian--------------------  
c$$$  ---input---   
      
      if (task .eq. 30 .or. task .eq. 31) then
         if (task .eq. 31) then
            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            no = temp_m
            temp_n = mxGetN(prhs(3))
            ni = temp_n
            allocate(x(1:no*ni))
            call getreal(temp_pr,x,no*ni)

            temp_pr = mxGetPr(prhs(18))
            call getinteger(temp_pr,jsd,1)
            
         else
            temp_m = mxGetM(prhs(4))
            no = temp_m

            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            nnz = temp_m
            allocate(xs(1:nnz))
            call getreal(temp_pr,xs,nnz)

            temp_pr = mxGetPr(prhs(18))
            allocate(irs(1:nnz))
            call getinteger(temp_pr,irs,nnz) 

            temp_pr = mxGetPr(prhs(19))
            temp_n = mxGetM(prhs(19))
            ni = temp_n - 1
            allocate(jcs(1:(ni+1)))
            call getinteger(temp_pr,jcs,(ni+1)) 

            temp_pr = mxGetPr(prhs(20))
            call getinteger(temp_pr,jsd,1)
         endif

         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,parm,1)

         temp_n = mxGetN(prhs(4))
         nr = temp_n
         temp_pr = mxGetPr(prhs(4))
         allocate(y(1:no*nr))
         call getreal(temp_pr,y,no*nr)
         
         temp_pr = mxGetPr(prhs(5))
         temp_m = mxGetM(prhs(5))
         temp_n = mxGetN(prhs(5))
         allocate(jd(temp_m*temp_n))
         call getinteger(temp_pr,jd,temp_m*temp_n)     
         
         temp_pr = mxGetPr(prhs(6))
         allocate(vp(1:ni))
         call getreal(temp_pr,vp,ni)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,ne,1)

         temp_pr = mxGetPr(prhs(8))
         call getinteger(temp_pr,nx,1)
         
         temp_pr = mxGetPr(prhs(9))
         call getinteger(temp_pr,nlam,1)

         temp_pr = mxGetPr(prhs(10))
         call getreal(temp_pr,flmin,1)     
         
         temp_pr = mxGetPr(prhs(11))
         temp_m = mxGetM(prhs(11))
         temp_n = mxGetN(prhs(11))
         allocate(ulam(1:temp_m * temp_n))
         call getreal(temp_pr,ulam,temp_m * temp_n)
         
         temp_pr = mxGetPr(prhs(12))
         call getreal(temp_pr,thr,1)
         
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,isd,1)

         temp_pr = mxGetPr(prhs(14))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)

         temp_pr = mxGetPr(prhs(15))
         allocate(cl(1:2*ni))
         call getreal(temp_pr,cl,2*ni)

         temp_pr = mxGetPr(prhs(16))
         call getinteger(temp_pr,intr,1)

         temp_pr = mxGetPr(prhs(17))
         call getinteger(temp_pr,maxit,1) 
      

c$$$  ---prepare output---

         allocate(ia(1:nx))
         call zerointeger(ia,nx)
         allocate(nin(1:nlam))
         call zerointeger(nin,nlam)
         allocate(alm(1:nlam))
         call zeroreal(alm,nlam)
         allocate(a0(1:nr*nlam))
         call zeroreal(a0,nr*nlam)        
         allocate(ca(1:nx*nr*nlam))
         call zeroreal(ca,nx*nr*nlam)        
         allocate(rsq(1:nlam))
         call zeroreal(rsq,nlam)

c$$$  ---computation----

         if (task .eq. 31) then    
            call multelnet(parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,
     $           flmin,ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,
     $           rsq,alm,nlp,jerr)
         else
            call multspelnet(parm,no,ni,nr,xs,jcs,irs,y,w,jd,vp,cl,ne,
     $           nx,nlam,flmin,ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,
     $           ia,nin,rsq,alm,nlp,jerr)
         endif

c$$$  ----output-----

         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putinteger(lmu,temp_pr,1)

         plhs(4) = mxCreateDoubleMatrix(nx,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putinteger(ia,temp_pr,nx)
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putinteger(nin,temp_pr,lmu)
         
         plhs(7) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putreal(alm,temp_pr,lmu)
         
         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putinteger(nlp,temp_pr,1)
         
         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putinteger(jerr,temp_pr,1)

         plhs(2) = mxCreateDoubleMatrix(nr,lmu,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(a0,temp_pr,nr*lmu)   
         
         plhs(3) = mxCreateDoubleMatrix(nx*nr,lmu,0)
         temp_pr = mxGetPr(plhs(3))
         call putreal(ca,temp_pr,nx*nr*lmu)  
         
         plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(rsq,temp_pr,lmu)  

         deallocate(y)
         deallocate(jd)
         deallocate(vp)
         deallocate(ulam)
         deallocate(a0)
         deallocate(ca)
         deallocate(ia)
         deallocate(nin)
         deallocate(alm)
         deallocate(w)
         deallocate(rsq)
         deallocate(cl)
         
         if (task .eq. 31) then
            deallocate(x)
         else
            deallocate(xs)
            deallocate(irs)
            deallocate(jcs)  
         endif
         return
      endif      
c$$$  --------------end of mGaussian---------------------------

c$$$  ---------------Logistic--------------------------
c$$$  ---input---   
      
      if (task .eq. 20 .or. task .eq. 21) then
         if (task .eq. 21) then
            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            no = temp_m
            temp_n = mxGetN(prhs(3))
            ni = temp_n
            allocate(x(1:no*ni))
            call getreal(temp_pr,x,no*ni)
            
         else
            temp_m = mxGetM(prhs(4))
            no = temp_m

            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            nnz = temp_m
            allocate(xs(1:nnz))
            call getreal(temp_pr,xs,nnz)

            temp_pr = mxGetPr(prhs(20))
            allocate(irs(1:nnz))
            call getinteger(temp_pr,irs,nnz) 

            temp_pr = mxGetPr(prhs(21))
            temp_n = mxGetM(prhs(21))
            ni = temp_n - 1
            allocate(jcs(1:(ni+1)))
            call getinteger(temp_pr,jcs,(ni+1)) 
         endif

         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,parm,1)

         temp_pr = mxGetPr(prhs(17))
         call getinteger(temp_pr,nc,1)

         temp_pr = mxGetPr(prhs(4))
         allocate(y(1:no*(max(2,nc))))
         call getreal(temp_pr,y,no*(max(2,nc)))
         
         temp_pr = mxGetPr(prhs(5))
         temp_m = mxGetM(prhs(5))
         temp_n = mxGetN(prhs(5))
         allocate(jd(temp_m*temp_n))
         call getinteger(temp_pr,jd,temp_m*temp_n)     
         
         temp_pr = mxGetPr(prhs(6))
         allocate(vp(1:ni))
         call getreal(temp_pr,vp,ni)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,ne,1)

         temp_pr = mxGetPr(prhs(8))
         call getinteger(temp_pr,nx,1)
         
         temp_pr = mxGetPr(prhs(9))
         call getinteger(temp_pr,nlam,1)

         temp_pr = mxGetPr(prhs(10))
         call getreal(temp_pr,flmin,1)     
         
         temp_pr = mxGetPr(prhs(11))
         temp_m = mxGetM(prhs(11))
         temp_n = mxGetN(prhs(11))
         allocate(ulam(1:temp_m * temp_n))
         call getreal(temp_pr,ulam,temp_m * temp_n)
         
         temp_pr = mxGetPr(prhs(12))
         call getreal(temp_pr,thr,1)
         
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,isd,1)

         temp_pr = mxGetPr(prhs(14))
         allocate(cl(1:2*ni))
         call getreal(temp_pr,cl,2*ni)

         temp_pr = mxGetPr(prhs(15))
         call getinteger(temp_pr,intr,1)

         temp_pr = mxGetPr(prhs(16))
         call getinteger(temp_pr,maxit,1)

         temp_pr = mxGetPr(prhs(18))
         call getinteger(temp_pr,kopt,1)

         temp_pr = mxGetPr(prhs(19))
         allocate(o(1:no*nc))
         call getreal(temp_pr,o,no*nc)

c$$$  ---prepare output---

         allocate(ia(1:nx))
         call zerointeger(ia,nx)
         allocate(nin(1:nlam))
         call zerointeger(nin,nlam)
         allocate(alm(1:nlam))
         call zeroreal(alm,nlam)
         allocate(a0(1:(nc*nlam)))
         call zeroreal(a0,nc*nlam)         
         allocate(ca(1:(nx*nc*nlam)))
         call zeroreal(ca,nx*nc*nlam)
         allocate(dev(1:nlam))
         call zeroreal(dev,nlam)

c$$$  ---computation---- 

         if (task .eq. 21) then    
            call lognet(parm,no,ni,nc,x,y,o,jd,vp,cl,ne,nx,nlam,flmin,
     $           ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,
     $           alm,nlp,jerr)
         else
            call splognet(parm,no,ni,nc,xs,jcs,irs,y,o,jd,vp,cl,ne,nx,
     $           nlam,flmin,ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,
     $           nin,dev0,dev,alm,nlp,jerr)
         endif

c$$$  ----output-----

         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putinteger(lmu,temp_pr,1)

         plhs(4) = mxCreateDoubleMatrix(nx,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putinteger(ia,temp_pr,nx)
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putinteger(nin,temp_pr,lmu)
         
         plhs(7) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putreal(alm,temp_pr,lmu)
         
         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putinteger(nlp,temp_pr,1)
         
         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putinteger(jerr,temp_pr,1)

         plhs(2) = mxCreateDoubleMatrix(nc,lmu,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(a0,temp_pr,nc*lmu)

         plhs(3) = mxCreateDoubleMatrix(nx*nc,lmu,0)
         temp_pr = mxGetPr(plhs(3))
         call putreal(ca,temp_pr,nx*nc*lmu)
         
         plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(dev,temp_pr,lmu)

         plhs(10) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(10))
         call putreal(dev0,temp_pr,1)

         plhs(11) = mxCreateDoubleMatrix(no,nc,0)
         temp_pr = mxGetPr(plhs(11))
         call putreal(o,temp_pr,no*nc)

         deallocate(y)
         deallocate(jd)
         deallocate(vp)
         deallocate(ulam)
         deallocate(a0)
         deallocate(ca)
         deallocate(ia)
         deallocate(nin)
         deallocate(alm)
         deallocate(cl)
         deallocate(o)
         deallocate(dev)
         
         if (task .eq. 21) then
            deallocate(x)
         else
            deallocate(xs)
            deallocate(irs)
            deallocate(jcs)  
         endif
         return
      endif

c$$$  --------------------end of Logistic------------------

c$$$  --------------------start of cox--------------------

      if (task .eq. 41) then
         temp_pr = mxGetPr(prhs(3))
         temp_m = mxGetM(prhs(3))
         no = temp_m
         temp_n = mxGetN(prhs(3))
         ni = temp_n
         allocate(x(1:no*ni))
         call getreal(temp_pr,x,no*ni)
         

         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,parm,1)

         temp_pr = mxGetPr(prhs(4))
         allocate(y(1:no))
         call getreal(temp_pr,y,no)
         
         temp_pr = mxGetPr(prhs(5))
         temp_m = mxGetM(prhs(5))
         temp_n = mxGetN(prhs(5))
         allocate(jd(temp_m*temp_n))
         call getinteger(temp_pr,jd,temp_m*temp_n)     
         
         temp_pr = mxGetPr(prhs(6))
         allocate(vp(1:ni))
         call getreal(temp_pr,vp,ni)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,ne,1)

         temp_pr = mxGetPr(prhs(8))
         call getinteger(temp_pr,nx,1)
         
         temp_pr = mxGetPr(prhs(9))
         call getinteger(temp_pr,nlam,1)

         temp_pr = mxGetPr(prhs(10))
         call getreal(temp_pr,flmin,1)     
         
         temp_pr = mxGetPr(prhs(11))
         temp_m = mxGetM(prhs(11))
         temp_n = mxGetN(prhs(11))
         allocate(ulam(1:temp_m * temp_n))
         call getreal(temp_pr,ulam,temp_m * temp_n)
         
         temp_pr = mxGetPr(prhs(12))
         call getreal(temp_pr,thr,1)
         
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,isd,1)

         temp_pr = mxGetPr(prhs(14))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)

         temp_pr = mxGetPr(prhs(15))
         allocate(d(1:no))
         call getreal(temp_pr,d,no)

         temp_pr = mxGetPr(prhs(16))
         allocate(cl(1:2*ni))
         call getreal(temp_pr,cl,2*ni)

         temp_pr = mxGetPr(prhs(17))
         call getinteger(temp_pr,maxit,1)

         temp_pr = mxGetPr(prhs(18))
         allocate(o(1:no))
         call getreal(temp_pr,o,no)       

c$$$  ---prepare output---

         allocate(ia(1:nx))
         call zerointeger(ia,nx)
         allocate(nin(1:nlam))
         call zerointeger(nin,nlam)
         allocate(alm(1:nlam))
         call zeroreal(alm,nlam)         
         allocate(ca(1:nx*nlam))
         call zeroreal(ca,nx*nlam)       
         allocate(dev(1:nlam))
         call zeroreal(dev,nlam)

c$$$  ---computation----

         call coxnet(parm,no,ni,x,y,d,o,w,jd,vp,cl,ne,nx,nlam,flmin,
     $        ulam,thr,maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)

c$$$  ----output-----
c$$$
         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putinteger(lmu,temp_pr,1)

         plhs(3) = mxCreateDoubleMatrix(nx,1,0)
         temp_pr = mxGetPr(plhs(3))
         call putinteger(ia,temp_pr,nx)
         
         plhs(4) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putinteger(nin,temp_pr,lmu)
         
         plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(alm,temp_pr,lmu)
         
         plhs(7) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putinteger(nlp,temp_pr,1)
         
         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putinteger(jerr,temp_pr,1) 
         
         plhs(2) = mxCreateDoubleMatrix(nx,lmu,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(ca,temp_pr,nx*lmu)  
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putreal(dev,temp_pr,lmu)  

         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putreal(dev0,temp_pr,1)

         plhs(10) = mxCreateDoubleMatrix(no,1,0)
         temp_pr = mxGetPr(plhs(10))
         call putreal(o,temp_pr,no)

         deallocate(y)
         deallocate(jd)
         deallocate(vp)
         deallocate(ulam)
         deallocate(ca)
         deallocate(ia)
         deallocate(nin)
         deallocate(alm)
         deallocate(w)
         deallocate(dev)
         deallocate(cl)
         deallocate(d)
         deallocate(o)
         deallocate(x)

         return
      endif      

c$$$  -------------------end of cox -----------------------

c$$$  --------------start of log-likelihood----------------
      if (task .eq. 42) then
c$$$  prepare input 
         temp_pr = mxGetPr(prhs(2))
         temp_m = mxGetM(prhs(2))
         no = temp_m
         temp_n = mxGetN(prhs(2))
         ni = temp_n
         allocate(x(1:no*ni))
         call getreal(temp_pr,x,no*ni)

         temp_pr = mxGetPr(prhs(3))
         allocate(y(1:no))
         call getreal(temp_pr,y,no)

         temp_pr = mxGetPr(prhs(4))
         allocate(d(1:no))
         call getreal(temp_pr,d,no)

         temp_pr = mxGetPr(prhs(5))
         allocate(o(1:no))
         call getreal(temp_pr,o,no)

         temp_pr = mxGetPr(prhs(6))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,nvec,1)

         temp_pr = mxGetPr(prhs(8))
         allocate(a(1:ni*nvec))
         call getreal(temp_pr,a,ni*nvec)

c$$$  prepare output
         allocate(flog(1:nvec))
         call zeroreal(flog,nvec)

c$$$  compututation

         call loglike(no,ni,x,y,d,o,w,nvec,a,flog,jerr)

c$$$  output
         plhs(1) = mxCreateDoubleMatrix(nvec,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putreal(flog,temp_pr,nvec)

         plhs(2) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(2))
         call putinteger(jerr,temp_pr,1) 
      endif


c$$$  ---------------Poisson--------------------------
c$$$  ---input---   
      
      if (task .eq. 50 .or. task .eq. 51) then
         if (task .eq. 51) then
            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            no = temp_m
            temp_n = mxGetN(prhs(3))
            ni = temp_n
            allocate(x(1:no*ni))
            call getreal(temp_pr,x,no*ni)
            
         else
            temp_m = mxGetM(prhs(4))
            no = temp_m

            temp_pr = mxGetPr(prhs(3))
            temp_m = mxGetM(prhs(3))
            nnz = temp_m
            allocate(xs(1:nnz))
            call getreal(temp_pr,xs,nnz)

            temp_pr = mxGetPr(prhs(19))
            allocate(irs(1:nnz))
            call getinteger(temp_pr,irs,nnz) 

            temp_pr = mxGetPr(prhs(20))
            temp_n = mxGetM(prhs(20))
            ni = temp_n - 1
            allocate(jcs(1:(ni+1)))
            call getinteger(temp_pr,jcs,(ni+1)) 
         endif

         temp_pr = mxGetPr(prhs(2))
         call getreal(temp_pr,parm,1)

         temp_pr = mxGetPr(prhs(4))
         allocate(y(1:no))
         call getreal(temp_pr,y,no)
         
         temp_pr = mxGetPr(prhs(5))
         temp_m = mxGetM(prhs(5))
         temp_n = mxGetN(prhs(5))
         allocate(jd(temp_m*temp_n))
         call getinteger(temp_pr,jd,temp_m*temp_n)     
         
         temp_pr = mxGetPr(prhs(6))
         allocate(vp(1:ni))
         call getreal(temp_pr,vp,ni)

         temp_pr = mxGetPr(prhs(7))
         call getinteger(temp_pr,ne,1)

         temp_pr = mxGetPr(prhs(8))
         call getinteger(temp_pr,nx,1)
         
         temp_pr = mxGetPr(prhs(9))
         call getinteger(temp_pr,nlam,1)

         temp_pr = mxGetPr(prhs(10))
         call getreal(temp_pr,flmin,1)     
         
         temp_pr = mxGetPr(prhs(11))
         temp_m = mxGetM(prhs(11))
         temp_n = mxGetN(prhs(11))
         allocate(ulam(1:temp_m * temp_n))
         call getreal(temp_pr,ulam,temp_m * temp_n)
         
         temp_pr = mxGetPr(prhs(12))
         call getreal(temp_pr,thr,1)
         
         temp_pr = mxGetPr(prhs(13))
         call getinteger(temp_pr,isd,1)

         temp_pr = mxGetPr(prhs(14))
         allocate(w(1:no))
         call getreal(temp_pr,w,no)

         temp_pr = mxGetPr(prhs(15))
         allocate(cl(1:2*ni))
         call getreal(temp_pr,cl,2*ni)

         temp_pr = mxGetPr(prhs(16))
         call getinteger(temp_pr,intr,1)

         temp_pr = mxGetPr(prhs(17))
         call getinteger(temp_pr,maxit,1)

         temp_pr = mxGetPr(prhs(18))
         allocate(o(1:no))
         call getreal(temp_pr,o,no)

c$$$  ---prepare output---

         allocate(ia(1:nx))
         call zerointeger(ia,nx)
         allocate(nin(1:nlam))
         call zerointeger(nin,nlam)
         allocate(alm(1:nlam))
         call zeroreal(alm,nlam)
         allocate(a0(1:nlam))
         call zeroreal(a0,nlam)
         allocate(ca(1:nx*nlam))
         call zeroreal(ca,nx*nlam)
         allocate(dev(1:nlam))
         call zeroreal(dev,nlam)

c$$$  ---computation----

         if (task .eq. 51) then    
            call fishnet(parm,no,ni,x,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,
     $           ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,
     $           nlp,jerr)
         else
            call spfishnet(parm,no,ni,xs,jcs,irs,y,o,w,jd,vp,cl,ne,nx,
     $           nlam,flmin,ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,
     $           nin,dev0,dev,alm,nlp,jerr)
         endif

c$$$  ----output-----

         plhs(1) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(1))
         call putinteger(lmu,temp_pr,1)

         plhs(4) = mxCreateDoubleMatrix(nx,1,0)
         temp_pr = mxGetPr(plhs(4))
         call putinteger(ia,temp_pr,nx)
         
         plhs(5) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(5))
         call putinteger(nin,temp_pr,lmu)
         
         plhs(7) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(7))
         call putreal(alm,temp_pr,lmu)
         
         plhs(8) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(8))
         call putinteger(nlp,temp_pr,1)
         
         plhs(9) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(9))
         call putinteger(jerr,temp_pr,1)

         plhs(2) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(2))
         call putreal(a0,temp_pr,lmu)

         plhs(3) = mxCreateDoubleMatrix(nx,lmu,0)
         temp_pr = mxGetPr(plhs(3))
         call putreal(ca,temp_pr,nx*lmu)
         
         plhs(6) = mxCreateDoubleMatrix(lmu,1,0)
         temp_pr = mxGetPr(plhs(6))
         call putreal(dev,temp_pr,lmu)

         plhs(10) = mxCreateDoubleMatrix(1,1,0)
         temp_pr = mxGetPr(plhs(10))
         call putreal(dev0,temp_pr,1)

         plhs(11) = mxCreateDoubleMatrix(no,1,0)
         temp_pr = mxGetPr(plhs(11))
         call putreal(o,temp_pr,no)

         deallocate(y)
         deallocate(jd)
         deallocate(vp)
         deallocate(ulam)
         deallocate(a0)
         deallocate(ca)
         deallocate(ia)
         deallocate(nin)
         deallocate(alm)
         deallocate(cl)
         deallocate(o)
         deallocate(dev)
         
         if (task .eq. 51) then
            deallocate(x)
         else
            deallocate(xs)
            deallocate(irs)
            deallocate(jcs)  
         endif
         return
      endif

c$$$  --------------------end of Poisson------------------

      return
      end

C     End of subroutine mexFunction
      
      subroutine real8toreal(x, y, size)
      integer size
      real*8 x(size)
      real y(size)
      do 10 i=1,size
         y(i)= x(i)
 10   continue
      return
      end

      subroutine realtoreal8(x, y, size)
      integer size
      real x(size)
      real*8 y(size)
      do 20 i=1,size
         y(i)= x(i)
 20   continue
      return
      end
      
      subroutine real8tointeger(x, y, size)
      integer size
      real*8 x(size)
      integer y(size)
      do 30 i=1,size
         y(i)= x(i)
 30   continue
      return
      end
      
      subroutine integertoreal8(x, y, size)
      integer size
      integer x(size)
      real*8 y(size)
      do 40 i=1,size
         y(i)= x(i)
 40   continue
      return
      end
      
      subroutine getreal(pr,x,size)
      mwpointer pr
      integer size
      real x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call mxCopyPtrToReal8(pr,temp,size)
      call real8toreal(temp,x,size)
      deallocate(temp)      
      return
      end
      
      subroutine getinteger(pr,x,size)
      mwpointer pr
      integer size
      integer x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call mxCopyPtrToReal8(pr,temp,size)
      call real8tointeger(temp,x,size)
      deallocate(temp)      
      return
      end      
      
      subroutine putreal(x,pr,size)
      mwpointer pr
      integer size
      real x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call realtoreal8(x,temp,size)
      call mxCopyReal8ToPtr(temp,pr,size)
      deallocate(temp)      
      return
      end
      
      subroutine putinteger(x,pr,size)
      mwpointer pr
      integer size
      integer x(size)
      real*8, dimension (:), allocatable :: temp
      allocate(temp(1:size))
      call integertoreal8(x,temp,size)
      call mxCopyReal8ToPtr(temp,pr,size)
      deallocate(temp)      
      return
      end            
      
      subroutine zeroreal(x,size)
      integer size
      real x(size)
      do 90 i=1,size
         x(i) = 0
 90   continue     
      return 
      end
      
      subroutine zerointeger(x,size)
      integer size
      integer x(size)
      do 100 i=1,size
         x(i) = 0
 100  continue
      return
      end

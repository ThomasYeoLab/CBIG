c
c                          newGLMnet (2/15/13)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c            intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = penalty member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   y(no) = response vector (overwritten)
c   w(no)= observation weights (overwritten)
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   cl(2,ni) = interval constraints on coefficient values (overwritten)
c      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
c      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum reduction in the criterion value
c      as a result of each parameter update over a single pass
c      is less than thr times the null criterion value.
c      (suggested value, thr=1.0e-5)
c   isd = predictor variable standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   intr = intercept flag
c      intr = 0/1 => don't/do include intercept in model
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested values, maxit = 100000)
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = actual number of passes over the data for all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c
c
c
c least-squares utility routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx = input to elnet
c    lmu,ca,ia,nin = output from elnet
c
c output:
c
c    b(ni,lmu) = all elnet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c                           Multiple response
c                  elastic net with squared-error loss
c
c dense predictor matrix:
c
c call multelnet(parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c                jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call multspelnet(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   nr = number of response variables
c   y(no,nr) = response data matrix (overwritten)
c   jsd = response variable standardization flag
c      jsd = 0 => regression using original response variables
c      jsd = 1 => regression using standardized response variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   all other inputs same as elnet/spelnet above
c
c output:
c
c   a0(nr,lmu) = intercept values for each solution
c   ca(nx,nr,lmu) = compressed coefficient values for each solution
c   all other outputs same as elnet/spelnet above
c   (jerr = 90000 => bounds adjustment non convergence)
c
c
c
c multiple response least-squares utility routines:
c
c
c uncompress coefficient matrix for all solutions:
c
c call multsolns(ni,nx,nr,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,nr = input to multelnet
c    lmu,ca,ia,nin = output from multelnet
c
c output:
c
c    b(ni,nr,lmu) = all multelnet returned solutions in uncompressed format
c
c
c uncompress coefficient matrix for particular solution:
c
c call multuncomp(ni,nr,nx,ca,ia,nin,a)
c
c input:
c
c    ni,nr,nx = input to multelnet
c    ca(nx,nr) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni,nr) =  uncompressed coefficient matrix
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call multmodval(nx,nr,a0,ca,ia,nin,n,x,f);
c
c input:
c
c    nx,nr = input to multelnet
c    a0(nr) = intercepts
c    ca(nx,nr) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(nr,n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call multcmodval(nx,nr,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nx,nr = input to multelnet
c    a0(nr) = intercepts
c    ca(nx,nr) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nr,n) = model predictions
c
c
c
c
c          Symmetric binomial/multinomial logistic elastic net
c
c
c dense predictor matrix:
c
c call lognet (parm,no,ni,nc,x,y,o,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,o,jd,vp,cl,ne,nx,nlam,flmin,
c      ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   parm,no,ni,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit
c    = same as elnet above.
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
c   o(no,nc) = observation off-sets for each class
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson (recommended)
c      kpot = 1 => modified Newton-Raphson (sometimes faster)
c      kpot = 2 => nonzero coefficients same for each class (nc > 1)
c
c
c output:
c
c   lmu,ia,nin,alm,nlp = same as elent above
c
c   a0(nc,lmu) = intercept values for each class at each solution
c   ca(nx,nc,lmu) = compressed coefficient values for each class at
c                each solution
c   dev0 = null deviance (intercept only model)
c   fdev(lmu) = fraction of devience explained by each solution
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8000 + k => null probability < 1.0e-5 for class k
c         jerr = 9000 + k => null probability for class k
c                            > 1.0 - 1.0e-5
c         jerr = 10000 => maxval(vp) <= 0.0
c         jerr = 90000 => bounds adjustment non convergence
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c         jerr = -20000-k => max(p*(1-p)) < 1.0e-6 at kth lamda value.
c    o(no,nc) = training data values for last (lmu_th) solution linear
c               combination.
c
c
c
c logistic/multinomial utilitity routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call lsolns(ni,nx,nc,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,nc = input to lognet
c    lmu,ca,ia,nin = output from lognet
c
c output:
c
c    b(ni,nc,lmu) = all lognet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call luncomp(ni,nx,nc,ca,ia,nin,a)
c
c input:
c
c    ni, nx, nc = same as above
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c    a(ni,nc) =  uncompressed coefficient vectors
c                 referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor vectors:
c
c call lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
c
c input:
c
c    nt = number of observations
c    x(nt,ni) = full (uncompressed) predictor vectors
c    nc, nx = same as above
c    a0(nc) = intercepts
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c ans(nc,nt) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nc, nx = same as above
c    a0(nc) = intercept
c    ca(nx,nc) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nc,n) = model predictions
c
c
c
c
c                        Poisson elastic net
c
c
c dense predictor matrix:
c
c call fishnet (parm,no,ni,x,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c               isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c sparse predictor matrix:
c
c call spfishnet (parm,no,ni,x,ix,jx,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c               isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c    x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit
c    = same as elnet above
c
c output:
c
c   lmu,a0,ca,ia,nin,alm = same as elnet above
c   dev0,fdev = same as lognet above
c   nlp = total number of passes over predictor variables
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => negative response count y values
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c
c Poisson utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c    call modval(a0,ca,ia,nin,n,x,f);
c    call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c compute deviance for given uncompressed data and set of uncompressed
c solutions
c
c call deviance(no,ni,x,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output:
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and set of uncompressed solutions
c
c call spdeviance(no,ni,x,ix,jx,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and compressed solutions
c
c call cspdeviance(no,x,ix,jx,y,o,w,nx,lmu,a0,ca,ia,nin,flog,jerr)
c
c input:
c
c   no = number of observations
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nx = input to spfishnet
c   lmu,a0(lmu),ca(nx,lmu),ia(nx),nin(lmu) = output from spfishnet
c
c output
c
c   flog(lmu) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c
c          Elastic net with Cox proportional hazards model
c
c
c dense predictor matrix:
c
c call coxnet (parm,no,ni,x,y,d,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c              maxit,isd,lmu,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c input:
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,maxit
c                = same as fishnet above
c
c output:
c
c   lmu,ca,ia,nin,dev0,fdev,alm,nlp = same as fishnet above
c   jerr = error flag
c      jerr = 0  => no error - output returned
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => all observations censored (d(i)=0.0)
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
c         jerr = 20000, 30000 => initialization numerical error
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c         jerr = -30000-k => numerical error at kth lambda value
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c
c
c coxnet utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call cxmodval(ca,ia,nin,n,x,f);
c
c input:
c
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c compute log-likelihood for given data set and vectors of coefficients
c
c call loglike(no,ni,x,y,d,o,w,nvec,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nvec = number of coefficient vectors
c   a(ni,nvec) = coefficient vectors (uncompressed)
c
c output
c
c   flog(nvec) = respective log-likelihood values
c   jerr = error flag - see coxnet above
c
c
c
c
c                Changing internal parameter values
c
c
c call chg_fract_dev(fdev)
c   fdev = minimum fractional change in deviance for stopping path
c      default = 1.0e-5
c
c call chg_dev_max(devmax)
c   devmax = maximum fraction of explained deviance for stopping path
c      default = 0.999
c
c call chg_min_flmin(eps)
c   eps = minimum value of flmin (see above). default= 1.0e-6
c
c call chg_big(big)
c   big = large floating point number. default = 9.9e35
c
c call chg_min_lambdas(mnlam)
c   mnlam = minimum number of path points (lambda values) allowed
c      default = 5
c
c call chg_min_null_prob(pmin)
c   pmin = minimum null probability for any class. default = 1.0e-5
c
c call chg _max_exp(exmx)
c   exmx = maximum allowed exponent. default = 250.0
c
c call chg_bnorm(prec,mxit)
c   prec = convergence threshold for multi response bounds adjustment
c          solution. default = 1.0e-10.
c   mxit = maximum iterations for multiresponse bounds adjustment solution
c          default = 100.
c
c
c             Obtain current internal parameter values
c
c call get_int_parms(fdev,eps,big,mnlam,devmax,pmin,exmx)
c call get_bnorm(prec,mxit);
c
c
c             
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          771
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0e-5,1.0e-6,9.9    773 
     *e35,5,0.999,1.0e-5,250.0/
      sml=sml0                                                              773
      eps=eps0                                                              773
      big=big0                                                              773
      mnlam=mnlam0                                                          773
      rsqmax=rsqmax0                                                        774
      pmin=pmin0                                                            774
      exmx=exmx0                                                            775
      return                                                                776
      entry chg_fract_dev(arg)                                              776
      sml0=arg                                                              776
      return                                                                777
      entry chg_dev_max(arg)                                                777
      rsqmax0=arg                                                           777
      return                                                                778
      entry chg_min_flmin(arg)                                              778
      eps0=arg                                                              778
      return                                                                779
      entry chg_big(arg)                                                    779
      big0=arg                                                              779
      return                                                                780
      entry chg_min_lambdas(irg)                                            780
      mnlam0=irg                                                            780
      return                                                                781
      entry chg_min_null_prob(arg)                                          781
      pmin0=arg                                                             781
      return                                                                782
      entry chg_max_exp(arg)                                                782
      exmx0=arg                                                             782
      return                                                                783
      end                                                                   784
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,u    787 
     *lam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)                 788
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          789
      integer jd(*),ia(nx),nin(nlam)                                        790
      real, dimension (:), allocatable :: vq; 
      if(maxval(vp) .gt. 0.0)goto 10021                                     793
      jerr=10000                                                            793
      return                                                                793
10021 continue                                                              794
      allocate(vq(1:ni),stat=jerr)                                          794
      if(jerr.ne.0) return                                                  795
      vq=max(0.0,vp)                                                        795
      vq=vq*ni/sum(vq)                                                      796
      if(ka .ne. 1)goto 10041                                               797
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,    800 
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            801
10041 continue                                                              802
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i    805 
     *sd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              806
10031 continue                                                              806
      deallocate(vq)                                                        807
      return                                                                808
      end                                                                   809
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ula    812 
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                  813
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         814
      integer jd(*),ia(nx),nin(nlam)                                        815
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           820
      allocate(xm(1:ni),stat=ierr)                                          820
      jerr=jerr+ierr                                                        821
      allocate(xs(1:ni),stat=ierr)                                          821
      jerr=jerr+ierr                                                        822
      allocate(ju(1:ni),stat=ierr)                                          822
      jerr=jerr+ierr                                                        823
      allocate(xv(1:ni),stat=ierr)                                          823
      jerr=jerr+ierr                                                        824
      allocate(vlam(1:nlam),stat=ierr)                                      824
      jerr=jerr+ierr                                                        825
      if(jerr.ne.0) return                                                  826
      call chkvars(no,ni,x,ju)                                              827
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  828
      if(maxval(ju) .gt. 0)goto 10071                                       828
      jerr=7777                                                             828
      return                                                                828
10071 continue                                                              829
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)          830
      if(jerr.ne.0) return                                                  831
      cl=cl/ys                                                              831
      if(isd .le. 0)goto 10091                                              831
10100 do 10101 j=1,ni                                                       831
      cl(:,j)=cl(:,j)*xs(j)                                                 831
10101 continue                                                              831
10102 continue                                                              831
10091 continue                                                              832
      if(flmin.ge.1.0) vlam=ulam/ys                                         833
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    835 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  836
10110 do 10111 k=1,lmu                                                      836
      alm(k)=ys*alm(k)                                                      836
      nk=nin(k)                                                             837
10120 do 10121 l=1,nk                                                       837
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          837
10121 continue                                                              837
10122 continue                                                              837
      a0(k)=0.0                                                             838
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           839
10111 continue                                                              840
10112 continue                                                              840
      deallocate(xm,xs,g,ju,xv,vlam)                                        841
      return                                                                842
      end                                                                   843
      subroutine standard (no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr    844 
     *)
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  844
      integer ju(ni)                                                        845
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           848
      if(jerr.ne.0) return                                                  849
      w=w/sum(w)                                                            849
      v=sqrt(w)                                                             850
      if(intr .ne. 0)goto 10141                                             850
      ym=0.0                                                                850
      y=v*y                                                                 851
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         851
      y=y/ys                                                                852
10150 do 10151 j=1,ni                                                       852
      if(ju(j).eq.0)goto 10151                                              852
      xm(j)=0.0                                                             852
      x(:,j)=v*x(:,j)                                                       853
      xv(j)=dot_product(x(:,j),x(:,j))                                      854
      if(isd .eq. 0)goto 10171                                              854
      xbq=dot_product(v,x(:,j))**2                                          854
      vc=xv(j)-xbq                                                          855
      xs(j)=sqrt(vc)                                                        855
      x(:,j)=x(:,j)/xs(j)                                                   855
      xv(j)=1.0+xbq/vc                                                      856
      goto 10181                                                            857
10171 continue                                                              857
      xs(j)=1.0                                                             857
10181 continue                                                              858
10161 continue                                                              858
10151 continue                                                              859
10152 continue                                                              859
      goto 10191                                                            860
10141 continue                                                              861
10200 do 10201 j=1,ni                                                       861
      if(ju(j).eq.0)goto 10201                                              862
      xm(j)=dot_product(w,x(:,j))                                           862
      x(:,j)=v*(x(:,j)-xm(j))                                               863
      xv(j)=dot_product(x(:,j),x(:,j))                                      863
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        864
10201 continue                                                              865
10202 continue                                                              865
      if(isd .ne. 0)goto 10221                                              865
      xs=1.0                                                                865
      goto 10231                                                            866
10221 continue                                                              867
10240 do 10241 j=1,ni                                                       867
      if(ju(j).eq.0)goto 10241                                              867
      x(:,j)=x(:,j)/xs(j)                                                   867
10241 continue                                                              868
10242 continue                                                              868
      xv=1.0                                                                869
10231 continue                                                              870
10211 continue                                                              870
      ym=dot_product(w,y)                                                   870
      y=v*(y-ym)                                                            870
      ys=sqrt(dot_product(y,y))                                             870
      y=y/ys                                                                871
10191 continue                                                              872
10131 continue                                                              872
      g=0.0                                                                 872
10250 do 10251 j=1,ni                                                       872
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             872
10251 continue                                                              873
10252 continue                                                              873
      deallocate(v)                                                         874
      return                                                                875
      end                                                                   876
      subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,t    878 
     *hr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    879 
     *nlam),xv(ni)
      real cl(2,ni)                                                         880
      integer ju(ni),ia(nx),kin(nlam)                                       881
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)                888
      allocate(a(1:ni),stat=ierr)                                           888
      jerr=jerr+ierr                                                        889
      allocate(mm(1:ni),stat=ierr)                                          889
      jerr=jerr+ierr                                                        890
      allocate(da(1:ni),stat=ierr)                                          890
      jerr=jerr+ierr                                                        891
      if(jerr.ne.0) return                                                  892
      bta=beta                                                              892
      omb=1.0-bta                                                           893
      if(flmin .ge. 1.0)goto 10271                                          893
      eqs=max(eps,flmin)                                                    893
      alf=eqs**(1.0/(nlam-1))                                               893
10271 continue                                                              894
      rsq=0.0                                                               894
      a=0.0                                                                 894
      mm=0                                                                  894
      nlp=0                                                                 894
      nin=nlp                                                               894
      iz=0                                                                  894
      mnl=min(mnlam,nlam)                                                   895
10280 do 10281 m=1,nlam                                                     896
      if(flmin .lt. 1.0)goto 10301                                          896
      alm=ulam(m)                                                           896
      goto 10291                                                            897
10301 if(m .le. 2)goto 10311                                                897
      alm=alm*alf                                                           897
      goto 10291                                                            898
10311 if(m .ne. 1)goto 10321                                                898
      alm=big                                                               898
      goto 10331                                                            899
10321 continue                                                              899
      alm=0.0                                                               900
10340 do 10341 j=1,ni                                                       900
      if(ju(j).eq.0)goto 10341                                              900
      if(vp(j).le.0.0)goto 10341                                            901
      alm=max(alm,abs(g(j))/vp(j))                                          902
10341 continue                                                              903
10342 continue                                                              903
      alm=alf*alm/max(bta,1.0e-3)                                           904
10331 continue                                                              905
10291 continue                                                              905
      dem=alm*omb                                                           905
      ab=alm*bta                                                            905
      rsq0=rsq                                                              905
      jz=1                                                                  906
10350 continue                                                              906
10351 continue                                                              906
      if(iz*jz.ne.0) go to 10360                                            906
      nlp=nlp+1                                                             906
      dlx=0.0                                                               907
10370 do 10371 k=1,ni                                                       907
      if(ju(k).eq.0)goto 10371                                              908
      ak=a(k)                                                               908
      u=g(k)+ak*xv(k)                                                       908
      v=abs(u)-vp(k)*ab                                                     908
      a(k)=0.0                                                              910
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    911 
     *em)))
      if(a(k).eq.ak)goto 10371                                              912
      if(mm(k) .ne. 0)goto 10391                                            912
      nin=nin+1                                                             912
      if(nin.gt.nx)goto 10372                                               913
10400 do 10401 j=1,ni                                                       913
      if(ju(j).eq.0)goto 10401                                              914
      if(mm(j) .eq. 0)goto 10421                                            914
      c(j,nin)=c(k,mm(j))                                                   914
      goto 10401                                                            914
10421 continue                                                              915
      if(j .ne. k)goto 10441                                                915
      c(j,nin)=xv(j)                                                        915
      goto 10401                                                            915
10441 continue                                                              916
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   917
10401 continue                                                              918
10402 continue                                                              918
      mm(k)=nin                                                             918
      ia(nin)=k                                                             919
10391 continue                                                              920
      del=a(k)-ak                                                           920
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      921
      dlx=max(xv(k)*del**2,dlx)                                             922
10450 do 10451 j=1,ni                                                       922
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               922
10451 continue                                                              923
10452 continue                                                              923
10371 continue                                                              924
10372 continue                                                              924
      if(dlx.lt.thr)goto 10352                                              924
      if(nin.gt.nx)goto 10352                                               925
      if(nlp .le. maxit)goto 10471                                          925
      jerr=-m                                                               925
      return                                                                925
10471 continue                                                              926
10360 continue                                                              926
      iz=1                                                                  926
      da(1:nin)=a(ia(1:nin))                                                927
10480 continue                                                              927
10481 continue                                                              927
      nlp=nlp+1                                                             927
      dlx=0.0                                                               928
10490 do 10491 l=1,nin                                                      928
      k=ia(l)                                                               928
      ak=a(k)                                                               928
      u=g(k)+ak*xv(k)                                                       928
      v=abs(u)-vp(k)*ab                                                     929
      a(k)=0.0                                                              931
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    932 
     *em)))
      if(a(k).eq.ak)goto 10491                                              933
      del=a(k)-ak                                                           933
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      934
      dlx=max(xv(k)*del**2,dlx)                                             935
10500 do 10501 j=1,nin                                                      935
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  935
10501 continue                                                              936
10502 continue                                                              936
10491 continue                                                              937
10492 continue                                                              937
      if(dlx.lt.thr)goto 10482                                              937
      if(nlp .le. maxit)goto 10521                                          937
      jerr=-m                                                               937
      return                                                                937
10521 continue                                                              938
      goto 10481                                                            939
10482 continue                                                              939
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      940
10530 do 10531 j=1,ni                                                       940
      if(mm(j).ne.0)goto 10531                                              941
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            942
10531 continue                                                              943
10532 continue                                                              943
      jz=0                                                                  944
      goto 10351                                                            945
10352 continue                                                              945
      if(nin .le. nx)goto 10551                                             945
      jerr=-10000-m                                                         945
      goto 10282                                                            945
10551 continue                                                              946
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 946
      kin(m)=nin                                                            947
      rsqo(m)=rsq                                                           947
      almo(m)=alm                                                           947
      lmu=m                                                                 948
      if(m.lt.mnl)goto 10281                                                948
      if(flmin.ge.1.0)goto 10281                                            949
      me=0                                                                  949
10560 do 10561 j=1,nin                                                      949
      if(ao(j,m).ne.0.0) me=me+1                                            949
10561 continue                                                              949
10562 continue                                                              949
      if(me.gt.ne)goto 10282                                                950
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                     950
      if(rsq.gt.rsqmax)goto 10282                                           951
10281 continue                                                              952
10282 continue                                                              952
      deallocate(a,mm,c,da)                                                 953
      return                                                                954
      end                                                                   955
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam    957 
     *,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)                  958
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         959
      integer jd(*),ia(nx),nin(nlam)                                        960
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          965
      allocate(xs(1:ni),stat=ierr)                                          965
      jerr=jerr+ierr                                                        966
      allocate(ju(1:ni),stat=ierr)                                          966
      jerr=jerr+ierr                                                        967
      allocate(xv(1:ni),stat=ierr)                                          967
      jerr=jerr+ierr                                                        968
      allocate(vlam(1:nlam),stat=ierr)                                      968
      jerr=jerr+ierr                                                        969
      if(jerr.ne.0) return                                                  970
      call chkvars(no,ni,x,ju)                                              971
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  972
      if(maxval(ju) .gt. 0)goto 10581                                       972
      jerr=7777                                                             972
      return                                                                972
10581 continue                                                              973
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)           974
      if(jerr.ne.0) return                                                  975
      cl=cl/ys                                                              975
      if(isd .le. 0)goto 10601                                              975
10610 do 10611 j=1,ni                                                       975
      cl(:,j)=cl(:,j)*xs(j)                                                 975
10611 continue                                                              975
10612 continue                                                              975
10601 continue                                                              976
      if(flmin.ge.1.0) vlam=ulam/ys                                         977
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    979 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  980
10620 do 10621 k=1,lmu                                                      980
      alm(k)=ys*alm(k)                                                      980
      nk=nin(k)                                                             981
10630 do 10631 l=1,nk                                                       981
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          981
10631 continue                                                              981
10632 continue                                                              981
      a0(k)=0.0                                                             982
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           983
10621 continue                                                              984
10622 continue                                                              984
      deallocate(xm,xs,ju,xv,vlam)                                          985
      return                                                                986
      end                                                                   987
      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)    988
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        988
      integer ju(ni)                                                        989
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           992
      if(jerr.ne.0) return                                                  993
      w=w/sum(w)                                                            993
      v=sqrt(w)                                                             994
      if(intr .ne. 0)goto 10651                                             994
      ym=0.0                                                                994
      y=v*y                                                                 995
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         995
      y=y/ys                                                                996
10660 do 10661 j=1,ni                                                       996
      if(ju(j).eq.0)goto 10661                                              996
      xm(j)=0.0                                                             996
      x(:,j)=v*x(:,j)                                                       997
      xv(j)=dot_product(x(:,j),x(:,j))                                      998
      if(isd .eq. 0)goto 10681                                              998
      xbq=dot_product(v,x(:,j))**2                                          998
      vc=xv(j)-xbq                                                          999
      xs(j)=sqrt(vc)                                                        999
      x(:,j)=x(:,j)/xs(j)                                                   999
      xv(j)=1.0+xbq/vc                                                     1000
      goto 10691                                                           1001
10681 continue                                                             1001
      xs(j)=1.0                                                            1001
10691 continue                                                             1002
10671 continue                                                             1002
10661 continue                                                             1003
10662 continue                                                             1003
      go to 10700                                                          1004
10651 continue                                                             1005
10710 do 10711 j=1,ni                                                      1005
      if(ju(j).eq.0)goto 10711                                             1006
      xm(j)=dot_product(w,x(:,j))                                          1006
      x(:,j)=v*(x(:,j)-xm(j))                                              1007
      xv(j)=dot_product(x(:,j),x(:,j))                                     1007
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1008
10711 continue                                                             1009
10712 continue                                                             1009
      if(isd .ne. 0)goto 10731                                             1009
      xs=1.0                                                               1009
      goto 10741                                                           1010
10731 continue                                                             1010
10750 do 10751 j=1,ni                                                      1010
      if(ju(j).eq.0)goto 10751                                             1010
      x(:,j)=x(:,j)/xs(j)                                                  1010
10751 continue                                                             1011
10752 continue                                                             1011
      xv=1.0                                                               1012
10741 continue                                                             1013
10721 continue                                                             1013
      ym=dot_product(w,y)                                                  1013
      y=v*(y-ym)                                                           1013
      ys=sqrt(dot_product(y,y))                                            1013
      y=y/ys                                                               1014
10700 continue                                                             1014
      deallocate(v)                                                        1015
      return                                                               1016
      end                                                                  1017
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th   1019 
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(   1020 
     *nlam),xv(ni)
      real cl(2,ni)                                                        1021
      integer ju(ni),ia(nx),kin(nlam)                                      1022
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1027
      allocate(a(1:ni),stat=jerr)                                          1028
      allocate(mm(1:ni),stat=ierr)                                         1028
      jerr=jerr+ierr                                                       1029
      allocate(g(1:ni),stat=ierr)                                          1029
      jerr=jerr+ierr                                                       1030
      allocate(ix(1:ni),stat=ierr)                                         1030
      jerr=jerr+ierr                                                       1031
      if(jerr.ne.0) return                                                 1032
      bta=beta                                                             1032
      omb=1.0-bta                                                          1032
      ix=0                                                                 1033
      if(flmin .ge. 1.0)goto 10771                                         1033
      eqs=max(eps,flmin)                                                   1033
      alf=eqs**(1.0/(nlam-1))                                              1033
10771 continue                                                             1034
      rsq=0.0                                                              1034
      a=0.0                                                                1034
      mm=0                                                                 1034
      nlp=0                                                                1034
      nin=nlp                                                              1034
      iz=0                                                                 1034
      mnl=min(mnlam,nlam)                                                  1034
      alm=0.0                                                              1035
10780 do 10781 j=1,ni                                                      1035
      if(ju(j).eq.0)goto 10781                                             1035
      g(j)=abs(dot_product(y,x(:,j)))                                      1035
10781 continue                                                             1036
10782 continue                                                             1036
10790 do 10791 m=1,nlam                                                    1036
      alm0=alm                                                             1037
      if(flmin .lt. 1.0)goto 10811                                         1037
      alm=ulam(m)                                                          1037
      goto 10801                                                           1038
10811 if(m .le. 2)goto 10821                                               1038
      alm=alm*alf                                                          1038
      goto 10801                                                           1039
10821 if(m .ne. 1)goto 10831                                               1039
      alm=big                                                              1039
      goto 10841                                                           1040
10831 continue                                                             1040
      alm0=0.0                                                             1041
10850 do 10851 j=1,ni                                                      1041
      if(ju(j).eq.0)goto 10851                                             1041
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1041
10851 continue                                                             1042
10852 continue                                                             1042
      alm0=alm0/max(bta,1.0e-3)                                            1042
      alm=alf*alm0                                                         1043
10841 continue                                                             1044
10801 continue                                                             1044
      dem=alm*omb                                                          1044
      ab=alm*bta                                                           1044
      rsq0=rsq                                                             1044
      jz=1                                                                 1045
      tlam=bta*(2.0*alm-alm0)                                              1046
10860 do 10861 k=1,ni                                                      1046
      if(ix(k).eq.1)goto 10861                                             1046
      if(ju(k).eq.0)goto 10861                                             1047
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1048
10861 continue                                                             1049
10862 continue                                                             1049
10870 continue                                                             1049
10871 continue                                                             1049
      if(iz*jz.ne.0) go to 10360                                           1050
10880 continue                                                             1050
      nlp=nlp+1                                                            1050
      dlx=0.0                                                              1051
10890 do 10891 k=1,ni                                                      1051
      if(ix(k).eq.0)goto 10891                                             1051
      gk=dot_product(y,x(:,k))                                             1052
      ak=a(k)                                                              1052
      u=gk+ak*xv(k)                                                        1052
      v=abs(u)-vp(k)*ab                                                    1052
      a(k)=0.0                                                             1054
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1055 
     *em)))
      if(a(k).eq.ak)goto 10891                                             1056
      if(mm(k) .ne. 0)goto 10911                                           1056
      nin=nin+1                                                            1056
      if(nin.gt.nx)goto 10892                                              1057
      mm(k)=nin                                                            1057
      ia(nin)=k                                                            1058
10911 continue                                                             1059
      del=a(k)-ak                                                          1059
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1060
      y=y-del*x(:,k)                                                       1060
      dlx=max(xv(k)*del**2,dlx)                                            1061
10891 continue                                                             1062
10892 continue                                                             1062
      if(nin.gt.nx)goto 10872                                              1063
      if(dlx .ge. thr)goto 10931                                           1063
      ixx=0                                                                1064
10940 do 10941 k=1,ni                                                      1064
      if(ix(k).eq.1)goto 10941                                             1064
      if(ju(k).eq.0)goto 10941                                             1065
      g(k)=abs(dot_product(y,x(:,k)))                                      1066
      if(g(k) .le. ab*vp(k))goto 10961                                     1066
      ix(k)=1                                                              1066
      ixx=1                                                                1066
10961 continue                                                             1067
10941 continue                                                             1068
10942 continue                                                             1068
      if(ixx.eq.1) go to 10880                                             1069
      goto 10872                                                           1070
10931 continue                                                             1071
      if(nlp .le. maxit)goto 10981                                         1071
      jerr=-m                                                              1071
      return                                                               1071
10981 continue                                                             1072
10360 continue                                                             1072
      iz=1                                                                 1073
10990 continue                                                             1073
10991 continue                                                             1073
      nlp=nlp+1                                                            1073
      dlx=0.0                                                              1074
11000 do 11001 l=1,nin                                                     1074
      k=ia(l)                                                              1074
      gk=dot_product(y,x(:,k))                                             1075
      ak=a(k)                                                              1075
      u=gk+ak*xv(k)                                                        1075
      v=abs(u)-vp(k)*ab                                                    1075
      a(k)=0.0                                                             1077
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1078 
     *em)))
      if(a(k).eq.ak)goto 11001                                             1079
      del=a(k)-ak                                                          1079
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1080
      y=y-del*x(:,k)                                                       1080
      dlx=max(xv(k)*del**2,dlx)                                            1081
11001 continue                                                             1082
11002 continue                                                             1082
      if(dlx.lt.thr)goto 10992                                             1082
      if(nlp .le. maxit)goto 11021                                         1082
      jerr=-m                                                              1082
      return                                                               1082
11021 continue                                                             1083
      goto 10991                                                           1084
10992 continue                                                             1084
      jz=0                                                                 1085
      goto 10871                                                           1086
10872 continue                                                             1086
      if(nin .le. nx)goto 11041                                            1086
      jerr=-10000-m                                                        1086
      goto 10792                                                           1086
11041 continue                                                             1087
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1087
      kin(m)=nin                                                           1088
      rsqo(m)=rsq                                                          1088
      almo(m)=alm                                                          1088
      lmu=m                                                                1089
      if(m.lt.mnl)goto 10791                                               1089
      if(flmin.ge.1.0)goto 10791                                           1090
      me=0                                                                 1090
11050 do 11051 j=1,nin                                                     1090
      if(ao(j,m).ne.0.0) me=me+1                                           1090
11051 continue                                                             1090
11052 continue                                                             1090
      if(me.gt.ne)goto 10792                                               1091
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                    1091
      if(rsq.gt.rsqmax)goto 10792                                          1092
10791 continue                                                             1093
10792 continue                                                             1093
      deallocate(a,mm,g,ix)                                                1094
      return                                                               1095
      end                                                                  1096
      subroutine chkvars(no,ni,x,ju)                                       1097
      real x(no,ni)                                                        1097
      integer ju(ni)                                                       1098
11060 do 11061 j=1,ni                                                      1098
      ju(j)=0                                                              1098
      t=x(1,j)                                                             1099
11070 do 11071 i=2,no                                                      1099
      if(x(i,j).eq.t)goto 11071                                            1099
      ju(j)=1                                                              1099
      goto 11072                                                           1099
11071 continue                                                             1100
11072 continue                                                             1100
11061 continue                                                             1101
11062 continue                                                             1101
      return                                                               1102
      end                                                                  1103
      subroutine uncomp(ni,ca,ia,nin,a)                                    1104
      real ca(*),a(ni)                                                     1104
      integer ia(*)                                                        1105
      a=0.0                                                                1105
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1106
      return                                                               1107
      end                                                                  1108
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1109
      real ca(nin),x(n,*),f(n)                                             1109
      integer ia(nin)                                                      1110
      f=a0                                                                 1110
      if(nin.le.0) return                                                  1111
11080 do 11081 i=1,n                                                       1111
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1111
11081 continue                                                             1112
11082 continue                                                             1112
      return                                                               1113
      end                                                                  1114
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   1117 
     *,flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                     1118
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1119
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1120
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 11101                                    1123
      jerr=10000                                                           1123
      return                                                               1123
11101 continue                                                             1124
      allocate(vq(1:ni),stat=jerr)                                         1124
      if(jerr.ne.0) return                                                 1125
      vq=max(0.0,vp)                                                       1125
      vq=vq*ni/sum(vq)                                                     1126
      if(ka .ne. 1)goto 11121                                              1127
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   1130 
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11131                                                           1131
11121 continue                                                             1132
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   1135 
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11131 continue                                                             1136
11111 continue                                                             1136
      deallocate(vq)                                                       1137
      return                                                               1138
      end                                                                  1139
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   1142 
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                     1143
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1144
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1145
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1150
      allocate(xm(1:ni),stat=ierr)                                         1150
      jerr=jerr+ierr                                                       1151
      allocate(xs(1:ni),stat=ierr)                                         1151
      jerr=jerr+ierr                                                       1152
      allocate(ju(1:ni),stat=ierr)                                         1152
      jerr=jerr+ierr                                                       1153
      allocate(xv(1:ni),stat=ierr)                                         1153
      jerr=jerr+ierr                                                       1154
      allocate(vlam(1:nlam),stat=ierr)                                     1154
      jerr=jerr+ierr                                                       1155
      if(jerr.ne.0) return                                                 1156
      call spchkvars(no,ni,x,ix,ju)                                        1157
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1158
      if(maxval(ju) .gt. 0)goto 11151                                      1158
      jerr=7777                                                            1158
      return                                                               1158
11151 continue                                                             1159
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer   1160 
     *r)
      if(jerr.ne.0) return                                                 1161
      cl=cl/ys                                                             1161
      if(isd .le. 0)goto 11171                                             1161
11180 do 11181 j=1,ni                                                      1161
      cl(:,j)=cl(:,j)*xs(j)                                                1161
11181 continue                                                             1161
11182 continue                                                             1161
11171 continue                                                             1162
      if(flmin.ge.1.0) vlam=ulam/ys                                        1163
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1165 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1166
11190 do 11191 k=1,lmu                                                     1166
      alm(k)=ys*alm(k)                                                     1166
      nk=nin(k)                                                            1167
11200 do 11201 l=1,nk                                                      1167
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1167
11201 continue                                                             1167
11202 continue                                                             1167
      a0(k)=0.0                                                            1168
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1169
11191 continue                                                             1170
11192 continue                                                             1170
      deallocate(xm,xs,g,ju,xv,vlam)                                       1171
      return                                                               1172
      end                                                                  1173
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys   1174 
     *,xv,jerr)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                     1174
      integer ix(*),jx(*),ju(ni)                                           1175
      w=w/sum(w)                                                           1176
      if(intr .ne. 0)goto 11221                                            1176
      ym=0.0                                                               1177
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1177
      y=y/ys                                                               1178
11230 do 11231 j=1,ni                                                      1178
      if(ju(j).eq.0)goto 11231                                             1178
      xm(j)=0.0                                                            1178
      jb=ix(j)                                                             1178
      je=ix(j+1)-1                                                         1179
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1180
      if(isd .eq. 0)goto 11251                                             1180
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1180
      vc=xv(j)-xbq                                                         1181
      xs(j)=sqrt(vc)                                                       1181
      xv(j)=1.0+xbq/vc                                                     1182
      goto 11261                                                           1183
11251 continue                                                             1183
      xs(j)=1.0                                                            1183
11261 continue                                                             1184
11241 continue                                                             1184
11231 continue                                                             1185
11232 continue                                                             1185
      goto 11271                                                           1186
11221 continue                                                             1187
11280 do 11281 j=1,ni                                                      1187
      if(ju(j).eq.0)goto 11281                                             1188
      jb=ix(j)                                                             1188
      je=ix(j+1)-1                                                         1188
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1189
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1190
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1191
11281 continue                                                             1192
11282 continue                                                             1192
      if(isd .ne. 0)goto 11301                                             1192
      xs=1.0                                                               1192
      goto 11311                                                           1192
11301 continue                                                             1192
      xv=1.0                                                               1192
11311 continue                                                             1193
11291 continue                                                             1193
      ym=dot_product(w,y)                                                  1193
      y=y-ym                                                               1193
      ys=sqrt(dot_product(w,y**2))                                         1193
      y=y/ys                                                               1194
11271 continue                                                             1195
11211 continue                                                             1195
      g=0.0                                                                1196
11320 do 11321 j=1,ni                                                      1196
      if(ju(j).eq.0)goto 11321                                             1196
      jb=ix(j)                                                             1196
      je=ix(j+1)-1                                                         1197
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1198
11321 continue                                                             1199
11322 continue                                                             1199
      return                                                               1200
      end                                                                  1201
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1203 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                              1204
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni),cl(2,n   1205 
     *i)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1206
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1213
      allocate(a(1:ni),stat=ierr)                                          1213
      jerr=jerr+ierr                                                       1214
      allocate(mm(1:ni),stat=ierr)                                         1214
      jerr=jerr+ierr                                                       1215
      allocate(da(1:ni),stat=ierr)                                         1215
      jerr=jerr+ierr                                                       1216
      if(jerr.ne.0) return                                                 1217
      bta=beta                                                             1217
      omb=1.0-bta                                                          1218
      if(flmin .ge. 1.0)goto 11341                                         1218
      eqs=max(eps,flmin)                                                   1218
      alf=eqs**(1.0/(nlam-1))                                              1218
11341 continue                                                             1219
      rsq=0.0                                                              1219
      a=0.0                                                                1219
      mm=0                                                                 1219
      nlp=0                                                                1219
      nin=nlp                                                              1219
      iz=0                                                                 1219
      mnl=min(mnlam,nlam)                                                  1220
11350 do 11351 m=1,nlam                                                    1221
      if(flmin .lt. 1.0)goto 11371                                         1221
      alm=ulam(m)                                                          1221
      goto 11361                                                           1222
11371 if(m .le. 2)goto 11381                                               1222
      alm=alm*alf                                                          1222
      goto 11361                                                           1223
11381 if(m .ne. 1)goto 11391                                               1223
      alm=big                                                              1223
      goto 11401                                                           1224
11391 continue                                                             1224
      alm=0.0                                                              1225
11410 do 11411 j=1,ni                                                      1225
      if(ju(j).eq.0)goto 11411                                             1225
      if(vp(j).le.0.0)goto 11411                                           1226
      alm=max(alm,abs(g(j))/vp(j))                                         1227
11411 continue                                                             1228
11412 continue                                                             1228
      alm=alf*alm/max(bta,1.0e-3)                                          1229
11401 continue                                                             1230
11361 continue                                                             1230
      dem=alm*omb                                                          1230
      ab=alm*bta                                                           1230
      rsq0=rsq                                                             1230
      jz=1                                                                 1231
11420 continue                                                             1231
11421 continue                                                             1231
      if(iz*jz.ne.0) go to 10360                                           1231
      nlp=nlp+1                                                            1231
      dlx=0.0                                                              1232
11430 do 11431 k=1,ni                                                      1232
      if(ju(k).eq.0)goto 11431                                             1233
      ak=a(k)                                                              1233
      u=g(k)+ak*xv(k)                                                      1233
      v=abs(u)-vp(k)*ab                                                    1233
      a(k)=0.0                                                             1235
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1236 
     *em)))
      if(a(k).eq.ak)goto 11431                                             1237
      if(mm(k) .ne. 0)goto 11451                                           1237
      nin=nin+1                                                            1237
      if(nin.gt.nx)goto 11432                                              1238
11460 do 11461 j=1,ni                                                      1238
      if(ju(j).eq.0)goto 11461                                             1239
      if(mm(j) .eq. 0)goto 11481                                           1239
      c(j,nin)=c(k,mm(j))                                                  1239
      goto 11461                                                           1239
11481 continue                                                             1240
      if(j .ne. k)goto 11501                                               1240
      c(j,nin)=xv(j)                                                       1240
      goto 11461                                                           1240
11501 continue                                                             1241
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1243
11461 continue                                                             1244
11462 continue                                                             1244
      mm(k)=nin                                                            1244
      ia(nin)=k                                                            1245
11451 continue                                                             1246
      del=a(k)-ak                                                          1246
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1247
      dlx=max(xv(k)*del**2,dlx)                                            1248
11510 do 11511 j=1,ni                                                      1248
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1248
11511 continue                                                             1249
11512 continue                                                             1249
11431 continue                                                             1250
11432 continue                                                             1250
      if(dlx.lt.thr)goto 11422                                             1250
      if(nin.gt.nx)goto 11422                                              1251
      if(nlp .le. maxit)goto 11531                                         1251
      jerr=-m                                                              1251
      return                                                               1251
11531 continue                                                             1252
10360 continue                                                             1252
      iz=1                                                                 1252
      da(1:nin)=a(ia(1:nin))                                               1253
11540 continue                                                             1253
11541 continue                                                             1253
      nlp=nlp+1                                                            1253
      dlx=0.0                                                              1254
11550 do 11551 l=1,nin                                                     1254
      k=ia(l)                                                              1255
      ak=a(k)                                                              1255
      u=g(k)+ak*xv(k)                                                      1255
      v=abs(u)-vp(k)*ab                                                    1255
      a(k)=0.0                                                             1257
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1258 
     *em)))
      if(a(k).eq.ak)goto 11551                                             1259
      del=a(k)-ak                                                          1259
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1260
      dlx=max(xv(k)*del**2,dlx)                                            1261
11560 do 11561 j=1,nin                                                     1261
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1261
11561 continue                                                             1262
11562 continue                                                             1262
11551 continue                                                             1263
11552 continue                                                             1263
      if(dlx.lt.thr)goto 11542                                             1263
      if(nlp .le. maxit)goto 11581                                         1263
      jerr=-m                                                              1263
      return                                                               1263
11581 continue                                                             1264
      goto 11541                                                           1265
11542 continue                                                             1265
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1266
11590 do 11591 j=1,ni                                                      1266
      if(mm(j).ne.0)goto 11591                                             1267
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1268
11591 continue                                                             1269
11592 continue                                                             1269
      jz=0                                                                 1270
      goto 11421                                                           1271
11422 continue                                                             1271
      if(nin .le. nx)goto 11611                                            1271
      jerr=-10000-m                                                        1271
      goto 11352                                                           1271
11611 continue                                                             1272
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1272
      kin(m)=nin                                                           1273
      rsqo(m)=rsq                                                          1273
      almo(m)=alm                                                          1273
      lmu=m                                                                1274
      if(m.lt.mnl)goto 11351                                               1274
      if(flmin.ge.1.0)goto 11351                                           1275
      me=0                                                                 1275
11620 do 11621 j=1,nin                                                     1275
      if(ao(j,m).ne.0.0) me=me+1                                           1275
11621 continue                                                             1275
11622 continue                                                             1275
      if(me.gt.ne)goto 11352                                               1276
      if(rsq-rsq0.lt.sml*rsq)goto 11352                                    1276
      if(rsq.gt.rsqmax)goto 11352                                          1277
11351 continue                                                             1278
11352 continue                                                             1278
      deallocate(a,mm,c,da)                                                1279
      return                                                               1280
      end                                                                  1281
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   1283 
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)                     1284
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1285
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1286
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1291
      allocate(xs(1:ni),stat=ierr)                                         1291
      jerr=jerr+ierr                                                       1292
      allocate(ju(1:ni),stat=ierr)                                         1292
      jerr=jerr+ierr                                                       1293
      allocate(xv(1:ni),stat=ierr)                                         1293
      jerr=jerr+ierr                                                       1294
      allocate(vlam(1:nlam),stat=ierr)                                     1294
      jerr=jerr+ierr                                                       1295
      if(jerr.ne.0) return                                                 1296
      call spchkvars(no,ni,x,ix,ju)                                        1297
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1298
      if(maxval(ju) .gt. 0)goto 11641                                      1298
      jerr=7777                                                            1298
      return                                                               1298
11641 continue                                                             1299
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr   1300 
     *)
      if(jerr.ne.0) return                                                 1301
      cl=cl/ys                                                             1301
      if(isd .le. 0)goto 11661                                             1301
11670 do 11671 j=1,ni                                                      1301
      cl(:,j)=cl(:,j)*xs(j)                                                1301
11671 continue                                                             1301
11672 continue                                                             1301
11661 continue                                                             1302
      if(flmin.ge.1.0) vlam=ulam/ys                                        1303
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1305 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1306
11680 do 11681 k=1,lmu                                                     1306
      alm(k)=ys*alm(k)                                                     1306
      nk=nin(k)                                                            1307
11690 do 11691 l=1,nk                                                      1307
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1307
11691 continue                                                             1307
11692 continue                                                             1307
      a0(k)=0.0                                                            1308
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1309
11681 continue                                                             1310
11682 continue                                                             1310
      deallocate(xm,xs,ju,xv,vlam)                                         1311
      return                                                               1312
      end                                                                  1313
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,   1314 
     *xv,jerr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1314
      integer ix(*),jx(*),ju(ni)                                           1315
      w=w/sum(w)                                                           1316
      if(intr .ne. 0)goto 11711                                            1316
      ym=0.0                                                               1317
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1317
      y=y/ys                                                               1318
11720 do 11721 j=1,ni                                                      1318
      if(ju(j).eq.0)goto 11721                                             1318
      xm(j)=0.0                                                            1318
      jb=ix(j)                                                             1318
      je=ix(j+1)-1                                                         1319
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1320
      if(isd .eq. 0)goto 11741                                             1320
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1320
      vc=xv(j)-xbq                                                         1321
      xs(j)=sqrt(vc)                                                       1321
      xv(j)=1.0+xbq/vc                                                     1322
      goto 11751                                                           1323
11741 continue                                                             1323
      xs(j)=1.0                                                            1323
11751 continue                                                             1324
11731 continue                                                             1324
11721 continue                                                             1325
11722 continue                                                             1325
      return                                                               1326
11711 continue                                                             1327
11760 do 11761 j=1,ni                                                      1327
      if(ju(j).eq.0)goto 11761                                             1328
      jb=ix(j)                                                             1328
      je=ix(j+1)-1                                                         1328
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1329
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1330
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1331
11761 continue                                                             1332
11762 continue                                                             1332
      if(isd .ne. 0)goto 11781                                             1332
      xs=1.0                                                               1332
      goto 11791                                                           1332
11781 continue                                                             1332
      xv=1.0                                                               1332
11791 continue                                                             1333
11771 continue                                                             1333
      ym=dot_product(w,y)                                                  1333
      y=y-ym                                                               1333
      ys=sqrt(dot_product(w,y**2))                                         1333
      y=y/ys                                                               1334
      return                                                               1335
      end                                                                  1336
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1338 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)                     1339
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1340
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1341
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1346
      allocate(a(1:ni),stat=jerr)                                          1347
      allocate(mm(1:ni),stat=ierr)                                         1347
      jerr=jerr+ierr                                                       1348
      allocate(g(1:ni),stat=ierr)                                          1348
      jerr=jerr+ierr                                                       1349
      allocate(iy(1:ni),stat=ierr)                                         1349
      jerr=jerr+ierr                                                       1350
      if(jerr.ne.0) return                                                 1351
      bta=beta                                                             1351
      omb=1.0-bta                                                          1351
      alm=0.0                                                              1351
      iy=0                                                                 1352
      if(flmin .ge. 1.0)goto 11811                                         1352
      eqs=max(eps,flmin)                                                   1352
      alf=eqs**(1.0/(nlam-1))                                              1352
11811 continue                                                             1353
      rsq=0.0                                                              1353
      a=0.0                                                                1353
      mm=0                                                                 1353
      o=0.0                                                                1353
      nlp=0                                                                1353
      nin=nlp                                                              1353
      iz=0                                                                 1353
      mnl=min(mnlam,nlam)                                                  1354
11820 do 11821 j=1,ni                                                      1354
      if(ju(j).eq.0)goto 11821                                             1355
      jb=ix(j)                                                             1355
      je=ix(j+1)-1                                                         1356
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1357
11821 continue                                                             1358
11822 continue                                                             1358
11830 do 11831 m=1,nlam                                                    1358
      alm0=alm                                                             1359
      if(flmin .lt. 1.0)goto 11851                                         1359
      alm=ulam(m)                                                          1359
      goto 11841                                                           1360
11851 if(m .le. 2)goto 11861                                               1360
      alm=alm*alf                                                          1360
      goto 11841                                                           1361
11861 if(m .ne. 1)goto 11871                                               1361
      alm=big                                                              1361
      goto 11881                                                           1362
11871 continue                                                             1362
      alm0=0.0                                                             1363
11890 do 11891 j=1,ni                                                      1363
      if(ju(j).eq.0)goto 11891                                             1363
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1363
11891 continue                                                             1364
11892 continue                                                             1364
      alm0=alm0/max(bta,1.0e-3)                                            1364
      alm=alf*alm0                                                         1365
11881 continue                                                             1366
11841 continue                                                             1366
      dem=alm*omb                                                          1366
      ab=alm*bta                                                           1366
      rsq0=rsq                                                             1366
      jz=1                                                                 1367
      tlam=bta*(2.0*alm-alm0)                                              1368
11900 do 11901 k=1,ni                                                      1368
      if(iy(k).eq.1)goto 11901                                             1368
      if(ju(k).eq.0)goto 11901                                             1369
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1370
11901 continue                                                             1371
11902 continue                                                             1371
11910 continue                                                             1371
11911 continue                                                             1371
      if(iz*jz.ne.0) go to 10360                                           1372
10880 continue                                                             1372
      nlp=nlp+1                                                            1372
      dlx=0.0                                                              1373
11920 do 11921 k=1,ni                                                      1373
      if(iy(k).eq.0)goto 11921                                             1373
      jb=ix(k)                                                             1373
      je=ix(k+1)-1                                                         1374
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1375
      ak=a(k)                                                              1375
      u=gk+ak*xv(k)                                                        1375
      v=abs(u)-vp(k)*ab                                                    1375
      a(k)=0.0                                                             1377
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1378 
     *em)))
      if(a(k).eq.ak)goto 11921                                             1379
      if(mm(k) .ne. 0)goto 11941                                           1379
      nin=nin+1                                                            1379
      if(nin.gt.nx)goto 11922                                              1380
      mm(k)=nin                                                            1380
      ia(nin)=k                                                            1381
11941 continue                                                             1382
      del=a(k)-ak                                                          1382
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1383
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1384
      o=o+del*xm(k)/xs(k)                                                  1384
      dlx=max(xv(k)*del**2,dlx)                                            1385
11921 continue                                                             1386
11922 continue                                                             1386
      if(nin.gt.nx)goto 11912                                              1387
      if(dlx .ge. thr)goto 11961                                           1387
      ixx=0                                                                1388
11970 do 11971 j=1,ni                                                      1388
      if(iy(j).eq.1)goto 11971                                             1388
      if(ju(j).eq.0)goto 11971                                             1389
      jb=ix(j)                                                             1389
      je=ix(j+1)-1                                                         1390
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1391
      if(g(j) .le. ab*vp(j))goto 11991                                     1391
      iy(j)=1                                                              1391
      ixx=1                                                                1391
11991 continue                                                             1392
11971 continue                                                             1393
11972 continue                                                             1393
      if(ixx.eq.1) go to 10880                                             1394
      goto 11912                                                           1395
11961 continue                                                             1396
      if(nlp .le. maxit)goto 12011                                         1396
      jerr=-m                                                              1396
      return                                                               1396
12011 continue                                                             1397
10360 continue                                                             1397
      iz=1                                                                 1398
12020 continue                                                             1398
12021 continue                                                             1398
      nlp=nlp+1                                                            1398
      dlx=0.0                                                              1399
12030 do 12031 l=1,nin                                                     1399
      k=ia(l)                                                              1399
      jb=ix(k)                                                             1399
      je=ix(k+1)-1                                                         1400
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1401
      ak=a(k)                                                              1401
      u=gk+ak*xv(k)                                                        1401
      v=abs(u)-vp(k)*ab                                                    1401
      a(k)=0.0                                                             1403
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1404 
     *em)))
      if(a(k).eq.ak)goto 12031                                             1405
      del=a(k)-ak                                                          1405
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1406
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1407
      o=o+del*xm(k)/xs(k)                                                  1407
      dlx=max(xv(k)*del**2,dlx)                                            1408
12031 continue                                                             1409
12032 continue                                                             1409
      if(dlx.lt.thr)goto 12022                                             1409
      if(nlp .le. maxit)goto 12051                                         1409
      jerr=-m                                                              1409
      return                                                               1409
12051 continue                                                             1410
      goto 12021                                                           1411
12022 continue                                                             1411
      jz=0                                                                 1412
      goto 11911                                                           1413
11912 continue                                                             1413
      if(nin .le. nx)goto 12071                                            1413
      jerr=-10000-m                                                        1413
      goto 11832                                                           1413
12071 continue                                                             1414
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1414
      kin(m)=nin                                                           1415
      rsqo(m)=rsq                                                          1415
      almo(m)=alm                                                          1415
      lmu=m                                                                1416
      if(m.lt.mnl)goto 11831                                               1416
      if(flmin.ge.1.0)goto 11831                                           1417
      me=0                                                                 1417
12080 do 12081 j=1,nin                                                     1417
      if(ao(j,m).ne.0.0) me=me+1                                           1417
12081 continue                                                             1417
12082 continue                                                             1417
      if(me.gt.ne)goto 11832                                               1418
      if(rsq-rsq0.lt.sml*rsq)goto 11832                                    1418
      if(rsq.gt.rsqmax)goto 11832                                          1419
11831 continue                                                             1420
11832 continue                                                             1420
      deallocate(a,mm,g,iy)                                                1421
      return                                                               1422
      end                                                                  1423
      subroutine spchkvars(no,ni,x,ix,ju)                                  1424
      real x(*)                                                            1424
      integer ix(*),ju(ni)                                                 1425
12090 do 12091 j=1,ni                                                      1425
      ju(j)=0                                                              1425
      jb=ix(j)                                                             1425
      nj=ix(j+1)-jb                                                        1425
      if(nj.eq.0)goto 12091                                                1426
      je=ix(j+1)-1                                                         1427
      if(nj .ge. no)goto 12111                                             1427
12120 do 12121 i=jb,je                                                     1427
      if(x(i).eq.0.0)goto 12121                                            1427
      ju(j)=1                                                              1427
      goto 12122                                                           1427
12121 continue                                                             1427
12122 continue                                                             1427
      goto 12131                                                           1428
12111 continue                                                             1428
      t=x(jb)                                                              1428
12140 do 12141 i=jb+1,je                                                   1428
      if(x(i).eq.t)goto 12141                                              1428
      ju(j)=1                                                              1428
      goto 12142                                                           1428
12141 continue                                                             1428
12142 continue                                                             1428
12131 continue                                                             1429
12101 continue                                                             1429
12091 continue                                                             1430
12092 continue                                                             1430
      return                                                               1431
      end                                                                  1432
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1433
      real ca(*),x(*),f(n)                                                 1433
      integer ia(*),ix(*),jx(*)                                            1434
      f=a0                                                                 1435
12150 do 12151 j=1,nin                                                     1435
      k=ia(j)                                                              1435
      kb=ix(k)                                                             1435
      ke=ix(k+1)-1                                                         1436
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1437
12151 continue                                                             1438
12152 continue                                                             1438
      return                                                               1439
      end                                                                  1440
      function row_prod(i,j,ia,ja,ra,w)                                    1441
      integer ia(*),ja(*)                                                  1441
      real ra(*),w(*)                                                      1442
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1444 
     *i),ia(j+1)-ia(j),w)
      return                                                               1445
      end                                                                  1446
      function dot(x,y,mx,my,nx,ny,w)                                      1447
      real x(*),y(*),w(*)                                                  1447
      integer mx(*),my(*)                                                  1448
      i=1                                                                  1448
      j=i                                                                  1448
      s=0.0                                                                1449
12160 continue                                                             1449
12161 continue                                                             1449
12170 continue                                                             1450
12171 if(mx(i).ge.my(j))goto 12172                                         1450
      i=i+1                                                                1450
      if(i.gt.nx) go to 12180                                              1450
      goto 12171                                                           1451
12172 continue                                                             1451
      if(mx(i).eq.my(j)) go to 12190                                       1452
12200 continue                                                             1452
12201 if(my(j).ge.mx(i))goto 12202                                         1452
      j=j+1                                                                1452
      if(j.gt.ny) go to 12180                                              1452
      goto 12201                                                           1453
12202 continue                                                             1453
      if(mx(i).eq.my(j)) go to 12190                                       1453
      goto 12161                                                           1454
12190 continue                                                             1454
      s=s+w(mx(i))*x(i)*y(j)                                               1455
      i=i+1                                                                1455
      if(i.gt.nx)goto 12162                                                1455
      j=j+1                                                                1455
      if(j.gt.ny)goto 12162                                                1456
      goto 12161                                                           1457
12162 continue                                                             1457
12180 continue                                                             1457
      dot=s                                                                1458
      return                                                               1459
      end                                                                  1460
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1462 
     *lam,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1463
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)         1464
      integer jd(*),ia(nx),nin(nlam)                                       1465
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 12221                                    1469
      jerr=10000                                                           1469
      return                                                               1469
12221 continue                                                             1470
      allocate(ww(1:no),stat=jerr)                                         1471
      allocate(ju(1:ni),stat=ierr)                                         1471
      jerr=jerr+ierr                                                       1472
      allocate(vq(1:ni),stat=ierr)                                         1472
      jerr=jerr+ierr                                                       1473
      allocate(xm(1:ni),stat=ierr)                                         1473
      jerr=jerr+ierr                                                       1474
      if(kopt .ne. 2)goto 12241                                            1474
      allocate(xv(1:ni),stat=ierr)                                         1474
      jerr=jerr+ierr                                                       1474
12241 continue                                                             1475
      if(isd .le. 0)goto 12261                                             1475
      allocate(xs(1:ni),stat=ierr)                                         1475
      jerr=jerr+ierr                                                       1475
12261 continue                                                             1476
      if(jerr.ne.0) return                                                 1477
      call chkvars(no,ni,x,ju)                                             1478
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1479
      if(maxval(ju) .gt. 0)goto 12281                                      1479
      jerr=7777                                                            1479
      return                                                               1479
12281 continue                                                             1480
      vq=max(0.0,vp)                                                       1480
      vq=vq*ni/sum(vq)                                                     1481
12290 do 12291 i=1,no                                                      1481
      ww(i)=sum(y(i,:))                                                    1481
      y(i,:)=y(i,:)/ww(i)                                                  1481
12291 continue                                                             1481
12292 continue                                                             1481
      sw=sum(ww)                                                           1481
      ww=ww/sw                                                             1482
      if(nc .ne. 1)goto 12311                                              1482
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1483
      if(isd .le. 0)goto 12331                                             1483
12340 do 12341 j=1,ni                                                      1483
      cl(:,j)=cl(:,j)*xs(j)                                                1483
12341 continue                                                             1483
12342 continue                                                             1483
12331 continue                                                             1484
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   1486 
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12301                                                           1487
12311 if(kopt .ne. 2)goto 12351                                            1487
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)                 1488
      if(isd .le. 0)goto 12371                                             1488
12380 do 12381 j=1,ni                                                      1488
      cl(:,j)=cl(:,j)*xs(j)                                                1488
12381 continue                                                             1488
12382 continue                                                             1488
12371 continue                                                             1489
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   1491 
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12391                                                           1492
12351 continue                                                             1492
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1493
      if(isd .le. 0)goto 12411                                             1493
12420 do 12421 j=1,ni                                                      1493
      cl(:,j)=cl(:,j)*xs(j)                                                1493
12421 continue                                                             1493
12422 continue                                                             1493
12411 continue                                                             1494
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   1496 
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12391 continue                                                             1497
12301 continue                                                             1497
      if(jerr.gt.0) return                                                 1497
      dev0=2.0*sw*dev0                                                     1498
12430 do 12431 k=1,lmu                                                     1498
      nk=nin(k)                                                            1499
12440 do 12441 ic=1,nc                                                     1499
      if(isd .le. 0)goto 12461                                             1499
12470 do 12471 l=1,nk                                                      1499
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1499
12471 continue                                                             1499
12472 continue                                                             1499
12461 continue                                                             1500
      if(intr .ne. 0)goto 12491                                            1500
      a0(ic,k)=0.0                                                         1500
      goto 12501                                                           1501
12491 continue                                                             1501
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1501
12501 continue                                                             1502
12481 continue                                                             1502
12441 continue                                                             1503
12442 continue                                                             1503
12431 continue                                                             1504
12432 continue                                                             1504
      deallocate(ww,ju,vq,xm)                                              1504
      if(isd.gt.0) deallocate(xs)                                          1505
      if(kopt.eq.2) deallocate(xv)                                         1506
      return                                                               1507
      end                                                                  1508
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  1509
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1509
      integer ju(ni)                                                       1510
      if(intr .ne. 0)goto 12521                                            1511
12530 do 12531 j=1,ni                                                      1511
      if(ju(j).eq.0)goto 12531                                             1511
      xm(j)=0.0                                                            1512
      if(isd .eq. 0)goto 12551                                             1512
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 1513
      xs(j)=sqrt(vc)                                                       1513
      x(:,j)=x(:,j)/xs(j)                                                  1514
12551 continue                                                             1515
12531 continue                                                             1516
12532 continue                                                             1516
      return                                                               1517
12521 continue                                                             1518
12560 do 12561 j=1,ni                                                      1518
      if(ju(j).eq.0)goto 12561                                             1519
      xm(j)=dot_product(w,x(:,j))                                          1519
      x(:,j)=x(:,j)-xm(j)                                                  1520
      if(isd .le. 0)goto 12581                                             1520
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1520
      x(:,j)=x(:,j)/xs(j)                                                  1520
12581 continue                                                             1521
12561 continue                                                             1522
12562 continue                                                             1522
      return                                                               1523
      end                                                                  1524
      subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv)           1525
      real x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                             1525
      integer ju(ni)                                                       1526
      if(intr .ne. 0)goto 12601                                            1527
12610 do 12611 j=1,ni                                                      1527
      if(ju(j).eq.0)goto 12611                                             1527
      xm(j)=0.0                                                            1528
      xv(j)=dot_product(w,x(:,j)**2)                                       1529
      if(isd .eq. 0)goto 12631                                             1529
      xbq=dot_product(w,x(:,j))**2                                         1529
      vc=xv(j)-xbq                                                         1530
      xs(j)=sqrt(vc)                                                       1530
      x(:,j)=x(:,j)/xs(j)                                                  1530
      xv(j)=1.0+xbq/vc                                                     1531
12631 continue                                                             1532
12611 continue                                                             1533
12612 continue                                                             1533
      return                                                               1534
12601 continue                                                             1535
12640 do 12641 j=1,ni                                                      1535
      if(ju(j).eq.0)goto 12641                                             1536
      xm(j)=dot_product(w,x(:,j))                                          1536
      x(:,j)=x(:,j)-xm(j)                                                  1537
      xv(j)=dot_product(w,x(:,j)**2)                                       1538
      if(isd .le. 0)goto 12661                                             1538
      xs(j)=sqrt(xv(j))                                                    1538
      x(:,j)=x(:,j)/xs(j)                                                  1538
      xv(j)=1.0                                                            1538
12661 continue                                                             1539
12641 continue                                                             1540
12642 continue                                                             1540
      return                                                               1541
      end                                                                  1542
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1544 
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)           1545
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1546
      integer ju(ni),m(nx),kin(nlam)                                       1547
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1552
      allocate(b(0:ni),stat=jerr)                                          1553
      allocate(xv(1:ni),stat=ierr)                                         1553
      jerr=jerr+ierr                                                       1554
      allocate(ga(1:ni),stat=ierr)                                         1554
      jerr=jerr+ierr                                                       1555
      allocate(bs(0:ni),stat=ierr)                                         1555
      jerr=jerr+ierr                                                       1556
      allocate(mm(1:ni),stat=ierr)                                         1556
      jerr=jerr+ierr                                                       1557
      allocate(ixx(1:ni),stat=ierr)                                        1557
      jerr=jerr+ierr                                                       1558
      allocate(r(1:no),stat=ierr)                                          1558
      jerr=jerr+ierr                                                       1559
      allocate(v(1:no),stat=ierr)                                          1559
      jerr=jerr+ierr                                                       1560
      allocate(q(1:no),stat=ierr)                                          1560
      jerr=jerr+ierr                                                       1561
      if(jerr.ne.0) return                                                 1562
      fmax=log(1.0/pmin-1.0)                                               1562
      fmin=-fmax                                                           1562
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1563
      bta=parm                                                             1563
      omb=1.0-bta                                                          1564
      q0=dot_product(w,y)                                                  1564
      if(q0 .gt. pmin)goto 12681                                           1564
      jerr=8001                                                            1564
      return                                                               1564
12681 continue                                                             1565
      if(q0 .lt. 1.0-pmin)goto 12701                                       1565
      jerr=9001                                                            1565
      return                                                               1565
12701 continue                                                             1566
      if(intr.eq.0.0) q0=0.5                                               1567
      ixx=0                                                                1567
      al=0.0                                                               1567
      bz=0.0                                                               1567
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    1568
      if(nonzero(no,g) .ne. 0)goto 12721                                   1568
      vi=q0*(1.0-q0)                                                       1568
      b(0)=bz                                                              1568
      v=vi*w                                                               1569
      r=w*(y-q0)                                                           1569
      q=q0                                                                 1569
      xmz=vi                                                               1569
      dev1=-(bz*q0+log(1.0-q0))                                            1570
      goto 12731                                                           1571
12721 continue                                                             1571
      b(0)=0.0                                                             1572
      if(intr .eq. 0)goto 12751                                            1572
      b(0)=azero(no,y,g,w,jerr)                                            1572
      if(jerr.ne.0) return                                                 1572
12751 continue                                                             1573
      q=1.0/(1.0+exp(-b(0)-g))                                             1573
      v=w*q*(1.0-q)                                                        1573
      r=w*(y-q)                                                            1573
      xmz=sum(v)                                                           1574
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1575
12731 continue                                                             1576
12711 continue                                                             1576
      if(kopt .le. 0)goto 12771                                            1577
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                            1577
      xv=0.25                                                              1577
      goto 12801                                                           1578
12791 continue                                                             1578
12810 do 12811 j=1,ni                                                      1578
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1578
12811 continue                                                             1578
12812 continue                                                             1578
12801 continue                                                             1579
12781 continue                                                             1579
12771 continue                                                             1580
      dev0=dev1                                                            1581
12820 do 12821 i=1,no                                                      1581
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1582
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1583
12821 continue                                                             1584
12822 continue                                                             1584
      if(flmin .ge. 1.0)goto 12841                                         1584
      eqs=max(eps,flmin)                                                   1584
      alf=eqs**(1.0/(nlam-1))                                              1584
12841 continue                                                             1585
      m=0                                                                  1585
      mm=0                                                                 1585
      nlp=0                                                                1585
      nin=nlp                                                              1585
      mnl=min(mnlam,nlam)                                                  1585
      bs=0.0                                                               1585
      b(1:ni)=0.0                                                          1586
      shr=shri*dev0                                                        1587
12850 do 12851 j=1,ni                                                      1587
      if(ju(j).eq.0)goto 12851                                             1587
      ga(j)=abs(dot_product(r,x(:,j)))                                     1587
12851 continue                                                             1588
12852 continue                                                             1588
12860 do 12861 ilm=1,nlam                                                  1588
      al0=al                                                               1589
      if(flmin .lt. 1.0)goto 12881                                         1589
      al=ulam(ilm)                                                         1589
      goto 12871                                                           1590
12881 if(ilm .le. 2)goto 12891                                             1590
      al=al*alf                                                            1590
      goto 12871                                                           1591
12891 if(ilm .ne. 1)goto 12901                                             1591
      al=big                                                               1591
      goto 12911                                                           1592
12901 continue                                                             1592
      al0=0.0                                                              1593
12920 do 12921 j=1,ni                                                      1593
      if(ju(j).eq.0)goto 12921                                             1593
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1593
12921 continue                                                             1594
12922 continue                                                             1594
      al0=al0/max(bta,1.0e-3)                                              1594
      al=alf*al0                                                           1595
12911 continue                                                             1596
12871 continue                                                             1596
      al2=al*omb                                                           1596
      al1=al*bta                                                           1596
      tlam=bta*(2.0*al-al0)                                                1597
12930 do 12931 k=1,ni                                                      1597
      if(ixx(k).eq.1)goto 12931                                            1597
      if(ju(k).eq.0)goto 12931                                             1598
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1599
12931 continue                                                             1600
12932 continue                                                             1600
10880 continue                                                             1601
12940 continue                                                             1601
12941 continue                                                             1601
      bs(0)=b(0)                                                           1601
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1602
      if(kopt .ne. 0)goto 12961                                            1603
12970 do 12971 j=1,ni                                                      1603
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1603
12971 continue                                                             1604
12972 continue                                                             1604
12961 continue                                                             1605
12980 continue                                                             1605
12981 continue                                                             1605
      nlp=nlp+1                                                            1605
      dlx=0.0                                                              1606
12990 do 12991 k=1,ni                                                      1606
      if(ixx(k).eq.0)goto 12991                                            1607
      bk=b(k)                                                              1607
      gk=dot_product(r,x(:,k))                                             1608
      u=gk+xv(k)*b(k)                                                      1608
      au=abs(u)-vp(k)*al1                                                  1609
      if(au .gt. 0.0)goto 13011                                            1609
      b(k)=0.0                                                             1609
      goto 13021                                                           1610
13011 continue                                                             1611
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1612
13021 continue                                                             1613
13001 continue                                                             1613
      d=b(k)-bk                                                            1613
      if(abs(d).le.0.0)goto 12991                                          1613
      dlx=max(dlx,xv(k)*d**2)                                              1614
      r=r-d*v*x(:,k)                                                       1615
      if(mm(k) .ne. 0)goto 13041                                           1615
      nin=nin+1                                                            1615
      if(nin.gt.nx)goto 12992                                              1616
      mm(k)=nin                                                            1616
      m(nin)=k                                                             1617
13041 continue                                                             1618
12991 continue                                                             1619
12992 continue                                                             1619
      if(nin.gt.nx)goto 12982                                              1620
      d=0.0                                                                1620
      if(intr.ne.0) d=sum(r)/xmz                                           1621
      if(d .eq. 0.0)goto 13061                                             1621
      b(0)=b(0)+d                                                          1621
      dlx=max(dlx,xmz*d**2)                                                1621
      r=r-d*v                                                              1621
13061 continue                                                             1622
      if(dlx.lt.shr)goto 12982                                             1622
      if(nlp .le. maxit)goto 13081                                         1622
      jerr=-ilm                                                            1622
      return                                                               1622
13081 continue                                                             1623
13090 continue                                                             1623
13091 continue                                                             1623
      nlp=nlp+1                                                            1623
      dlx=0.0                                                              1624
13100 do 13101 l=1,nin                                                     1624
      k=m(l)                                                               1624
      bk=b(k)                                                              1625
      gk=dot_product(r,x(:,k))                                             1626
      u=gk+xv(k)*b(k)                                                      1626
      au=abs(u)-vp(k)*al1                                                  1627
      if(au .gt. 0.0)goto 13121                                            1627
      b(k)=0.0                                                             1627
      goto 13131                                                           1628
13121 continue                                                             1629
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1630
13131 continue                                                             1631
13111 continue                                                             1631
      d=b(k)-bk                                                            1631
      if(abs(d).le.0.0)goto 13101                                          1631
      dlx=max(dlx,xv(k)*d**2)                                              1632
      r=r-d*v*x(:,k)                                                       1633
13101 continue                                                             1634
13102 continue                                                             1634
      d=0.0                                                                1634
      if(intr.ne.0) d=sum(r)/xmz                                           1635
      if(d .eq. 0.0)goto 13151                                             1635
      b(0)=b(0)+d                                                          1635
      dlx=max(dlx,xmz*d**2)                                                1635
      r=r-d*v                                                              1635
13151 continue                                                             1636
      if(dlx.lt.shr)goto 13092                                             1636
      if(nlp .le. maxit)goto 13171                                         1636
      jerr=-ilm                                                            1636
      return                                                               1636
13171 continue                                                             1637
      goto 13091                                                           1638
13092 continue                                                             1638
      goto 12981                                                           1639
12982 continue                                                             1639
      if(nin.gt.nx)goto 12942                                              1640
13180 do 13181 i=1,no                                                      1640
      fi=b(0)+g(i)                                                         1641
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1642
      if(fi .ge. fmin)goto 13201                                           1642
      q(i)=0.0                                                             1642
      goto 13191                                                           1642
13201 if(fi .le. fmax)goto 13211                                           1642
      q(i)=1.0                                                             1642
      goto 13221                                                           1643
13211 continue                                                             1643
      q(i)=1.0/(1.0+exp(-fi))                                              1643
13221 continue                                                             1644
13191 continue                                                             1644
13181 continue                                                             1645
13182 continue                                                             1645
      v=w*q*(1.0-q)                                                        1645
      xmz=sum(v)                                                           1645
      if(xmz.le.vmin)goto 12942                                            1645
      r=w*(y-q)                                                            1646
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                           1646
      ix=0                                                                 1647
13250 do 13251 j=1,nin                                                     1647
      k=m(j)                                                               1648
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                           1648
      ix=1                                                                 1648
      goto 13252                                                           1649
13251 continue                                                             1650
13252 continue                                                             1650
      if(ix .ne. 0)goto 13271                                              1651
13280 do 13281 k=1,ni                                                      1651
      if(ixx(k).eq.1)goto 13281                                            1651
      if(ju(k).eq.0)goto 13281                                             1652
      ga(k)=abs(dot_product(r,x(:,k)))                                     1653
      if(ga(k) .le. al1*vp(k))goto 13301                                   1653
      ixx(k)=1                                                             1653
      ix=1                                                                 1653
13301 continue                                                             1654
13281 continue                                                             1655
13282 continue                                                             1655
      if(ix.eq.1) go to 10880                                              1656
      goto 12942                                                           1657
13271 continue                                                             1658
13241 continue                                                             1659
      goto 12941                                                           1660
12942 continue                                                             1660
      if(nin .le. nx)goto 13321                                            1660
      jerr=-10000-ilm                                                      1660
      goto 12862                                                           1660
13321 continue                                                             1661
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1661
      kin(ilm)=nin                                                         1662
      a0(ilm)=b(0)                                                         1662
      alm(ilm)=al                                                          1662
      lmu=ilm                                                              1663
      devi=dev2(no,w,y,q,pmin)                                             1664
      dev(ilm)=(dev1-devi)/dev0                                            1664
      if(xmz.le.vmin)goto 12862                                            1665
      if(ilm.lt.mnl)goto 12861                                             1665
      if(flmin.ge.1.0)goto 12861                                           1666
      me=0                                                                 1666
13330 do 13331 j=1,nin                                                     1666
      if(a(j,ilm).ne.0.0) me=me+1                                          1666
13331 continue                                                             1666
13332 continue                                                             1666
      if(me.gt.ne)goto 12862                                               1667
      if(dev(ilm).gt.devmax)goto 12862                                     1667
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                             1668
12861 continue                                                             1669
12862 continue                                                             1669
      g=log(q/(1.0-q))                                                     1670
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1671
      return                                                               1672
      end                                                                  1673
      function dev2(n,w,y,p,pmin)                                          1674
      real w(n),y(n),p(n)                                                  1675
      pmax=1.0-pmin                                                        1675
      s=0.0                                                                1676
13340 do 13341 i=1,n                                                       1676
      pi=min(max(pmin,p(i)),pmax)                                          1677
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1678
13341 continue                                                             1679
13342 continue                                                             1679
      dev2=s                                                               1680
      return                                                               1681
      end                                                                  1682
      function azero(n,y,g,q,jerr)                                         1683
      parameter(eps=1.0e-7)                                                1684
      real y(n),g(n),q(n)                                                  1685
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1689
      allocate(p(1:n),stat=ierr)                                           1689
      jerr=jerr+ierr                                                       1690
      allocate(w(1:n),stat=ierr)                                           1690
      jerr=jerr+ierr                                                       1691
      if(jerr.ne.0) return                                                 1692
      az=0.0                                                               1692
      e=exp(-g)                                                            1692
      qy=dot_product(q,y)                                                  1692
      p=1.0/(1.0+e)                                                        1693
13350 continue                                                             1693
13351 continue                                                             1693
      w=q*p*(1.0-p)                                                        1694
      d=(qy-dot_product(q,p))/sum(w)                                       1694
      az=az+d                                                              1694
      if(abs(d).lt.eps)goto 13352                                          1695
      ea0=exp(-az)                                                         1695
      p=1.0/(1.0+ea0*e)                                                    1696
      goto 13351                                                           1697
13352 continue                                                             1697
      azero=az                                                             1698
      deallocate(e,p,w)                                                    1699
      return                                                               1700
      end                                                                  1701
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1703 
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1704
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          1705
      integer ju(ni),m(nx),kin(nlam)                                       1706
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1717
      exmn=-exmx                                                           1718
      allocate(r(1:no),stat=ierr)                                          1718
      jerr=jerr+ierr                                                       1719
      allocate(v(1:no),stat=ierr)                                          1719
      jerr=jerr+ierr                                                       1720
      allocate(mm(1:ni),stat=ierr)                                         1720
      jerr=jerr+ierr                                                       1721
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1721
      jerr=jerr+ierr                                                       1722
      allocate(sxp(1:no),stat=ierr)                                        1722
      jerr=jerr+ierr                                                       1723
      allocate(sxpl(1:no),stat=ierr)                                       1723
      jerr=jerr+ierr                                                       1724
      allocate(di(1:no),stat=ierr)                                         1724
      jerr=jerr+ierr                                                       1725
      allocate(ga(1:ni),stat=ierr)                                         1725
      jerr=jerr+ierr                                                       1726
      allocate(ixx(1:ni),stat=ierr)                                        1726
      jerr=jerr+ierr                                                       1727
      if(jerr.ne.0) return                                                 1728
      pmax=1.0-pmin                                                        1728
      emin=pmin/pmax                                                       1728
      emax=1.0/emin                                                        1729
      pfm=(1.0+pmin)*pmin                                                  1729
      pfx=(1.0-pmin)*pmax                                                  1729
      vmin=pfm*pmax                                                        1730
      bta=parm                                                             1730
      omb=1.0-bta                                                          1730
      dev1=0.0                                                             1730
      dev0=0.0                                                             1731
13360 do 13361 ic=1,nc                                                     1731
      q0=dot_product(w,y(:,ic))                                            1732
      if(q0 .gt. pmin)goto 13381                                           1732
      jerr =8000+ic                                                        1732
      return                                                               1732
13381 continue                                                             1733
      if(q0 .lt. 1.0-pmin)goto 13401                                       1733
      jerr =9000+ic                                                        1733
      return                                                               1733
13401 continue                                                             1734
      if(intr .ne. 0)goto 13421                                            1734
      q0=1.0/nc                                                            1734
      b(0,ic)=0.0                                                          1734
      goto 13431                                                           1735
13421 continue                                                             1735
      b(0,ic)=log(q0)                                                      1735
      dev1=dev1-q0*b(0,ic)                                                 1735
13431 continue                                                             1736
13411 continue                                                             1736
      b(1:ni,ic)=0.0                                                       1737
13361 continue                                                             1738
13362 continue                                                             1738
      if(intr.eq.0) dev1=log(float(nc))                                    1738
      ixx=0                                                                1738
      al=0.0                                                               1739
      if(nonzero(no*nc,g) .ne. 0)goto 13451                                1740
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1740
      sxp=0.0                                                              1741
13460 do 13461 ic=1,nc                                                     1741
      q(:,ic)=exp(b(0,ic))                                                 1741
      sxp=sxp+q(:,ic)                                                      1741
13461 continue                                                             1742
13462 continue                                                             1742
      goto 13471                                                           1743
13451 continue                                                             1743
13480 do 13481 i=1,no                                                      1743
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1743
13481 continue                                                             1743
13482 continue                                                             1743
      sxp=0.0                                                              1744
      if(intr .ne. 0)goto 13501                                            1744
      b(0,:)=0.0                                                           1744
      goto 13511                                                           1745
13501 continue                                                             1745
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1745
      if(jerr.ne.0) return                                                 1745
13511 continue                                                             1746
13491 continue                                                             1746
      dev1=0.0                                                             1747
13520 do 13521 ic=1,nc                                                     1747
      q(:,ic)=b(0,ic)+g(:,ic)                                              1748
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1749
      q(:,ic)=exp(q(:,ic))                                                 1749
      sxp=sxp+q(:,ic)                                                      1750
13521 continue                                                             1751
13522 continue                                                             1751
      sxpl=w*log(sxp)                                                      1751
13530 do 13531 ic=1,nc                                                     1751
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1751
13531 continue                                                             1752
13532 continue                                                             1752
13471 continue                                                             1753
13441 continue                                                             1753
13540 do 13541 ic=1,nc                                                     1753
13550 do 13551 i=1,no                                                      1753
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1753
13551 continue                                                             1753
13552 continue                                                             1753
13541 continue                                                             1754
13542 continue                                                             1754
      dev0=dev0+dev1                                                       1755
      if(kopt .le. 0)goto 13571                                            1756
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                            1756
      xv=0.25                                                              1756
      goto 13601                                                           1757
13591 continue                                                             1757
13610 do 13611 j=1,ni                                                      1757
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1757
13611 continue                                                             1757
13612 continue                                                             1757
13601 continue                                                             1758
13581 continue                                                             1758
13571 continue                                                             1759
      if(flmin .ge. 1.0)goto 13631                                         1759
      eqs=max(eps,flmin)                                                   1759
      alf=eqs**(1.0/(nlam-1))                                              1759
13631 continue                                                             1760
      m=0                                                                  1760
      mm=0                                                                 1760
      nin=0                                                                1760
      nlp=0                                                                1760
      mnl=min(mnlam,nlam)                                                  1760
      bs=0.0                                                               1760
      shr=shri*dev0                                                        1761
      ga=0.0                                                               1762
13640 do 13641 ic=1,nc                                                     1762
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1763
13650 do 13651 j=1,ni                                                      1763
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1763
13651 continue                                                             1764
13652 continue                                                             1764
13641 continue                                                             1765
13642 continue                                                             1765
13660 do 13661 ilm=1,nlam                                                  1765
      al0=al                                                               1766
      if(flmin .lt. 1.0)goto 13681                                         1766
      al=ulam(ilm)                                                         1766
      goto 13671                                                           1767
13681 if(ilm .le. 2)goto 13691                                             1767
      al=al*alf                                                            1767
      goto 13671                                                           1768
13691 if(ilm .ne. 1)goto 13701                                             1768
      al=big                                                               1768
      goto 13711                                                           1769
13701 continue                                                             1769
      al0=0.0                                                              1770
13720 do 13721 j=1,ni                                                      1770
      if(ju(j).eq.0)goto 13721                                             1770
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1770
13721 continue                                                             1771
13722 continue                                                             1771
      al0=al0/max(bta,1.0e-3)                                              1771
      al=alf*al0                                                           1772
13711 continue                                                             1773
13671 continue                                                             1773
      al2=al*omb                                                           1773
      al1=al*bta                                                           1773
      tlam=bta*(2.0*al-al0)                                                1774
13730 do 13731 k=1,ni                                                      1774
      if(ixx(k).eq.1)goto 13731                                            1774
      if(ju(k).eq.0)goto 13731                                             1775
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1776
13731 continue                                                             1777
13732 continue                                                             1777
10880 continue                                                             1778
13740 continue                                                             1778
13741 continue                                                             1778
      ix=0                                                                 1778
      jx=ix                                                                1778
      ig=0                                                                 1779
13750 do 13751 ic=1,nc                                                     1779
      bs(0,ic)=b(0,ic)                                                     1780
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1781
      xmz=0.0                                                              1782
13760 do 13761 i=1,no                                                      1782
      pic=q(i,ic)/sxp(i)                                                   1783
      if(pic .ge. pfm)goto 13781                                           1783
      pic=0.0                                                              1783
      v(i)=0.0                                                             1783
      goto 13771                                                           1784
13781 if(pic .le. pfx)goto 13791                                           1784
      pic=1.0                                                              1784
      v(i)=0.0                                                             1784
      goto 13801                                                           1785
13791 continue                                                             1785
      v(i)=w(i)*pic*(1.0-pic)                                              1785
      xmz=xmz+v(i)                                                         1785
13801 continue                                                             1786
13771 continue                                                             1786
      r(i)=w(i)*(y(i,ic)-pic)                                              1787
13761 continue                                                             1788
13762 continue                                                             1788
      if(xmz.le.vmin)goto 13751                                            1788
      ig=1                                                                 1789
      if(kopt .ne. 0)goto 13821                                            1790
13830 do 13831 j=1,ni                                                      1790
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1790
13831 continue                                                             1791
13832 continue                                                             1791
13821 continue                                                             1792
13840 continue                                                             1792
13841 continue                                                             1792
      nlp=nlp+1                                                            1792
      dlx=0.0                                                              1793
13850 do 13851 k=1,ni                                                      1793
      if(ixx(k).eq.0)goto 13851                                            1794
      bk=b(k,ic)                                                           1794
      gk=dot_product(r,x(:,k))                                             1795
      u=gk+xv(k,ic)*b(k,ic)                                                1795
      au=abs(u)-vp(k)*al1                                                  1796
      if(au .gt. 0.0)goto 13871                                            1796
      b(k,ic)=0.0                                                          1796
      goto 13881                                                           1797
13871 continue                                                             1798
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1800 
     *)
13881 continue                                                             1801
13861 continue                                                             1801
      d=b(k,ic)-bk                                                         1801
      if(abs(d).le.0.0)goto 13851                                          1802
      dlx=max(dlx,xv(k,ic)*d**2)                                           1802
      r=r-d*v*x(:,k)                                                       1803
      if(mm(k) .ne. 0)goto 13901                                           1803
      nin=nin+1                                                            1804
      if(nin .le. nx)goto 13921                                            1804
      jx=1                                                                 1804
      goto 13852                                                           1804
13921 continue                                                             1805
      mm(k)=nin                                                            1805
      m(nin)=k                                                             1806
13901 continue                                                             1807
13851 continue                                                             1808
13852 continue                                                             1808
      if(jx.gt.0)goto 13842                                                1809
      d=0.0                                                                1809
      if(intr.ne.0) d=sum(r)/xmz                                           1810
      if(d .eq. 0.0)goto 13941                                             1810
      b(0,ic)=b(0,ic)+d                                                    1810
      dlx=max(dlx,xmz*d**2)                                                1810
      r=r-d*v                                                              1810
13941 continue                                                             1811
      if(dlx.lt.shr)goto 13842                                             1812
      if(nlp .le. maxit)goto 13961                                         1812
      jerr=-ilm                                                            1812
      return                                                               1812
13961 continue                                                             1813
13970 continue                                                             1813
13971 continue                                                             1813
      nlp=nlp+1                                                            1813
      dlx=0.0                                                              1814
13980 do 13981 l=1,nin                                                     1814
      k=m(l)                                                               1814
      bk=b(k,ic)                                                           1815
      gk=dot_product(r,x(:,k))                                             1816
      u=gk+xv(k,ic)*b(k,ic)                                                1816
      au=abs(u)-vp(k)*al1                                                  1817
      if(au .gt. 0.0)goto 14001                                            1817
      b(k,ic)=0.0                                                          1817
      goto 14011                                                           1818
14001 continue                                                             1819
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1821 
     *)
14011 continue                                                             1822
13991 continue                                                             1822
      d=b(k,ic)-bk                                                         1822
      if(abs(d).le.0.0)goto 13981                                          1823
      dlx=max(dlx,xv(k,ic)*d**2)                                           1823
      r=r-d*v*x(:,k)                                                       1824
13981 continue                                                             1825
13982 continue                                                             1825
      d=0.0                                                                1825
      if(intr.ne.0) d=sum(r)/xmz                                           1826
      if(d .eq. 0.0)goto 14031                                             1826
      b(0,ic)=b(0,ic)+d                                                    1827
      dlx=max(dlx,xmz*d**2)                                                1827
      r=r-d*v                                                              1828
14031 continue                                                             1829
      if(dlx.lt.shr)goto 13972                                             1829
      if(nlp .le. maxit)goto 14051                                         1829
      jerr=-ilm                                                            1829
      return                                                               1829
14051 continue                                                             1830
      goto 13971                                                           1831
13972 continue                                                             1831
      goto 13841                                                           1832
13842 continue                                                             1832
      if(jx.gt.0)goto 13752                                                1833
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1834
      if(ix .ne. 0)goto 14071                                              1835
14080 do 14081 j=1,nin                                                     1835
      k=m(j)                                                               1836
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101                1836
      ix=1                                                                 1836
      goto 14082                                                           1836
14101 continue                                                             1837
14081 continue                                                             1838
14082 continue                                                             1838
14071 continue                                                             1839
14110 do 14111 i=1,no                                                      1839
      fi=b(0,ic)+g(i,ic)                                                   1841
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1842
      fi=min(max(exmn,fi),exmx)                                            1842
      sxp(i)=sxp(i)-q(i,ic)                                                1843
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1844
      sxp(i)=sxp(i)+q(i,ic)                                                1845
14111 continue                                                             1846
14112 continue                                                             1846
13751 continue                                                             1847
13752 continue                                                             1847
      s=-sum(b(0,:))/nc                                                    1847
      b(0,:)=b(0,:)+s                                                      1847
      di=s                                                                 1848
14120 do 14121 j=1,nin                                                     1848
      l=m(j)                                                               1849
      if(vp(l) .gt. 0.0)goto 14141                                         1849
      s=sum(b(l,:))/nc                                                     1849
      goto 14151                                                           1850
14141 continue                                                             1850
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1850
14151 continue                                                             1851
14131 continue                                                             1851
      b(l,:)=b(l,:)-s                                                      1851
      di=di-s*x(:,l)                                                       1852
14121 continue                                                             1853
14122 continue                                                             1853
      di=exp(di)                                                           1853
      sxp=sxp*di                                                           1853
14160 do 14161 ic=1,nc                                                     1853
      q(:,ic)=q(:,ic)*di                                                   1853
14161 continue                                                             1854
14162 continue                                                             1854
      if(jx.gt.0)goto 13742                                                1854
      if(ig.eq.0)goto 13742                                                1855
      if(ix .ne. 0)goto 14181                                              1856
14190 do 14191 k=1,ni                                                      1856
      if(ixx(k).eq.1)goto 14191                                            1856
      if(ju(k).eq.0)goto 14191                                             1856
      ga(k)=0.0                                                            1856
14191 continue                                                             1857
14192 continue                                                             1857
14200 do 14201 ic=1,nc                                                     1857
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1858
14210 do 14211 k=1,ni                                                      1858
      if(ixx(k).eq.1)goto 14211                                            1858
      if(ju(k).eq.0)goto 14211                                             1859
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1860
14211 continue                                                             1861
14212 continue                                                             1861
14201 continue                                                             1862
14202 continue                                                             1862
14220 do 14221 k=1,ni                                                      1862
      if(ixx(k).eq.1)goto 14221                                            1862
      if(ju(k).eq.0)goto 14221                                             1863
      if(ga(k) .le. al1*vp(k))goto 14241                                   1863
      ixx(k)=1                                                             1863
      ix=1                                                                 1863
14241 continue                                                             1864
14221 continue                                                             1865
14222 continue                                                             1865
      if(ix.eq.1) go to 10880                                              1866
      goto 13742                                                           1867
14181 continue                                                             1868
      goto 13741                                                           1869
13742 continue                                                             1869
      if(jx .le. 0)goto 14261                                              1869
      jerr=-10000-ilm                                                      1869
      goto 13662                                                           1869
14261 continue                                                             1869
      devi=0.0                                                             1870
14270 do 14271 ic=1,nc                                                     1871
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1871
      a0(ic,ilm)=b(0,ic)                                                   1872
14280 do 14281 i=1,no                                                      1872
      if(y(i,ic).le.0.0)goto 14281                                         1873
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1874
14281 continue                                                             1875
14282 continue                                                             1875
14271 continue                                                             1876
14272 continue                                                             1876
      kin(ilm)=nin                                                         1876
      alm(ilm)=al                                                          1876
      lmu=ilm                                                              1877
      dev(ilm)=(dev1-devi)/dev0                                            1877
      if(ig.eq.0)goto 13662                                                1878
      if(ilm.lt.mnl)goto 13661                                             1878
      if(flmin.ge.1.0)goto 13661                                           1879
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662             1880
      if(dev(ilm).gt.devmax)goto 13662                                     1880
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                             1881
13661 continue                                                             1882
13662 continue                                                             1882
      g=log(q)                                                             1882
14290 do 14291 i=1,no                                                      1882
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1882
14291 continue                                                             1883
14292 continue                                                             1883
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1884
      return                                                               1885
      end                                                                  1886
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1887
      parameter(eps=1.0e-7)                                                1888
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1889
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1894
      jerr=jerr+ierr                                                       1895
      if(jerr.ne.0) return                                                 1896
      az=0.0                                                               1896
      e=exp(g)                                                             1896
14300 do 14301 i=1,n                                                       1896
      s(i)=sum(e(i,:))                                                     1896
14301 continue                                                             1897
14302 continue                                                             1897
14310 continue                                                             1897
14311 continue                                                             1897
      dm=0.0                                                               1898
14320 do 14321 k=1,kk                                                      1898
      t=0.0                                                                1898
      u=t                                                                  1899
14330 do 14331 i=1,n                                                       1899
      pik=e(i,k)/s(i)                                                      1900
      t=t+q(i)*(y(i,k)-pik)                                                1900
      u=u+q(i)*pik*(1.0-pik)                                               1901
14331 continue                                                             1902
14332 continue                                                             1902
      d=t/u                                                                1902
      az(k)=az(k)+d                                                        1902
      ed=exp(d)                                                            1902
      dm=max(dm,abs(d))                                                    1903
14340 do 14341 i=1,n                                                       1903
      z=e(i,k)                                                             1903
      e(i,k)=z*ed                                                          1903
      s(i)=s(i)-z+e(i,k)                                                   1903
14341 continue                                                             1904
14342 continue                                                             1904
14321 continue                                                             1905
14322 continue                                                             1905
      if(dm.lt.eps)goto 14312                                              1905
      goto 14311                                                           1906
14312 continue                                                             1906
      az=az-sum(az)/kk                                                     1907
      deallocate(e,s)                                                      1908
      return                                                               1909
      end                                                                  1910
      function elc(parm,n,cl,a,m)                                          1911
      real a(n),cl(2)                                                      1911
      integer m(n)                                                         1912
      fn=n                                                                 1912
      am=sum(a)/fn                                                         1913
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14361                       1913
      elc=am                                                               1913
      go to 14370                                                          1913
14361 continue                                                             1914
14380 do 14381 i=1,n                                                       1914
      m(i)=i                                                               1914
14381 continue                                                             1914
14382 continue                                                             1914
      call psort7(a,m,1,n)                                                 1915
      if(a(m(1)) .ne. a(m(n)))goto 14401                                   1915
      elc=a(1)                                                             1915
      go to 14370                                                          1915
14401 continue                                                             1916
      if(mod(n,2) .ne. 1)goto 14421                                        1916
      ad=a(m(n/2+1))                                                       1916
      goto 14431                                                           1917
14421 continue                                                             1917
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1917
14431 continue                                                             1918
14411 continue                                                             1918
      if(parm .ne. 1.0)goto 14451                                          1918
      elc=ad                                                               1918
      go to 14370                                                          1918
14451 continue                                                             1919
      b1=min(am,ad)                                                        1919
      b2=max(am,ad)                                                        1919
      k2=1                                                                 1920
14460 continue                                                             1920
14461 if(a(m(k2)).gt.b1)goto 14462                                         1920
      k2=k2+1                                                              1920
      goto 14461                                                           1920
14462 continue                                                             1920
      k1=k2-1                                                              1921
14470 continue                                                             1921
14471 if(a(m(k2)).ge.b2)goto 14472                                         1921
      k2=k2+1                                                              1921
      goto 14471                                                           1922
14472 continue                                                             1922
      r=parm/((1.0-parm)*fn)                                               1922
      is=0                                                                 1922
      sm=n-2*(k1-1)                                                        1923
14480 do 14481 k=k1,k2-1                                                   1923
      sm=sm-2.0                                                            1923
      s=r*sm+am                                                            1924
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14501                   1924
      is=k                                                                 1924
      goto 14482                                                           1924
14501 continue                                                             1925
14481 continue                                                             1926
14482 continue                                                             1926
      if(is .eq. 0)goto 14521                                              1926
      elc=s                                                                1926
      go to 14370                                                          1926
14521 continue                                                             1926
      r2=2.0*r                                                             1926
      s1=a(m(k1))                                                          1926
      am2=2.0*am                                                           1927
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1927
      elc=s1                                                               1928
14530 do 14531 k=k1+1,k2                                                   1928
      s=a(m(k))                                                            1928
      if(s.eq.s1)goto 14531                                                1929
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1930
      if(c .ge. cri)goto 14551                                             1930
      cri=c                                                                1930
      elc=s                                                                1930
14551 continue                                                             1930
      s1=s                                                                 1931
14531 continue                                                             1932
14532 continue                                                             1932
14370 continue                                                             1932
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    1933
      return                                                               1934
      end                                                                  1935
      function nintot(ni,nx,nc,a,m,nin,is)                                 1936
      real a(nx,nc)                                                        1936
      integer m(nx),is(ni)                                                 1937
      is=0                                                                 1937
      nintot=0                                                             1938
14560 do 14561 ic=1,nc                                                     1938
14570 do 14571 j=1,nin                                                     1938
      k=m(j)                                                               1938
      if(is(k).ne.0)goto 14571                                             1939
      if(a(j,ic).eq.0.0)goto 14571                                         1939
      is(k)=k                                                              1939
      nintot=nintot+1                                                      1940
14571 continue                                                             1940
14572 continue                                                             1940
14561 continue                                                             1941
14562 continue                                                             1941
      return                                                               1942
      end                                                                  1943
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1944
      real ca(nx,nc),a(ni,nc)                                              1944
      integer ia(nx)                                                       1945
      a=0.0                                                                1946
14580 do 14581 ic=1,nc                                                     1946
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1946
14581 continue                                                             1947
14582 continue                                                             1947
      return                                                               1948
      end                                                                  1949
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1950
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1950
      integer ia(nx)                                                       1951
14590 do 14591 i=1,nt                                                      1951
14600 do 14601 ic=1,nc                                                     1951
      ans(ic,i)=a0(ic)                                                     1953
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1954 
     *:nin)))
14601 continue                                                             1954
14602 continue                                                             1954
14591 continue                                                             1955
14592 continue                                                             1955
      return                                                               1956
      end                                                                  1957
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   1959 
     *,flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1960
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)         1961
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1962
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14621                                    1966
      jerr=10000                                                           1966
      return                                                               1966
14621 continue                                                             1967
      allocate(ww(1:no),stat=jerr)                                         1968
      allocate(ju(1:ni),stat=ierr)                                         1968
      jerr=jerr+ierr                                                       1969
      allocate(vq(1:ni),stat=ierr)                                         1969
      jerr=jerr+ierr                                                       1970
      allocate(xm(1:ni),stat=ierr)                                         1970
      jerr=jerr+ierr                                                       1971
      allocate(xs(1:ni),stat=ierr)                                         1971
      jerr=jerr+ierr                                                       1972
      if(kopt .ne. 2)goto 14641                                            1972
      allocate(xv(1:ni),stat=ierr)                                         1972
      jerr=jerr+ierr                                                       1972
14641 continue                                                             1973
      if(jerr.ne.0) return                                                 1974
      call spchkvars(no,ni,x,ix,ju)                                        1975
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1976
      if(maxval(ju) .gt. 0)goto 14661                                      1976
      jerr=7777                                                            1976
      return                                                               1976
14661 continue                                                             1977
      vq=max(0.0,vp)                                                       1977
      vq=vq*ni/sum(vq)                                                     1978
14670 do 14671 i=1,no                                                      1978
      ww(i)=sum(y(i,:))                                                    1978
      y(i,:)=y(i,:)/ww(i)                                                  1978
14671 continue                                                             1978
14672 continue                                                             1978
      sw=sum(ww)                                                           1978
      ww=ww/sw                                                             1979
      if(nc .ne. 1)goto 14691                                              1979
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                1980
      if(isd .le. 0)goto 14711                                             1980
14720 do 14721 j=1,ni                                                      1980
      cl(:,j)=cl(:,j)*xs(j)                                                1980
14721 continue                                                             1980
14722 continue                                                             1980
14711 continue                                                             1981
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   1984 
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14681                                                           1985
14691 if(kopt .ne. 2)goto 14731                                            1986
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)         1987
      if(isd .le. 0)goto 14751                                             1987
14760 do 14761 j=1,ni                                                      1987
      cl(:,j)=cl(:,j)*xs(j)                                                1987
14761 continue                                                             1987
14762 continue                                                             1987
14751 continue                                                             1988
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   1990 
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14771                                                           1991
14731 continue                                                             1991
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                1992
      if(isd .le. 0)goto 14791                                             1992
14800 do 14801 j=1,ni                                                      1992
      cl(:,j)=cl(:,j)*xs(j)                                                1992
14801 continue                                                             1992
14802 continue                                                             1992
14791 continue                                                             1993
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   1996 
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14771 continue                                                             1997
14681 continue                                                             1997
      if(jerr.gt.0) return                                                 1997
      dev0=2.0*sw*dev0                                                     1998
14810 do 14811 k=1,lmu                                                     1998
      nk=nin(k)                                                            1999
14820 do 14821 ic=1,nc                                                     1999
      if(isd .le. 0)goto 14841                                             1999
14850 do 14851 l=1,nk                                                      1999
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1999
14851 continue                                                             1999
14852 continue                                                             1999
14841 continue                                                             2000
      if(intr .ne. 0)goto 14871                                            2000
      a0(ic,k)=0.0                                                         2000
      goto 14881                                                           2001
14871 continue                                                             2001
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            2001
14881 continue                                                             2002
14861 continue                                                             2002
14821 continue                                                             2003
14822 continue                                                             2003
14811 continue                                                             2004
14812 continue                                                             2004
      deallocate(ww,ju,vq,xm,xs)                                           2004
      if(kopt.eq.2) deallocate(xv)                                         2005
      return                                                               2006
      end                                                                  2007
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv)    2008
      real x(*),w(no),xm(ni),xs(ni),xv(ni)                                 2008
      integer ix(*),jx(*),ju(ni)                                           2009
      if(intr .ne. 0)goto 14901                                            2010
14910 do 14911 j=1,ni                                                      2010
      if(ju(j).eq.0)goto 14911                                             2010
      xm(j)=0.0                                                            2010
      jb=ix(j)                                                             2010
      je=ix(j+1)-1                                                         2011
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          2012
      if(isd .eq. 0)goto 14931                                             2012
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            2012
      vc=xv(j)-xbq                                                         2013
      xs(j)=sqrt(vc)                                                       2013
      xv(j)=1.0+vbq/vc                                                     2014
      goto 14941                                                           2015
14931 continue                                                             2015
      xs(j)=1.0                                                            2015
14941 continue                                                             2016
14921 continue                                                             2016
14911 continue                                                             2017
14912 continue                                                             2017
      return                                                               2018
14901 continue                                                             2019
14950 do 14951 j=1,ni                                                      2019
      if(ju(j).eq.0)goto 14951                                             2019
      jb=ix(j)                                                             2019
      je=ix(j+1)-1                                                         2020
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2021
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 2022
      if(isd .le. 0)goto 14971                                             2022
      xs(j)=sqrt(xv(j))                                                    2022
      xv(j)=1.0                                                            2022
14971 continue                                                             2023
14951 continue                                                             2024
14952 continue                                                             2024
      if(isd.eq.0) xs=1.0                                                  2025
      return                                                               2026
      end                                                                  2027
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)           2028
      real x(*),w(no),xm(ni),xs(ni)                                        2028
      integer ix(*),jx(*),ju(ni)                                           2029
      if(intr .ne. 0)goto 14991                                            2030
15000 do 15001 j=1,ni                                                      2030
      if(ju(j).eq.0)goto 15001                                             2030
      xm(j)=0.0                                                            2030
      jb=ix(j)                                                             2030
      je=ix(j+1)-1                                                         2031
      if(isd .eq. 0)goto 15021                                             2032
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)   2034 
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                       2035
      goto 15031                                                           2036
15021 continue                                                             2036
      xs(j)=1.0                                                            2036
15031 continue                                                             2037
15011 continue                                                             2037
15001 continue                                                             2038
15002 continue                                                             2038
      return                                                               2039
14991 continue                                                             2040
15040 do 15041 j=1,ni                                                      2040
      if(ju(j).eq.0)goto 15041                                             2040
      jb=ix(j)                                                             2040
      je=ix(j+1)-1                                                         2041
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2042
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   2043 
     *)**2)
15041 continue                                                             2044
15042 continue                                                             2044
      if(isd.eq.0) xs=1.0                                                  2045
      return                                                               2046
      end                                                                  2047
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   2050 
     *am,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,de
     *v0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)               2051
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         2052
      real xb(ni),xs(ni)                                                   2052
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2053
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2058
      allocate(b(0:ni),stat=jerr)                                          2059
      allocate(xm(0:ni),stat=ierr)                                         2059
      jerr=jerr+ierr                                                       2060
      allocate(xv(1:ni),stat=ierr)                                         2060
      jerr=jerr+ierr                                                       2061
      allocate(bs(0:ni),stat=ierr)                                         2061
      jerr=jerr+ierr                                                       2062
      allocate(ga(1:ni),stat=ierr)                                         2062
      jerr=jerr+ierr                                                       2063
      allocate(mm(1:ni),stat=ierr)                                         2063
      jerr=jerr+ierr                                                       2064
      allocate(ixx(1:ni),stat=ierr)                                        2064
      jerr=jerr+ierr                                                       2065
      allocate(q(1:no),stat=ierr)                                          2065
      jerr=jerr+ierr                                                       2066
      allocate(r(1:no),stat=ierr)                                          2066
      jerr=jerr+ierr                                                       2067
      allocate(v(1:no),stat=ierr)                                          2067
      jerr=jerr+ierr                                                       2068
      allocate(sc(1:no),stat=ierr)                                         2068
      jerr=jerr+ierr                                                       2069
      if(jerr.ne.0) return                                                 2070
      fmax=log(1.0/pmin-1.0)                                               2070
      fmin=-fmax                                                           2070
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      2071
      bta=parm                                                             2071
      omb=1.0-bta                                                          2072
      q0=dot_product(w,y)                                                  2072
      if(q0 .gt. pmin)goto 15061                                           2072
      jerr=8001                                                            2072
      return                                                               2072
15061 continue                                                             2073
      if(q0 .lt. 1.0-pmin)goto 15081                                       2073
      jerr=9001                                                            2073
      return                                                               2073
15081 continue                                                             2074
      if(intr.eq.0) q0=0.5                                                 2074
      bz=0.0                                                               2074
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    2075
      if(nonzero(no,g) .ne. 0)goto 15101                                   2075
      vi=q0*(1.0-q0)                                                       2075
      b(0)=bz                                                              2075
      v=vi*w                                                               2076
      r=w*(y-q0)                                                           2076
      q=q0                                                                 2076
      xm(0)=vi                                                             2076
      dev1=-(bz*q0+log(1.0-q0))                                            2077
      goto 15111                                                           2078
15101 continue                                                             2078
      b(0)=0.0                                                             2079
      if(intr .eq. 0)goto 15131                                            2079
      b(0)=azero(no,y,g,w,jerr)                                            2079
      if(jerr.ne.0) return                                                 2079
15131 continue                                                             2080
      q=1.0/(1.0+exp(-b(0)-g))                                             2080
      v=w*q*(1.0-q)                                                        2080
      r=w*(y-q)                                                            2080
      xm(0)=sum(v)                                                         2081
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        2082
15111 continue                                                             2083
15091 continue                                                             2083
      if(kopt .le. 0)goto 15151                                            2084
      if(isd .le. 0 .or. intr .eq. 0)goto 15171                            2084
      xv=0.25                                                              2084
      goto 15181                                                           2085
15171 continue                                                             2086
15190 do 15191 j=1,ni                                                      2086
      if(ju(j).eq.0)goto 15191                                             2086
      jb=ix(j)                                                             2086
      je=ix(j+1)-1                                                         2087
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          2088
15191 continue                                                             2089
15192 continue                                                             2089
15181 continue                                                             2090
15161 continue                                                             2090
15151 continue                                                             2091
      b(1:ni)=0.0                                                          2091
      dev0=dev1                                                            2092
15200 do 15201 i=1,no                                                      2092
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        2093
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              2094
15201 continue                                                             2095
15202 continue                                                             2095
      if(flmin .ge. 1.0)goto 15221                                         2095
      eqs=max(eps,flmin)                                                   2095
      alf=eqs**(1.0/(nlam-1))                                              2095
15221 continue                                                             2096
      m=0                                                                  2096
      mm=0                                                                 2096
      nin=0                                                                2096
      o=0.0                                                                2096
      svr=o                                                                2096
      mnl=min(mnlam,nlam)                                                  2096
      bs=0.0                                                               2096
      nlp=0                                                                2096
      nin=nlp                                                              2097
      shr=shri*dev0                                                        2097
      al=0.0                                                               2097
      ixx=0                                                                2098
15230 do 15231 j=1,ni                                                      2098
      if(ju(j).eq.0)goto 15231                                             2099
      jb=ix(j)                                                             2099
      je=ix(j+1)-1                                                         2099
      jn=ix(j+1)-ix(j)                                                     2100
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2101
      gj=dot_product(sc(1:jn),x(jb:je))                                    2102
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2103
15231 continue                                                             2104
15232 continue                                                             2104
15240 do 15241 ilm=1,nlam                                                  2104
      al0=al                                                               2105
      if(flmin .lt. 1.0)goto 15261                                         2105
      al=ulam(ilm)                                                         2105
      goto 15251                                                           2106
15261 if(ilm .le. 2)goto 15271                                             2106
      al=al*alf                                                            2106
      goto 15251                                                           2107
15271 if(ilm .ne. 1)goto 15281                                             2107
      al=big                                                               2107
      goto 15291                                                           2108
15281 continue                                                             2108
      al0=0.0                                                              2109
15300 do 15301 j=1,ni                                                      2109
      if(ju(j).eq.0)goto 15301                                             2109
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2109
15301 continue                                                             2110
15302 continue                                                             2110
      al0=al0/max(bta,1.0e-3)                                              2110
      al=alf*al0                                                           2111
15291 continue                                                             2112
15251 continue                                                             2112
      al2=al*omb                                                           2112
      al1=al*bta                                                           2112
      tlam=bta*(2.0*al-al0)                                                2113
15310 do 15311 k=1,ni                                                      2113
      if(ixx(k).eq.1)goto 15311                                            2113
      if(ju(k).eq.0)goto 15311                                             2114
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2115
15311 continue                                                             2116
15312 continue                                                             2116
10880 continue                                                             2117
15320 continue                                                             2117
15321 continue                                                             2117
      bs(0)=b(0)                                                           2117
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                2118
15330 do 15331 j=1,ni                                                      2118
      if(ixx(j).eq.0)goto 15331                                            2119
      jb=ix(j)                                                             2119
      je=ix(j+1)-1                                                         2119
      jn=ix(j+1)-ix(j)                                                     2120
      sc(1:jn)=v(jx(jb:je))                                                2121
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 2122
      if(kopt .ne. 0)goto 15351                                            2123
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              2124
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                2125
15351 continue                                                             2126
15331 continue                                                             2127
15332 continue                                                             2127
15360 continue                                                             2127
15361 continue                                                             2127
      nlp=nlp+1                                                            2127
      dlx=0.0                                                              2128
15370 do 15371 k=1,ni                                                      2128
      if(ixx(k).eq.0)goto 15371                                            2129
      jb=ix(k)                                                             2129
      je=ix(k+1)-1                                                         2129
      jn=ix(k+1)-ix(k)                                                     2129
      bk=b(k)                                                              2130
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2131
      gk=dot_product(sc(1:jn),x(jb:je))                                    2132
      gk=(gk-svr*xb(k))/xs(k)                                              2133
      u=gk+xv(k)*b(k)                                                      2133
      au=abs(u)-vp(k)*al1                                                  2134
      if(au .gt. 0.0)goto 15391                                            2134
      b(k)=0.0                                                             2134
      goto 15401                                                           2135
15391 continue                                                             2136
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2137
15401 continue                                                             2138
15381 continue                                                             2138
      d=b(k)-bk                                                            2138
      if(abs(d).le.0.0)goto 15371                                          2138
      dlx=max(dlx,xv(k)*d**2)                                              2139
      if(mm(k) .ne. 0)goto 15421                                           2139
      nin=nin+1                                                            2139
      if(nin.gt.nx)goto 15372                                              2140
      mm(k)=nin                                                            2140
      m(nin)=k                                                             2140
      sc(1:jn)=v(jx(jb:je))                                                2141
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 2142
15421 continue                                                             2143
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2144
      o=o+d*(xb(k)/xs(k))                                                  2145
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2146
15371 continue                                                             2147
15372 continue                                                             2147
      if(nin.gt.nx)goto 15362                                              2148
      d=0.0                                                                2148
      if(intr.ne.0) d=svr/xm(0)                                            2149
      if(d .eq. 0.0)goto 15441                                             2149
      b(0)=b(0)+d                                                          2149
      dlx=max(dlx,xm(0)*d**2)                                              2149
      r=r-d*v                                                              2150
      svr=svr-d*xm(0)                                                      2151
15441 continue                                                             2152
      if(dlx.lt.shr)goto 15362                                             2153
      if(nlp .le. maxit)goto 15461                                         2153
      jerr=-ilm                                                            2153
      return                                                               2153
15461 continue                                                             2154
15470 continue                                                             2154
15471 continue                                                             2154
      nlp=nlp+1                                                            2154
      dlx=0.0                                                              2155
15480 do 15481 l=1,nin                                                     2155
      k=m(l)                                                               2155
      jb=ix(k)                                                             2155
      je=ix(k+1)-1                                                         2156
      jn=ix(k+1)-ix(k)                                                     2156
      bk=b(k)                                                              2157
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2158
      gk=dot_product(sc(1:jn),x(jb:je))                                    2159
      gk=(gk-svr*xb(k))/xs(k)                                              2160
      u=gk+xv(k)*b(k)                                                      2160
      au=abs(u)-vp(k)*al1                                                  2161
      if(au .gt. 0.0)goto 15501                                            2161
      b(k)=0.0                                                             2161
      goto 15511                                                           2162
15501 continue                                                             2163
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2164
15511 continue                                                             2165
15491 continue                                                             2165
      d=b(k)-bk                                                            2165
      if(abs(d).le.0.0)goto 15481                                          2165
      dlx=max(dlx,xv(k)*d**2)                                              2166
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2167
      o=o+d*(xb(k)/xs(k))                                                  2168
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2169
15481 continue                                                             2170
15482 continue                                                             2170
      d=0.0                                                                2170
      if(intr.ne.0) d=svr/xm(0)                                            2171
      if(d .eq. 0.0)goto 15531                                             2171
      b(0)=b(0)+d                                                          2171
      dlx=max(dlx,xm(0)*d**2)                                              2171
      r=r-d*v                                                              2172
      svr=svr-d*xm(0)                                                      2173
15531 continue                                                             2174
      if(dlx.lt.shr)goto 15472                                             2175
      if(nlp .le. maxit)goto 15551                                         2175
      jerr=-ilm                                                            2175
      return                                                               2175
15551 continue                                                             2176
      goto 15471                                                           2177
15472 continue                                                             2177
      goto 15361                                                           2178
15362 continue                                                             2178
      if(nin.gt.nx)goto 15322                                              2179
      sc=b(0)                                                              2179
      b0=0.0                                                               2180
15560 do 15561 j=1,nin                                                     2180
      l=m(j)                                                               2180
      jb=ix(l)                                                             2180
      je=ix(l+1)-1                                                         2181
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      2182
      b0=b0-b(l)*xb(l)/xs(l)                                               2183
15561 continue                                                             2184
15562 continue                                                             2184
      sc=sc+b0                                                             2185
15570 do 15571 i=1,no                                                      2185
      fi=sc(i)+g(i)                                                        2186
      if(fi .ge. fmin)goto 15591                                           2186
      q(i)=0.0                                                             2186
      goto 15581                                                           2186
15591 if(fi .le. fmax)goto 15601                                           2186
      q(i)=1.0                                                             2186
      goto 15611                                                           2187
15601 continue                                                             2187
      q(i)=1.0/(1.0+exp(-fi))                                              2187
15611 continue                                                             2188
15581 continue                                                             2188
15571 continue                                                             2189
15572 continue                                                             2189
      v=w*q*(1.0-q)                                                        2189
      xm(0)=sum(v)                                                         2189
      if(xm(0).lt.vmin)goto 15322                                          2190
      r=w*(y-q)                                                            2190
      svr=sum(r)                                                           2190
      o=0.0                                                                2191
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15631                         2191
      kx=0                                                                 2192
15640 do 15641 j=1,nin                                                     2192
      k=m(j)                                                               2193
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15641                           2193
      kx=1                                                                 2193
      goto 15642                                                           2194
15641 continue                                                             2195
15642 continue                                                             2195
      if(kx .ne. 0)goto 15661                                              2196
15670 do 15671 j=1,ni                                                      2196
      if(ixx(j).eq.1)goto 15671                                            2196
      if(ju(j).eq.0)goto 15671                                             2197
      jb=ix(j)                                                             2197
      je=ix(j+1)-1                                                         2197
      jn=ix(j+1)-ix(j)                                                     2198
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2199
      gj=dot_product(sc(1:jn),x(jb:je))                                    2200
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2201
      if(ga(j) .le. al1*vp(j))goto 15691                                   2201
      ixx(j)=1                                                             2201
      kx=1                                                                 2201
15691 continue                                                             2202
15671 continue                                                             2203
15672 continue                                                             2203
      if(kx.eq.1) go to 10880                                              2204
      goto 15322                                                           2205
15661 continue                                                             2206
15631 continue                                                             2207
      goto 15321                                                           2208
15322 continue                                                             2208
      if(nin .le. nx)goto 15711                                            2208
      jerr=-10000-ilm                                                      2208
      goto 15242                                                           2208
15711 continue                                                             2209
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2209
      kin(ilm)=nin                                                         2210
      a0(ilm)=b(0)                                                         2210
      alm(ilm)=al                                                          2210
      lmu=ilm                                                              2211
      devi=dev2(no,w,y,q,pmin)                                             2212
      dev(ilm)=(dev1-devi)/dev0                                            2213
      if(ilm.lt.mnl)goto 15241                                             2213
      if(flmin.ge.1.0)goto 15241                                           2214
      me=0                                                                 2214
15720 do 15721 j=1,nin                                                     2214
      if(a(j,ilm).ne.0.0) me=me+1                                          2214
15721 continue                                                             2214
15722 continue                                                             2214
      if(me.gt.ne)goto 15242                                               2215
      if(dev(ilm).gt.devmax)goto 15242                                     2215
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15242                             2216
      if(xm(0).lt.vmin)goto 15242                                          2217
15241 continue                                                             2218
15242 continue                                                             2218
      g=log(q/(1.0-q))                                                     2219
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2220
      return                                                               2221
      end                                                                  2222
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   2224 
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    2225
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          2226
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2227
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2238
      exmn=-exmx                                                           2239
      allocate(xm(0:ni),stat=ierr)                                         2239
      jerr=jerr+ierr                                                       2240
      allocate(r(1:no),stat=ierr)                                          2240
      jerr=jerr+ierr                                                       2241
      allocate(v(1:no),stat=ierr)                                          2241
      jerr=jerr+ierr                                                       2242
      allocate(mm(1:ni),stat=ierr)                                         2242
      jerr=jerr+ierr                                                       2243
      allocate(ga(1:ni),stat=ierr)                                         2243
      jerr=jerr+ierr                                                       2244
      allocate(iy(1:ni),stat=ierr)                                         2244
      jerr=jerr+ierr                                                       2245
      allocate(is(1:max(nc,ni)),stat=ierr)                                 2245
      jerr=jerr+ierr                                                       2246
      allocate(sxp(1:no),stat=ierr)                                        2246
      jerr=jerr+ierr                                                       2247
      allocate(sxpl(1:no),stat=ierr)                                       2247
      jerr=jerr+ierr                                                       2248
      allocate(sc(1:no),stat=ierr)                                         2248
      jerr=jerr+ierr                                                       2249
      if(jerr.ne.0) return                                                 2250
      pmax=1.0-pmin                                                        2250
      emin=pmin/pmax                                                       2250
      emax=1.0/emin                                                        2251
      pfm=(1.0+pmin)*pmin                                                  2251
      pfx=(1.0-pmin)*pmax                                                  2251
      vmin=pfm*pmax                                                        2252
      bta=parm                                                             2252
      omb=1.0-bta                                                          2252
      dev1=0.0                                                             2252
      dev0=0.0                                                             2253
15730 do 15731 ic=1,nc                                                     2253
      q0=dot_product(w,y(:,ic))                                            2254
      if(q0 .gt. pmin)goto 15751                                           2254
      jerr =8000+ic                                                        2254
      return                                                               2254
15751 continue                                                             2255
      if(q0 .lt. 1.0-pmin)goto 15771                                       2255
      jerr =9000+ic                                                        2255
      return                                                               2255
15771 continue                                                             2256
      if(intr.eq.0) q0=1.0/nc                                              2257
      b(1:ni,ic)=0.0                                                       2257
      b(0,ic)=0.0                                                          2258
      if(intr .eq. 0)goto 15791                                            2258
      b(0,ic)=log(q0)                                                      2258
      dev1=dev1-q0*b(0,ic)                                                 2258
15791 continue                                                             2259
15731 continue                                                             2260
15732 continue                                                             2260
      if(intr.eq.0) dev1=log(float(nc))                                    2260
      iy=0                                                                 2260
      al=0.0                                                               2261
      if(nonzero(no*nc,g) .ne. 0)goto 15811                                2262
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2262
      sxp=0.0                                                              2263
15820 do 15821 ic=1,nc                                                     2263
      q(:,ic)=exp(b(0,ic))                                                 2263
      sxp=sxp+q(:,ic)                                                      2263
15821 continue                                                             2264
15822 continue                                                             2264
      goto 15831                                                           2265
15811 continue                                                             2265
15840 do 15841 i=1,no                                                      2265
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2265
15841 continue                                                             2265
15842 continue                                                             2265
      sxp=0.0                                                              2266
      if(intr .ne. 0)goto 15861                                            2266
      b(0,:)=0.0                                                           2266
      goto 15871                                                           2267
15861 continue                                                             2267
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2267
      if(jerr.ne.0) return                                                 2267
15871 continue                                                             2268
15851 continue                                                             2268
      dev1=0.0                                                             2269
15880 do 15881 ic=1,nc                                                     2269
      q(:,ic)=b(0,ic)+g(:,ic)                                              2270
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2271
      q(:,ic)=exp(q(:,ic))                                                 2271
      sxp=sxp+q(:,ic)                                                      2272
15881 continue                                                             2273
15882 continue                                                             2273
      sxpl=w*log(sxp)                                                      2273
15890 do 15891 ic=1,nc                                                     2273
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2273
15891 continue                                                             2274
15892 continue                                                             2274
15831 continue                                                             2275
15801 continue                                                             2275
15900 do 15901 ic=1,nc                                                     2275
15910 do 15911 i=1,no                                                      2275
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2275
15911 continue                                                             2275
15912 continue                                                             2275
15901 continue                                                             2276
15902 continue                                                             2276
      dev0=dev0+dev1                                                       2277
      if(kopt .le. 0)goto 15931                                            2278
      if(isd .le. 0 .or. intr .eq. 0)goto 15951                            2278
      xv=0.25                                                              2278
      goto 15961                                                           2279
15951 continue                                                             2280
15970 do 15971 j=1,ni                                                      2280
      if(ju(j).eq.0)goto 15971                                             2280
      jb=ix(j)                                                             2280
      je=ix(j+1)-1                                                         2281
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2282
15971 continue                                                             2283
15972 continue                                                             2283
15961 continue                                                             2284
15941 continue                                                             2284
15931 continue                                                             2285
      if(flmin .ge. 1.0)goto 15991                                         2285
      eqs=max(eps,flmin)                                                   2285
      alf=eqs**(1.0/(nlam-1))                                              2285
15991 continue                                                             2286
      m=0                                                                  2286
      mm=0                                                                 2286
      nin=0                                                                2286
      nlp=0                                                                2286
      mnl=min(mnlam,nlam)                                                  2286
      bs=0.0                                                               2286
      svr=0.0                                                              2286
      o=0.0                                                                2287
      shr=shri*dev0                                                        2287
      ga=0.0                                                               2288
16000 do 16001 ic=1,nc                                                     2288
      v=q(:,ic)/sxp                                                        2288
      r=w*(y(:,ic)-v)                                                      2288
      v=w*v*(1.0-v)                                                        2289
16010 do 16011 j=1,ni                                                      2289
      if(ju(j).eq.0)goto 16011                                             2290
      jb=ix(j)                                                             2290
      je=ix(j+1)-1                                                         2290
      jn=ix(j+1)-ix(j)                                                     2291
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2292
      gj=dot_product(sc(1:jn),x(jb:je))                                    2293
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2294
16011 continue                                                             2295
16012 continue                                                             2295
16001 continue                                                             2296
16002 continue                                                             2296
16020 do 16021 ilm=1,nlam                                                  2296
      al0=al                                                               2297
      if(flmin .lt. 1.0)goto 16041                                         2297
      al=ulam(ilm)                                                         2297
      goto 16031                                                           2298
16041 if(ilm .le. 2)goto 16051                                             2298
      al=al*alf                                                            2298
      goto 16031                                                           2299
16051 if(ilm .ne. 1)goto 16061                                             2299
      al=big                                                               2299
      goto 16071                                                           2300
16061 continue                                                             2300
      al0=0.0                                                              2301
16080 do 16081 j=1,ni                                                      2301
      if(ju(j).eq.0)goto 16081                                             2301
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2301
16081 continue                                                             2302
16082 continue                                                             2302
      al0=al0/max(bta,1.0e-3)                                              2302
      al=alf*al0                                                           2303
16071 continue                                                             2304
16031 continue                                                             2304
      al2=al*omb                                                           2304
      al1=al*bta                                                           2304
      tlam=bta*(2.0*al-al0)                                                2305
16090 do 16091 k=1,ni                                                      2305
      if(iy(k).eq.1)goto 16091                                             2305
      if(ju(k).eq.0)goto 16091                                             2306
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2307
16091 continue                                                             2308
16092 continue                                                             2308
10880 continue                                                             2309
16100 continue                                                             2309
16101 continue                                                             2309
      ixx=0                                                                2309
      jxx=ixx                                                              2309
      ig=0                                                                 2310
16110 do 16111 ic=1,nc                                                     2310
      bs(0,ic)=b(0,ic)                                                     2311
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2312
      xm(0)=0.0                                                            2312
      svr=0.0                                                              2312
      o=0.0                                                                2313
16120 do 16121 i=1,no                                                      2313
      pic=q(i,ic)/sxp(i)                                                   2314
      if(pic .ge. pfm)goto 16141                                           2314
      pic=0.0                                                              2314
      v(i)=0.0                                                             2314
      goto 16131                                                           2315
16141 if(pic .le. pfx)goto 16151                                           2315
      pic=1.0                                                              2315
      v(i)=0.0                                                             2315
      goto 16161                                                           2316
16151 continue                                                             2316
      v(i)=w(i)*pic*(1.0-pic)                                              2316
      xm(0)=xm(0)+v(i)                                                     2316
16161 continue                                                             2317
16131 continue                                                             2317
      r(i)=w(i)*(y(i,ic)-pic)                                              2317
      svr=svr+r(i)                                                         2318
16121 continue                                                             2319
16122 continue                                                             2319
      if(xm(0).le.vmin)goto 16111                                          2319
      ig=1                                                                 2320
16170 do 16171 j=1,ni                                                      2320
      if(iy(j).eq.0)goto 16171                                             2321
      jb=ix(j)                                                             2321
      je=ix(j+1)-1                                                         2322
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2323
      if(kopt .ne. 0)goto 16191                                            2324
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2325
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2326
16191 continue                                                             2327
16171 continue                                                             2328
16172 continue                                                             2328
16200 continue                                                             2328
16201 continue                                                             2328
      nlp=nlp+1                                                            2328
      dlx=0.0                                                              2329
16210 do 16211 k=1,ni                                                      2329
      if(iy(k).eq.0)goto 16211                                             2330
      jb=ix(k)                                                             2330
      je=ix(k+1)-1                                                         2330
      jn=ix(k+1)-ix(k)                                                     2330
      bk=b(k,ic)                                                           2331
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2332
      gk=dot_product(sc(1:jn),x(jb:je))                                    2333
      gk=(gk-svr*xb(k))/xs(k)                                              2334
      u=gk+xv(k,ic)*b(k,ic)                                                2334
      au=abs(u)-vp(k)*al1                                                  2335
      if(au .gt. 0.0)goto 16231                                            2335
      b(k,ic)=0.0                                                          2335
      goto 16241                                                           2336
16231 continue                                                             2337
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2339 
     *)
16241 continue                                                             2340
16221 continue                                                             2340
      d=b(k,ic)-bk                                                         2340
      if(abs(d).le.0.0)goto 16211                                          2341
      dlx=max(dlx,xv(k,ic)*d**2)                                           2342
      if(mm(k) .ne. 0)goto 16261                                           2342
      nin=nin+1                                                            2343
      if(nin .le. nx)goto 16281                                            2343
      jxx=1                                                                2343
      goto 16212                                                           2343
16281 continue                                                             2344
      mm(k)=nin                                                            2344
      m(nin)=k                                                             2345
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2346
16261 continue                                                             2347
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2348
      o=o+d*(xb(k)/xs(k))                                                  2349
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2350
16211 continue                                                             2351
16212 continue                                                             2351
      if(jxx.gt.0)goto 16202                                               2352
      d=0.0                                                                2352
      if(intr.ne.0) d=svr/xm(0)                                            2353
      if(d .eq. 0.0)goto 16301                                             2353
      b(0,ic)=b(0,ic)+d                                                    2353
      dlx=max(dlx,xm(0)*d**2)                                              2354
      r=r-d*v                                                              2354
      svr=svr-d*xm(0)                                                      2355
16301 continue                                                             2356
      if(dlx.lt.shr)goto 16202                                             2356
      if(nlp .le. maxit)goto 16321                                         2356
      jerr=-ilm                                                            2356
      return                                                               2356
16321 continue                                                             2357
16330 continue                                                             2357
16331 continue                                                             2357
      nlp=nlp+1                                                            2357
      dlx=0.0                                                              2358
16340 do 16341 l=1,nin                                                     2358
      k=m(l)                                                               2358
      jb=ix(k)                                                             2358
      je=ix(k+1)-1                                                         2359
      jn=ix(k+1)-ix(k)                                                     2359
      bk=b(k,ic)                                                           2360
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2361
      gk=dot_product(sc(1:jn),x(jb:je))                                    2362
      gk=(gk-svr*xb(k))/xs(k)                                              2363
      u=gk+xv(k,ic)*b(k,ic)                                                2363
      au=abs(u)-vp(k)*al1                                                  2364
      if(au .gt. 0.0)goto 16361                                            2364
      b(k,ic)=0.0                                                          2364
      goto 16371                                                           2365
16361 continue                                                             2366
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2368 
     *)
16371 continue                                                             2369
16351 continue                                                             2369
      d=b(k,ic)-bk                                                         2369
      if(abs(d).le.0.0)goto 16341                                          2370
      dlx=max(dlx,xv(k,ic)*d**2)                                           2371
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2372
      o=o+d*(xb(k)/xs(k))                                                  2373
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2374
16341 continue                                                             2375
16342 continue                                                             2375
      d=0.0                                                                2375
      if(intr.ne.0) d=svr/xm(0)                                            2376
      if(d .eq. 0.0)goto 16391                                             2376
      b(0,ic)=b(0,ic)+d                                                    2376
      dlx=max(dlx,xm(0)*d**2)                                              2377
      r=r-d*v                                                              2377
      svr=svr-d*xm(0)                                                      2378
16391 continue                                                             2379
      if(dlx.lt.shr)goto 16332                                             2379
      if(nlp .le. maxit)goto 16411                                         2379
      jerr=-ilm                                                            2379
      return                                                               2379
16411 continue                                                             2380
      goto 16331                                                           2381
16332 continue                                                             2381
      goto 16201                                                           2382
16202 continue                                                             2382
      if(jxx.gt.0)goto 16112                                               2383
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2384
      if(ixx .ne. 0)goto 16431                                             2385
16440 do 16441 j=1,nin                                                     2385
      k=m(j)                                                               2386
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16461                2386
      ixx=1                                                                2386
      goto 16442                                                           2386
16461 continue                                                             2387
16441 continue                                                             2388
16442 continue                                                             2388
16431 continue                                                             2389
      sc=b(0,ic)+g(:,ic)                                                   2389
      b0=0.0                                                               2390
16470 do 16471 j=1,nin                                                     2390
      l=m(j)                                                               2390
      jb=ix(l)                                                             2390
      je=ix(l+1)-1                                                         2391
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2392
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2393
16471 continue                                                             2394
16472 continue                                                             2394
      sc=min(max(exmn,sc+b0),exmx)                                         2395
      sxp=sxp-q(:,ic)                                                      2396
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2397
      sxp=sxp+q(:,ic)                                                      2398
16111 continue                                                             2399
16112 continue                                                             2399
      s=-sum(b(0,:))/nc                                                    2399
      b(0,:)=b(0,:)+s                                                      2399
      sc=s                                                                 2399
      b0=0.0                                                               2400
16480 do 16481 j=1,nin                                                     2400
      l=m(j)                                                               2401
      if(vp(l) .gt. 0.0)goto 16501                                         2401
      s=sum(b(l,:))/nc                                                     2401
      goto 16511                                                           2402
16501 continue                                                             2402
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     2402
16511 continue                                                             2403
16491 continue                                                             2403
      b(l,:)=b(l,:)-s                                                      2404
      jb=ix(l)                                                             2404
      je=ix(l+1)-1                                                         2405
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2406
      b0=b0+s*xb(l)/xs(l)                                                  2407
16481 continue                                                             2408
16482 continue                                                             2408
      sc=sc+b0                                                             2408
      sc=exp(sc)                                                           2408
      sxp=sxp*sc                                                           2408
16520 do 16521 ic=1,nc                                                     2408
      q(:,ic)=q(:,ic)*sc                                                   2408
16521 continue                                                             2409
16522 continue                                                             2409
      if(jxx.gt.0)goto 16102                                               2409
      if(ig.eq.0)goto 16102                                                2410
      if(ixx .ne. 0)goto 16541                                             2411
16550 do 16551 j=1,ni                                                      2411
      if(iy(j).eq.1)goto 16551                                             2411
      if(ju(j).eq.0)goto 16551                                             2411
      ga(j)=0.0                                                            2411
16551 continue                                                             2412
16552 continue                                                             2412
16560 do 16561 ic=1,nc                                                     2412
      v=q(:,ic)/sxp                                                        2412
      r=w*(y(:,ic)-v)                                                      2412
      v=w*v*(1.0-v)                                                        2413
16570 do 16571 j=1,ni                                                      2413
      if(iy(j).eq.1)goto 16571                                             2413
      if(ju(j).eq.0)goto 16571                                             2414
      jb=ix(j)                                                             2414
      je=ix(j+1)-1                                                         2414
      jn=ix(j+1)-ix(j)                                                     2415
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2416
      gj=dot_product(sc(1:jn),x(jb:je))                                    2417
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2418
16571 continue                                                             2419
16572 continue                                                             2419
16561 continue                                                             2420
16562 continue                                                             2420
16580 do 16581 k=1,ni                                                      2420
      if(iy(k).eq.1)goto 16581                                             2420
      if(ju(k).eq.0)goto 16581                                             2421
      if(ga(k) .le. al1*vp(k))goto 16601                                   2421
      iy(k)=1                                                              2421
      ixx=1                                                                2421
16601 continue                                                             2422
16581 continue                                                             2423
16582 continue                                                             2423
      if(ixx.eq.1) go to 10880                                             2424
      goto 16102                                                           2425
16541 continue                                                             2426
      goto 16101                                                           2427
16102 continue                                                             2427
      if(jxx .le. 0)goto 16621                                             2427
      jerr=-10000-ilm                                                      2427
      goto 16022                                                           2427
16621 continue                                                             2427
      devi=0.0                                                             2428
16630 do 16631 ic=1,nc                                                     2429
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2429
      a0(ic,ilm)=b(0,ic)                                                   2430
16640 do 16641 i=1,no                                                      2430
      if(y(i,ic).le.0.0)goto 16641                                         2431
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2432
16641 continue                                                             2433
16642 continue                                                             2433
16631 continue                                                             2434
16632 continue                                                             2434
      kin(ilm)=nin                                                         2434
      alm(ilm)=al                                                          2434
      lmu=ilm                                                              2435
      dev(ilm)=(dev1-devi)/dev0                                            2435
      if(ig.eq.0)goto 16022                                                2436
      if(ilm.lt.mnl)goto 16021                                             2436
      if(flmin.ge.1.0)goto 16021                                           2437
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16022             2438
      if(dev(ilm).gt.devmax)goto 16022                                     2438
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16022                             2439
16021 continue                                                             2440
16022 continue                                                             2440
      g=log(q)                                                             2440
16650 do 16651 i=1,no                                                      2440
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2440
16651 continue                                                             2441
16652 continue                                                             2441
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2442
      return                                                               2443
      end                                                                  2444
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2445
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2445
      integer ia(*),ix(*),jx(*)                                            2446
16660 do 16661 ic=1,nc                                                     2446
      f(ic,:)=a0(ic)                                                       2446
16661 continue                                                             2447
16662 continue                                                             2447
16670 do 16671 j=1,nin                                                     2447
      k=ia(j)                                                              2447
      kb=ix(k)                                                             2447
      ke=ix(k+1)-1                                                         2448
16680 do 16681 ic=1,nc                                                     2448
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2448
16681 continue                                                             2449
16682 continue                                                             2449
16671 continue                                                             2450
16672 continue                                                             2450
      return                                                               2451
      end                                                                  2452
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   2454 
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2455
      real ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)                        2456
      integer jd(*),ia(nx),nin(nlam)                                       2457
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16701                                    2461
      jerr=10000                                                           2461
      return                                                               2461
16701 continue                                                             2462
      allocate(ww(1:no),stat=jerr)                                         2463
      allocate(ju(1:ni),stat=ierr)                                         2463
      jerr=jerr+ierr                                                       2464
      allocate(vq(1:ni),stat=ierr)                                         2464
      jerr=jerr+ierr                                                       2465
      if(isd .le. 0)goto 16721                                             2465
      allocate(xs(1:ni),stat=ierr)                                         2465
      jerr=jerr+ierr                                                       2465
16721 continue                                                             2466
      if(jerr.ne.0) return                                                 2467
      call chkvars(no,ni,x,ju)                                             2468
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2469
      if(maxval(ju) .gt. 0)goto 16741                                      2469
      jerr=7777                                                            2469
      return                                                               2469
16741 continue                                                             2470
      vq=max(0.0,vp)                                                       2470
      vq=vq*ni/sum(vq)                                                     2471
      ww=max(0.0,w)                                                        2471
      sw=sum(ww)                                                           2472
      if(sw .gt. 0.0)goto 16761                                            2472
      jerr=9999                                                            2472
      return                                                               2472
16761 continue                                                             2472
      ww=ww/sw                                                             2473
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2474
      if(isd .le. 0)goto 16781                                             2474
16790 do 16791 j=1,ni                                                      2474
      cl(:,j)=cl(:,j)*xs(j)                                                2474
16791 continue                                                             2474
16792 continue                                                             2474
16781 continue                                                             2475
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   2477 
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2477
      dev0=2.0*sw*dev0                                                     2478
      if(isd .le. 0)goto 16811                                             2478
16820 do 16821 k=1,lmu                                                     2478
      nk=nin(k)                                                            2478
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2478
16821 continue                                                             2478
16822 continue                                                             2478
16811 continue                                                             2479
      deallocate(ww,ju,vq)                                                 2479
      if(isd.gt.0) deallocate(xs)                                          2480
      return                                                               2481
      end                                                                  2482
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2483
      real x(no,ni),w(no),xs(ni)                                           2483
      integer ju(ni)                                                       2484
16830 do 16831 j=1,ni                                                      2484
      if(ju(j).eq.0)goto 16831                                             2485
      xm=dot_product(w,x(:,j))                                             2485
      x(:,j)=x(:,j)-xm                                                     2486
      if(isd .le. 0)goto 16851                                             2486
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2486
      x(:,j)=x(:,j)/xs(j)                                                  2486
16851 continue                                                             2487
16831 continue                                                             2488
16832 continue                                                             2488
      return                                                               2489
      end                                                                  2490
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   2492 
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2493
      real ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)                        2494
      integer ju(ni),m(nx),kin(nlam)                                       2495
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2501
      sml=sml*100.0                                                        2501
      devmax=devmax*0.99/0.999                                             2502
      allocate(e(1:no),stat=jerr)                                          2503
      allocate(uu(1:no),stat=ierr)                                         2503
      jerr=jerr+ierr                                                       2504
      allocate(f(1:no),stat=ierr)                                          2504
      jerr=jerr+ierr                                                       2505
      allocate(w(1:no),stat=ierr)                                          2505
      jerr=jerr+ierr                                                       2506
      allocate(v(1:ni),stat=ierr)                                          2506
      jerr=jerr+ierr                                                       2507
      allocate(a(1:ni),stat=ierr)                                          2507
      jerr=jerr+ierr                                                       2508
      allocate(as(1:ni),stat=ierr)                                         2508
      jerr=jerr+ierr                                                       2509
      allocate(xs(1:ni),stat=ierr)                                         2509
      jerr=jerr+ierr                                                       2510
      allocate(ga(1:ni),stat=ierr)                                         2510
      jerr=jerr+ierr                                                       2511
      allocate(ixx(1:ni),stat=ierr)                                        2511
      jerr=jerr+ierr                                                       2512
      allocate(jp(1:no),stat=ierr)                                         2512
      jerr=jerr+ierr                                                       2513
      allocate(kp(1:no),stat=ierr)                                         2513
      jerr=jerr+ierr                                                       2514
      allocate(dk(1:no),stat=ierr)                                         2514
      jerr=jerr+ierr                                                       2515
      allocate(wr(1:no),stat=ierr)                                         2515
      jerr=jerr+ierr                                                       2516
      allocate(dq(1:no),stat=ierr)                                         2516
      jerr=jerr+ierr                                                       2517
      allocate(mm(1:ni),stat=ierr)                                         2517
      jerr=jerr+ierr                                                       2518
      if(jerr.ne.0)go to 12180                                             2519
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2520
      if(jerr.ne.0) go to 12180                                            2520
      alpha=parm                                                           2521
      oma=1.0-alpha                                                        2521
      nlm=0                                                                2521
      ixx=0                                                                2521
      al=0.0                                                               2522
      dq=d*q                                                               2522
      call died(no,nk,dq,kp,jp,dk)                                         2523
      a=0.0                                                                2523
      f(1)=0.0                                                             2523
      fmax=log(huge(f(1))*0.1)                                             2524
      if(nonzero(no,g) .eq. 0)goto 16871                                   2524
      f=g-dot_product(q,g)                                                 2525
      e=q*exp(sign(min(abs(f),fmax),f))                                    2526
      goto 16881                                                           2527
16871 continue                                                             2527
      f=0.0                                                                2527
      e=q                                                                  2527
16881 continue                                                             2528
16861 continue                                                             2528
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2529
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2529
      dev0=rr                                                              2530
16890 do 16891 i=1,no                                                      2530
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 16911                   2530
      w(i)=0.0                                                             2530
      wr(i)=w(i)                                                           2530
16911 continue                                                             2530
16891 continue                                                             2531
16892 continue                                                             2531
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2532
      if(jerr.ne.0) go to 12180                                            2533
      if(flmin .ge. 1.0)goto 16931                                         2533
      eqs=max(eps,flmin)                                                   2533
      alf=eqs**(1.0/(nlam-1))                                              2533
16931 continue                                                             2534
      m=0                                                                  2534
      mm=0                                                                 2534
      nlp=0                                                                2534
      nin=nlp                                                              2534
      mnl=min(mnlam,nlam)                                                  2534
      as=0.0                                                               2534
      cthr=cthri*dev0                                                      2535
16940 do 16941 j=1,ni                                                      2535
      if(ju(j).eq.0)goto 16941                                             2535
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2535
16941 continue                                                             2536
16942 continue                                                             2536
16950 do 16951 ilm=1,nlam                                                  2536
      al0=al                                                               2537
      if(flmin .lt. 1.0)goto 16971                                         2537
      al=ulam(ilm)                                                         2537
      goto 16961                                                           2538
16971 if(ilm .le. 2)goto 16981                                             2538
      al=al*alf                                                            2538
      goto 16961                                                           2539
16981 if(ilm .ne. 1)goto 16991                                             2539
      al=big                                                               2539
      goto 17001                                                           2540
16991 continue                                                             2540
      al0=0.0                                                              2541
17010 do 17011 j=1,ni                                                      2541
      if(ju(j).eq.0)goto 17011                                             2541
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2541
17011 continue                                                             2542
17012 continue                                                             2542
      al0=al0/max(parm,1.0e-3)                                             2542
      al=alf*al0                                                           2543
17001 continue                                                             2544
16961 continue                                                             2544
      sa=alpha*al                                                          2544
      omal=oma*al                                                          2544
      tlam=alpha*(2.0*al-al0)                                              2545
17020 do 17021 k=1,ni                                                      2545
      if(ixx(k).eq.1)goto 17021                                            2545
      if(ju(k).eq.0)goto 17021                                             2546
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2547
17021 continue                                                             2548
17022 continue                                                             2548
10880 continue                                                             2549
17030 continue                                                             2549
17031 continue                                                             2549
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2550
      call vars(no,ni,x,w,ixx,v)                                           2551
17040 continue                                                             2551
17041 continue                                                             2551
      nlp=nlp+1                                                            2551
      dli=0.0                                                              2552
17050 do 17051 j=1,ni                                                      2552
      if(ixx(j).eq.0)goto 17051                                            2553
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2554
      if(abs(u) .gt. vp(j)*sa)goto 17071                                   2554
      at=0.0                                                               2554
      goto 17081                                                           2555
17071 continue                                                             2555
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2557 
     *mal)))
17081 continue                                                             2558
17061 continue                                                             2558
      if(at .eq. a(j))goto 17101                                           2558
      del=at-a(j)                                                          2558
      a(j)=at                                                              2558
      dli=max(dli,v(j)*del**2)                                             2559
      wr=wr-del*w*x(:,j)                                                   2559
      f=f+del*x(:,j)                                                       2560
      if(mm(j) .ne. 0)goto 17121                                           2560
      nin=nin+1                                                            2560
      if(nin.gt.nx)goto 17052                                              2561
      mm(j)=nin                                                            2561
      m(nin)=j                                                             2562
17121 continue                                                             2563
17101 continue                                                             2564
17051 continue                                                             2565
17052 continue                                                             2565
      if(nin.gt.nx)goto 17042                                              2565
      if(dli.lt.cthr)goto 17042                                            2566
      if(nlp .le. maxit)goto 17141                                         2566
      jerr=-ilm                                                            2566
      return                                                               2566
17141 continue                                                             2567
17150 continue                                                             2567
17151 continue                                                             2567
      nlp=nlp+1                                                            2567
      dli=0.0                                                              2568
17160 do 17161 l=1,nin                                                     2568
      j=m(l)                                                               2569
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2570
      if(abs(u) .gt. vp(j)*sa)goto 17181                                   2570
      at=0.0                                                               2570
      goto 17191                                                           2571
17181 continue                                                             2571
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2573 
     *mal)))
17191 continue                                                             2574
17171 continue                                                             2574
      if(at .eq. a(j))goto 17211                                           2574
      del=at-a(j)                                                          2574
      a(j)=at                                                              2574
      dli=max(dli,v(j)*del**2)                                             2575
      wr=wr-del*w*x(:,j)                                                   2575
      f=f+del*x(:,j)                                                       2576
17211 continue                                                             2577
17161 continue                                                             2578
17162 continue                                                             2578
      if(dli.lt.cthr)goto 17152                                            2578
      if(nlp .le. maxit)goto 17231                                         2578
      jerr=-ilm                                                            2578
      return                                                               2578
17231 continue                                                             2579
      goto 17151                                                           2580
17152 continue                                                             2580
      goto 17041                                                           2581
17042 continue                                                             2581
      if(nin.gt.nx)goto 17032                                              2582
      e=q*exp(sign(min(abs(f),fmax),f))                                    2583
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2584
      if(jerr .eq. 0)goto 17251                                            2584
      jerr=jerr-ilm                                                        2584
      go to 12180                                                          2584
17251 continue                                                             2585
      ix=0                                                                 2586
17260 do 17261 j=1,nin                                                     2586
      k=m(j)                                                               2587
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17261                           2587
      ix=1                                                                 2587
      goto 17262                                                           2587
17261 continue                                                             2588
17262 continue                                                             2588
      if(ix .ne. 0)goto 17281                                              2589
17290 do 17291 k=1,ni                                                      2589
      if(ixx(k).eq.1)goto 17291                                            2589
      if(ju(k).eq.0)goto 17291                                             2590
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2591
      if(ga(k) .le. sa*vp(k))goto 17311                                    2591
      ixx(k)=1                                                             2591
      ix=1                                                                 2591
17311 continue                                                             2592
17291 continue                                                             2593
17292 continue                                                             2593
      if(ix.eq.1) go to 10880                                              2594
      goto 17032                                                           2595
17281 continue                                                             2596
      goto 17031                                                           2597
17032 continue                                                             2597
      if(nin .le. nx)goto 17331                                            2597
      jerr=-10000-ilm                                                      2597
      goto 16952                                                           2597
17331 continue                                                             2598
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2598
      kin(ilm)=nin                                                         2599
      alm(ilm)=al                                                          2599
      lmu=ilm                                                              2600
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2601
      if(ilm.lt.mnl)goto 16951                                             2601
      if(flmin.ge.1.0)goto 16951                                           2602
      me=0                                                                 2602
17340 do 17341 j=1,nin                                                     2602
      if(ao(j,ilm).ne.0.0) me=me+1                                         2602
17341 continue                                                             2602
17342 continue                                                             2602
      if(me.gt.ne)goto 16952                                               2603
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16952              2604
      if(dev(ilm).gt.devmax)goto 16952                                     2605
16951 continue                                                             2606
16952 continue                                                             2606
      g=f                                                                  2607
12180 continue                                                             2607
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2608
      return                                                               2609
      end                                                                  2610
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2611
      real ca(nin),x(n,*),f(n)                                             2611
      integer ia(nin)                                                      2612
      f=0.0                                                                2612
      if(nin.le.0) return                                                  2613
17350 do 17351 i=1,n                                                       2613
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2613
17351 continue                                                             2614
17352 continue                                                             2614
      return                                                               2615
      end                                                                  2616
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2617
      real y(no),d(no),q(no)                                               2617
      integer jp(no),kp(*)                                                 2618
17360 do 17361 j=1,no                                                      2618
      jp(j)=j                                                              2618
17361 continue                                                             2618
17362 continue                                                             2618
      call psort7(y,jp,1,no)                                               2619
      nj=0                                                                 2619
17370 do 17371 j=1,no                                                      2619
      if(q(jp(j)).le.0.0)goto 17371                                        2619
      nj=nj+1                                                              2619
      jp(nj)=jp(j)                                                         2619
17371 continue                                                             2620
17372 continue                                                             2620
      if(nj .ne. 0)goto 17391                                              2620
      jerr=20000                                                           2620
      return                                                               2620
17391 continue                                                             2621
      j=1                                                                  2621
17400 continue                                                             2621
17401 if(d(jp(j)).gt.0.0)goto 17402                                        2621
      j=j+1                                                                2621
      if(j.gt.nj)goto 17402                                                2621
      goto 17401                                                           2622
17402 continue                                                             2622
      if(j .lt. nj-1)goto 17421                                            2622
      jerr=30000                                                           2622
      return                                                               2622
17421 continue                                                             2623
      t0=y(jp(j))                                                          2623
      j0=j-1                                                               2624
      if(j0 .le. 0)goto 17441                                              2625
17450 continue                                                             2625
17451 if(y(jp(j0)).lt.t0)goto 17452                                        2625
      j0=j0-1                                                              2625
      if(j0.eq.0)goto 17452                                                2625
      goto 17451                                                           2626
17452 continue                                                             2626
      if(j0 .le. 0)goto 17471                                              2626
      nj=nj-j0                                                             2626
17480 do 17481 j=1,nj                                                      2626
      jp(j)=jp(j+j0)                                                       2626
17481 continue                                                             2626
17482 continue                                                             2626
17471 continue                                                             2627
17441 continue                                                             2628
      jerr=0                                                               2628
      nk=0                                                                 2628
      yk=t0                                                                2628
      j=2                                                                  2629
17490 continue                                                             2629
17491 continue                                                             2629
17500 continue                                                             2630
17501 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17502                     2630
      j=j+1                                                                2630
      if(j.gt.nj)goto 17502                                                2630
      goto 17501                                                           2631
17502 continue                                                             2631
      nk=nk+1                                                              2631
      kp(nk)=j-1                                                           2631
      if(j.gt.nj)goto 17492                                                2632
      if(j .ne. nj)goto 17521                                              2632
      nk=nk+1                                                              2632
      kp(nk)=nj                                                            2632
      goto 17492                                                           2632
17521 continue                                                             2633
      yk=y(jp(j))                                                          2633
      j=j+1                                                                2634
      goto 17491                                                           2635
17492 continue                                                             2635
      return                                                               2636
      end                                                                  2637
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2638
      real d(no),dk(nk),wr(no),w(no)                                       2639
      real e(no),u(no),b,c                                                 2639
      integer kp(nk),jp(no)                                                2640
      call usk(no,nk,kp,jp,e,u)                                            2641
      b=dk(1)/u(1)                                                         2641
      c=dk(1)/u(1)**2                                                      2641
      jerr=0                                                               2642
17530 do 17531 j=1,kp(1)                                                   2642
      i=jp(j)                                                              2643
      w(i)=e(i)*(b-e(i)*c)                                                 2643
      if(w(i) .gt. 0.0)goto 17551                                          2643
      jerr=-30000                                                          2643
      return                                                               2643
17551 continue                                                             2644
      wr(i)=d(i)-e(i)*b                                                    2645
17531 continue                                                             2646
17532 continue                                                             2646
17560 do 17561 k=2,nk                                                      2646
      j1=kp(k-1)+1                                                         2646
      j2=kp(k)                                                             2647
      b=b+dk(k)/u(k)                                                       2647
      c=c+dk(k)/u(k)**2                                                    2648
17570 do 17571 j=j1,j2                                                     2648
      i=jp(j)                                                              2649
      w(i)=e(i)*(b-e(i)*c)                                                 2649
      if(w(i) .gt. 0.0)goto 17591                                          2649
      jerr=-30000                                                          2649
      return                                                               2649
17591 continue                                                             2650
      wr(i)=d(i)-e(i)*b                                                    2651
17571 continue                                                             2652
17572 continue                                                             2652
17561 continue                                                             2653
17562 continue                                                             2653
      return                                                               2654
      end                                                                  2655
      subroutine vars(no,ni,x,w,ixx,v)                                     2656
      real x(no,ni),w(no),v(ni)                                            2656
      integer ixx(ni)                                                      2657
17600 do 17601 j=1,ni                                                      2657
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2657
17601 continue                                                             2658
17602 continue                                                             2658
      return                                                               2659
      end                                                                  2660
      subroutine died(no,nk,d,kp,jp,dk)                                    2661
      real d(no),dk(nk)                                                    2661
      integer kp(nk),jp(no)                                                2662
      dk(1)=sum(d(jp(1:kp(1))))                                            2663
17610 do 17611 k=2,nk                                                      2663
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2663
17611 continue                                                             2664
17612 continue                                                             2664
      return                                                               2665
      end                                                                  2666
      subroutine usk(no,nk,kp,jp,e,u)                                      2667
      real e(no),u(nk),h                                                   2667
      integer kp(nk),jp(no)                                                2668
      h=0.0                                                                2669
17620 do 17621 k=nk,1,-1                                                   2669
      j2=kp(k)                                                             2670
      j1=1                                                                 2670
      if(k.gt.1) j1=kp(k-1)+1                                              2671
17630 do 17631 j=j2,j1,-1                                                  2671
      h=h+e(jp(j))                                                         2671
17631 continue                                                             2672
17632 continue                                                             2672
      u(k)=h                                                               2673
17621 continue                                                             2674
17622 continue                                                             2674
      return                                                               2675
      end                                                                  2676
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2677
      real d(no),dk(nk),f(no)                                              2678
      integer kp(nk),jp(no)                                                2678
      real e(no),u(nk),s                                                   2679
      call usk(no,nk,kp,jp,e,u)                                            2679
      u=log(u)                                                             2680
      risk=dot_product(d,f)-dot_product(dk,u)                              2681
      return                                                               2682
      end                                                                  2683
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2684
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2685
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2691
      allocate(q(1:no),stat=ierr)                                          2691
      jerr=jerr+ierr                                                       2692
      allocate(uu(1:no),stat=ierr)                                         2692
      jerr=jerr+ierr                                                       2693
      allocate(f(1:no),stat=ierr)                                          2693
      jerr=jerr+ierr                                                       2694
      allocate(dk(1:no),stat=ierr)                                         2694
      jerr=jerr+ierr                                                       2695
      allocate(jp(1:no),stat=ierr)                                         2695
      jerr=jerr+ierr                                                       2696
      allocate(kp(1:no),stat=ierr)                                         2696
      jerr=jerr+ierr                                                       2697
      allocate(dq(1:no),stat=ierr)                                         2697
      jerr=jerr+ierr                                                       2698
      allocate(xm(1:ni),stat=ierr)                                         2698
      jerr=jerr+ierr                                                       2699
      if(jerr.ne.0) go to 12180                                            2700
      q=max(0.0,w)                                                         2700
      sw=sum(q)                                                            2701
      if(sw .gt. 0.0)goto 17651                                            2701
      jerr=9999                                                            2701
      go to 12180                                                          2701
17651 continue                                                             2702
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2703
      if(jerr.ne.0) go to 12180                                            2703
      fmax=log(huge(e(1))*0.1)                                             2704
      dq=d*q                                                               2704
      call died(no,nk,dq,kp,jp,dk)                                         2704
      gm=dot_product(q,g)/sw                                               2705
17660 do 17661 j=1,ni                                                      2705
      xm(j)=dot_product(q,x(:,j))/sw                                       2705
17661 continue                                                             2706
17662 continue                                                             2706
17670 do 17671 lam=1,nlam                                                  2707
17680 do 17681 i=1,no                                                      2707
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2708
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2709
17681 continue                                                             2710
17682 continue                                                             2710
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2711
17671 continue                                                             2712
17672 continue                                                             2712
12180 continue                                                             2712
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2713
      return                                                               2714
      end                                                                  2715
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   2717 
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2718
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2719
      integer jd(*),ia(nx),nin(nlam)                                       2720
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17701                                    2724
      jerr=10000                                                           2724
      return                                                               2724
17701 continue                                                             2725
      if(minval(y) .ge. 0.0)goto 17721                                     2725
      jerr=8888                                                            2725
      return                                                               2725
17721 continue                                                             2726
      allocate(ww(1:no),stat=jerr)                                         2727
      allocate(ju(1:ni),stat=ierr)                                         2727
      jerr=jerr+ierr                                                       2728
      allocate(vq(1:ni),stat=ierr)                                         2728
      jerr=jerr+ierr                                                       2729
      allocate(xm(1:ni),stat=ierr)                                         2729
      jerr=jerr+ierr                                                       2730
      if(isd .le. 0)goto 17741                                             2730
      allocate(xs(1:ni),stat=ierr)                                         2730
      jerr=jerr+ierr                                                       2730
17741 continue                                                             2731
      if(jerr.ne.0) return                                                 2732
      call chkvars(no,ni,x,ju)                                             2733
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2734
      if(maxval(ju) .gt. 0)goto 17761                                      2734
      jerr=7777                                                            2734
      go to 12180                                                          2734
17761 continue                                                             2735
      vq=max(0.0,vp)                                                       2735
      vq=vq*ni/sum(vq)                                                     2736
      ww=max(0.0,w)                                                        2736
      sw=sum(ww)                                                           2736
      if(sw .gt. 0.0)goto 17781                                            2736
      jerr=9999                                                            2736
      go to 12180                                                          2736
17781 continue                                                             2737
      ww=ww/sw                                                             2738
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        2739
      if(isd .le. 0)goto 17801                                             2739
17810 do 17811 j=1,ni                                                      2739
      cl(:,j)=cl(:,j)*xs(j)                                                2739
17811 continue                                                             2739
17812 continue                                                             2739
17801 continue                                                             2740
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   2742 
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12180                                            2742
      dev0=2.0*sw*dev0                                                     2743
17820 do 17821 k=1,lmu                                                     2743
      nk=nin(k)                                                            2744
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2745
      if(intr .ne. 0)goto 17841                                            2745
      a0(k)=0.0                                                            2745
      goto 17851                                                           2746
17841 continue                                                             2746
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2746
17851 continue                                                             2747
17831 continue                                                             2747
17821 continue                                                             2748
17822 continue                                                             2748
12180 continue                                                             2748
      deallocate(ww,ju,vq,xm)                                              2748
      if(isd.gt.0) deallocate(xs)                                          2749
      return                                                               2750
      end                                                                  2751
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2753 
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2754
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2755
      integer ju(ni),m(nx),kin(nlam)                                       2756
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2760
      sml=sml*10.0                                                         2761
      allocate(a(1:ni),stat=jerr)                                          2762
      allocate(as(1:ni),stat=ierr)                                         2762
      jerr=jerr+ierr                                                       2763
      allocate(t(1:no),stat=ierr)                                          2763
      jerr=jerr+ierr                                                       2764
      allocate(mm(1:ni),stat=ierr)                                         2764
      jerr=jerr+ierr                                                       2765
      allocate(ga(1:ni),stat=ierr)                                         2765
      jerr=jerr+ierr                                                       2766
      allocate(ixx(1:ni),stat=ierr)                                        2766
      jerr=jerr+ierr                                                       2767
      allocate(wr(1:no),stat=ierr)                                         2767
      jerr=jerr+ierr                                                       2768
      allocate(v(1:ni),stat=ierr)                                          2768
      jerr=jerr+ierr                                                       2769
      allocate(w(1:no),stat=ierr)                                          2769
      jerr=jerr+ierr                                                       2770
      allocate(f(1:no),stat=ierr)                                          2770
      jerr=jerr+ierr                                                       2771
      if(jerr.ne.0) return                                                 2772
      bta=parm                                                             2772
      omb=1.0-bta                                                          2773
      t=q*y                                                                2773
      yb=sum(t)                                                            2773
      fmax=log(huge(bta)*0.1)                                              2774
      if(nonzero(no,g) .ne. 0)goto 17871                                   2775
      if(intr .eq. 0)goto 17891                                            2775
      w=q*yb                                                               2775
      az=log(yb)                                                           2775
      f=az                                                                 2775
      dv0=yb*(az-1.0)                                                      2775
      goto 17901                                                           2776
17891 continue                                                             2776
      w=q                                                                  2776
      az=0.0                                                               2776
      f=az                                                                 2776
      dv0=-1.0                                                             2776
17901 continue                                                             2777
17881 continue                                                             2777
      goto 17911                                                           2778
17871 continue                                                             2778
      w=q*exp(sign(min(abs(g),fmax),g))                                    2778
      v0=sum(w)                                                            2779
      if(intr .eq. 0)goto 17931                                            2779
      eaz=yb/v0                                                            2779
      w=eaz*w                                                              2779
      az=log(eaz)                                                          2779
      f=az+g                                                               2780
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2781
      goto 17941                                                           2782
17931 continue                                                             2782
      az=0.0                                                               2782
      f=g                                                                  2782
      dv0=dot_product(t,g)-v0                                              2782
17941 continue                                                             2783
17921 continue                                                             2783
17911 continue                                                             2784
17861 continue                                                             2784
      a=0.0                                                                2784
      as=0.0                                                               2784
      wr=t-w                                                               2784
      v0=1.0                                                               2784
      if(intr.ne.0) v0=yb                                                  2784
      dvr=-yb                                                              2785
17950 do 17951 i=1,no                                                      2785
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2785
17951 continue                                                             2785
17952 continue                                                             2785
      dvr=dvr-dv0                                                          2785
      dev0=dvr                                                             2786
      if(flmin .ge. 1.0)goto 17971                                         2786
      eqs=max(eps,flmin)                                                   2786
      alf=eqs**(1.0/(nlam-1))                                              2786
17971 continue                                                             2787
      m=0                                                                  2787
      mm=0                                                                 2787
      nlp=0                                                                2787
      nin=nlp                                                              2787
      mnl=min(mnlam,nlam)                                                  2787
      shr=shri*dev0                                                        2787
      ixx=0                                                                2787
      al=0.0                                                               2788
17980 do 17981 j=1,ni                                                      2788
      if(ju(j).eq.0)goto 17981                                             2788
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2788
17981 continue                                                             2789
17982 continue                                                             2789
17990 do 17991 ilm=1,nlam                                                  2789
      al0=al                                                               2790
      if(flmin .lt. 1.0)goto 18011                                         2790
      al=ulam(ilm)                                                         2790
      goto 18001                                                           2791
18011 if(ilm .le. 2)goto 18021                                             2791
      al=al*alf                                                            2791
      goto 18001                                                           2792
18021 if(ilm .ne. 1)goto 18031                                             2792
      al=big                                                               2792
      goto 18041                                                           2793
18031 continue                                                             2793
      al0=0.0                                                              2794
18050 do 18051 j=1,ni                                                      2794
      if(ju(j).eq.0)goto 18051                                             2794
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2794
18051 continue                                                             2795
18052 continue                                                             2795
      al0=al0/max(bta,1.0e-3)                                              2795
      al=alf*al0                                                           2796
18041 continue                                                             2797
18001 continue                                                             2797
      al2=al*omb                                                           2797
      al1=al*bta                                                           2797
      tlam=bta*(2.0*al-al0)                                                2798
18060 do 18061 k=1,ni                                                      2798
      if(ixx(k).eq.1)goto 18061                                            2798
      if(ju(k).eq.0)goto 18061                                             2799
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2800
18061 continue                                                             2801
18062 continue                                                             2801
10880 continue                                                             2802
18070 continue                                                             2802
18071 continue                                                             2802
      az0=az                                                               2803
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2804
18080 do 18081 j=1,ni                                                      2804
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2804
18081 continue                                                             2805
18082 continue                                                             2805
18090 continue                                                             2805
18091 continue                                                             2805
      nlp=nlp+1                                                            2805
      dlx=0.0                                                              2806
18100 do 18101 k=1,ni                                                      2806
      if(ixx(k).eq.0)goto 18101                                            2806
      ak=a(k)                                                              2807
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2807
      au=abs(u)-vp(k)*al1                                                  2808
      if(au .gt. 0.0)goto 18121                                            2808
      a(k)=0.0                                                             2808
      goto 18131                                                           2809
18121 continue                                                             2810
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2811
18131 continue                                                             2812
18111 continue                                                             2812
      if(a(k).eq.ak)goto 18101                                             2812
      d=a(k)-ak                                                            2812
      dlx=max(dlx,v(k)*d**2)                                               2813
      wr=wr-d*w*x(:,k)                                                     2813
      f=f+d*x(:,k)                                                         2814
      if(mm(k) .ne. 0)goto 18151                                           2814
      nin=nin+1                                                            2814
      if(nin.gt.nx)goto 18102                                              2815
      mm(k)=nin                                                            2815
      m(nin)=k                                                             2816
18151 continue                                                             2817
18101 continue                                                             2818
18102 continue                                                             2818
      if(nin.gt.nx)goto 18092                                              2819
      if(intr .eq. 0)goto 18171                                            2819
      d=sum(wr)/v0                                                         2820
      az=az+d                                                              2820
      dlx=max(dlx,v0*d**2)                                                 2820
      wr=wr-d*w                                                            2820
      f=f+d                                                                2821
18171 continue                                                             2822
      if(dlx.lt.shr)goto 18092                                             2822
      if(nlp .le. maxit)goto 18191                                         2822
      jerr=-ilm                                                            2822
      return                                                               2822
18191 continue                                                             2823
18200 continue                                                             2823
18201 continue                                                             2823
      nlp=nlp+1                                                            2823
      dlx=0.0                                                              2824
18210 do 18211 l=1,nin                                                     2824
      k=m(l)                                                               2824
      ak=a(k)                                                              2825
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2825
      au=abs(u)-vp(k)*al1                                                  2826
      if(au .gt. 0.0)goto 18231                                            2826
      a(k)=0.0                                                             2826
      goto 18241                                                           2827
18231 continue                                                             2828
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2829
18241 continue                                                             2830
18221 continue                                                             2830
      if(a(k).eq.ak)goto 18211                                             2830
      d=a(k)-ak                                                            2830
      dlx=max(dlx,v(k)*d**2)                                               2831
      wr=wr-d*w*x(:,k)                                                     2831
      f=f+d*x(:,k)                                                         2833
18211 continue                                                             2833
18212 continue                                                             2833
      if(intr .eq. 0)goto 18261                                            2833
      d=sum(wr)/v0                                                         2833
      az=az+d                                                              2834
      dlx=max(dlx,v0*d**2)                                                 2834
      wr=wr-d*w                                                            2834
      f=f+d                                                                2835
18261 continue                                                             2836
      if(dlx.lt.shr)goto 18202                                             2836
      if(nlp .le. maxit)goto 18281                                         2836
      jerr=-ilm                                                            2836
      return                                                               2836
18281 continue                                                             2837
      goto 18201                                                           2838
18202 continue                                                             2838
      goto 18091                                                           2839
18092 continue                                                             2839
      if(nin.gt.nx)goto 18072                                              2840
      w=q*exp(sign(min(abs(f),fmax),f))                                    2840
      v0=sum(w)                                                            2840
      wr=t-w                                                               2841
      if(v0*(az-az0)**2 .ge. shr)goto 18301                                2841
      ix=0                                                                 2842
18310 do 18311 j=1,nin                                                     2842
      k=m(j)                                                               2843
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18311                            2843
      ix=1                                                                 2843
      goto 18312                                                           2844
18311 continue                                                             2845
18312 continue                                                             2845
      if(ix .ne. 0)goto 18331                                              2846
18340 do 18341 k=1,ni                                                      2846
      if(ixx(k).eq.1)goto 18341                                            2846
      if(ju(k).eq.0)goto 18341                                             2847
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2848
      if(ga(k) .le. al1*vp(k))goto 18361                                   2848
      ixx(k)=1                                                             2848
      ix=1                                                                 2848
18361 continue                                                             2849
18341 continue                                                             2850
18342 continue                                                             2850
      if(ix.eq.1) go to 10880                                              2851
      goto 18072                                                           2852
18331 continue                                                             2853
18301 continue                                                             2854
      goto 18071                                                           2855
18072 continue                                                             2855
      if(nin .le. nx)goto 18381                                            2855
      jerr=-10000-ilm                                                      2855
      goto 17992                                                           2855
18381 continue                                                             2856
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2856
      kin(ilm)=nin                                                         2857
      a0(ilm)=az                                                           2857
      alm(ilm)=al                                                          2857
      lmu=ilm                                                              2858
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2859
      if(ilm.lt.mnl)goto 17991                                             2859
      if(flmin.ge.1.0)goto 17991                                           2860
      me=0                                                                 2860
18390 do 18391 j=1,nin                                                     2860
      if(ca(j,ilm).ne.0.0) me=me+1                                         2860
18391 continue                                                             2860
18392 continue                                                             2860
      if(me.gt.ne)goto 17992                                               2861
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17992              2862
      if(dev(ilm).gt.devmax)goto 17992                                     2863
17991 continue                                                             2864
17992 continue                                                             2864
      g=f                                                                  2865
12180 continue                                                             2865
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2866
      return                                                               2867
      end                                                                  2868
      function nonzero(n,v)                                                2869
      real v(n)                                                            2870
      nonzero=0                                                            2870
18400 do 18401 i=1,n                                                       2870
      if(v(i) .eq. 0.0)goto 18421                                          2870
      nonzero=1                                                            2870
      return                                                               2870
18421 continue                                                             2870
18401 continue                                                             2871
18402 continue                                                             2871
      return                                                               2872
      end                                                                  2873
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2874
      real a(nx,lmu),b(ni,lmu)                                             2874
      integer ia(nx),nin(lmu)                                              2875
18430 do 18431 lam=1,lmu                                                   2875
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2875
18431 continue                                                             2876
18432 continue                                                             2876
      return                                                               2877
      end                                                                  2878
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2879
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2879
      integer ia(nx),nin(lmu)                                              2880
18440 do 18441 lam=1,lmu                                                   2880
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2880
18441 continue                                                             2881
18442 continue                                                             2881
      return                                                               2882
      end                                                                  2883
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2884
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2885
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 18461                                     2888
      jerr=8888                                                            2888
      return                                                               2888
18461 continue                                                             2889
      allocate(w(1:no),stat=jerr)                                          2889
      if(jerr.ne.0) return                                                 2890
      w=max(0.0,q)                                                         2890
      sw=sum(w)                                                            2890
      if(sw .gt. 0.0)goto 18481                                            2890
      jerr=9999                                                            2890
      go to 12180                                                          2890
18481 continue                                                             2891
      yb=dot_product(w,y)/sw                                               2891
      fmax=log(huge(y(1))*0.1)                                             2892
18490 do 18491 lam=1,nlam                                                  2892
      s=0.0                                                                2893
18500 do 18501 i=1,no                                                      2893
      if(w(i).le.0.0)goto 18501                                            2894
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2895
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2896
18501 continue                                                             2897
18502 continue                                                             2897
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2898
18491 continue                                                             2899
18492 continue                                                             2899
12180 continue                                                             2899
      deallocate(w)                                                        2900
      return                                                               2901
      end                                                                  2902
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   2904 
     *,flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)               2905
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2906
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2907
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 18521                                    2911
      jerr=10000                                                           2911
      return                                                               2911
18521 continue                                                             2912
      if(minval(y) .ge. 0.0)goto 18541                                     2912
      jerr=8888                                                            2912
      return                                                               2912
18541 continue                                                             2913
      allocate(ww(1:no),stat=jerr)                                         2914
      allocate(ju(1:ni),stat=ierr)                                         2914
      jerr=jerr+ierr                                                       2915
      allocate(vq(1:ni),stat=ierr)                                         2915
      jerr=jerr+ierr                                                       2916
      allocate(xm(1:ni),stat=ierr)                                         2916
      jerr=jerr+ierr                                                       2917
      allocate(xs(1:ni),stat=ierr)                                         2917
      jerr=jerr+ierr                                                       2918
      if(jerr.ne.0) return                                                 2919
      call spchkvars(no,ni,x,ix,ju)                                        2920
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2921
      if(maxval(ju) .gt. 0)goto 18561                                      2921
      jerr=7777                                                            2921
      go to 12180                                                          2921
18561 continue                                                             2922
      vq=max(0.0,vp)                                                       2922
      vq=vq*ni/sum(vq)                                                     2923
      ww=max(0.0,w)                                                        2923
      sw=sum(ww)                                                           2923
      if(sw .gt. 0.0)goto 18581                                            2923
      jerr=9999                                                            2923
      go to 12180                                                          2923
18581 continue                                                             2924
      ww=ww/sw                                                             2925
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2926
      if(isd .le. 0)goto 18601                                             2926
18610 do 18611 j=1,ni                                                      2926
      cl(:,j)=cl(:,j)*xs(j)                                                2926
18611 continue                                                             2926
18612 continue                                                             2926
18601 continue                                                             2927
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   2929 
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12180                                            2929
      dev0=2.0*sw*dev0                                                     2930
18620 do 18621 k=1,lmu                                                     2930
      nk=nin(k)                                                            2931
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2932
      if(intr .ne. 0)goto 18641                                            2932
      a0(k)=0.0                                                            2932
      goto 18651                                                           2933
18641 continue                                                             2933
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2933
18651 continue                                                             2934
18631 continue                                                             2934
18621 continue                                                             2935
18622 continue                                                             2935
12180 continue                                                             2935
      deallocate(ww,ju,vq,xm,xs)                                           2936
      return                                                               2937
      end                                                                  2938
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   2940 
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2941
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2942
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2943
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2947
      sml=sml*10.0                                                         2948
      allocate(a(1:ni),stat=jerr)                                          2949
      allocate(as(1:ni),stat=ierr)                                         2949
      jerr=jerr+ierr                                                       2950
      allocate(t(1:no),stat=ierr)                                          2950
      jerr=jerr+ierr                                                       2951
      allocate(mm(1:ni),stat=ierr)                                         2951
      jerr=jerr+ierr                                                       2952
      allocate(ga(1:ni),stat=ierr)                                         2952
      jerr=jerr+ierr                                                       2953
      allocate(ixx(1:ni),stat=ierr)                                        2953
      jerr=jerr+ierr                                                       2954
      allocate(wr(1:no),stat=ierr)                                         2954
      jerr=jerr+ierr                                                       2955
      allocate(v(1:ni),stat=ierr)                                          2955
      jerr=jerr+ierr                                                       2956
      allocate(xm(1:ni),stat=ierr)                                         2956
      jerr=jerr+ierr                                                       2957
      allocate(w(1:no),stat=ierr)                                          2957
      jerr=jerr+ierr                                                       2958
      allocate(qy(1:no),stat=ierr)                                         2958
      jerr=jerr+ierr                                                       2959
      if(jerr.ne.0) return                                                 2960
      bta=parm                                                             2960
      omb=1.0-bta                                                          2960
      fmax=log(huge(bta)*0.1)                                              2961
      qy=q*y                                                               2961
      yb=sum(qy)                                                           2962
      if(nonzero(no,g) .ne. 0)goto 18671                                   2962
      t=0.0                                                                2963
      if(intr .eq. 0)goto 18691                                            2963
      w=q*yb                                                               2963
      az=log(yb)                                                           2963
      uu=az                                                                2964
      xm=yb*xb                                                             2964
      dv0=yb*(az-1.0)                                                      2965
      goto 18701                                                           2966
18691 continue                                                             2966
      w=q                                                                  2966
      xm=0.0                                                               2966
      uu=0.0                                                               2966
      az=uu                                                                2966
      dv0=-1.0                                                             2966
18701 continue                                                             2967
18681 continue                                                             2967
      goto 18711                                                           2968
18671 continue                                                             2968
      w=q*exp(sign(min(abs(g),fmax),g))                                    2968
      ww=sum(w)                                                            2968
      t=g                                                                  2969
      if(intr .eq. 0)goto 18731                                            2969
      eaz=yb/ww                                                            2970
      w=eaz*w                                                              2970
      az=log(eaz)                                                          2970
      uu=az                                                                2970
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2971
      goto 18741                                                           2972
18731 continue                                                             2972
      uu=0.0                                                               2972
      az=uu                                                                2972
      dv0=dot_product(qy,g)-ww                                             2972
18741 continue                                                             2973
18721 continue                                                             2973
18750 do 18751 j=1,ni                                                      2973
      if(ju(j).eq.0)goto 18751                                             2973
      jb=ix(j)                                                             2973
      je=ix(j+1)-1                                                         2974
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2975
18751 continue                                                             2976
18752 continue                                                             2976
18711 continue                                                             2977
18661 continue                                                             2977
      tt=yb*uu                                                             2977
      ww=1.0                                                               2977
      if(intr.ne.0) ww=yb                                                  2977
      wr=qy-q*(yb*(1.0-uu))                                                2977
      a=0.0                                                                2977
      as=0.0                                                               2978
      dvr=-yb                                                              2979
18760 do 18761 i=1,no                                                      2979
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2979
18761 continue                                                             2979
18762 continue                                                             2979
      dvr=dvr-dv0                                                          2979
      dev0=dvr                                                             2980
      if(flmin .ge. 1.0)goto 18781                                         2980
      eqs=max(eps,flmin)                                                   2980
      alf=eqs**(1.0/(nlam-1))                                              2980
18781 continue                                                             2981
      m=0                                                                  2981
      mm=0                                                                 2981
      nlp=0                                                                2981
      nin=nlp                                                              2981
      mnl=min(mnlam,nlam)                                                  2981
      shr=shri*dev0                                                        2981
      al=0.0                                                               2981
      ixx=0                                                                2982
18790 do 18791 j=1,ni                                                      2982
      if(ju(j).eq.0)goto 18791                                             2983
      jb=ix(j)                                                             2983
      je=ix(j+1)-1                                                         2984
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2986 
     *)-xb(j)*tt)/xs(j)
18791 continue                                                             2987
18792 continue                                                             2987
18800 do 18801 ilm=1,nlam                                                  2987
      al0=al                                                               2988
      if(flmin .lt. 1.0)goto 18821                                         2988
      al=ulam(ilm)                                                         2988
      goto 18811                                                           2989
18821 if(ilm .le. 2)goto 18831                                             2989
      al=al*alf                                                            2989
      goto 18811                                                           2990
18831 if(ilm .ne. 1)goto 18841                                             2990
      al=big                                                               2990
      goto 18851                                                           2991
18841 continue                                                             2991
      al0=0.0                                                              2992
18860 do 18861 j=1,ni                                                      2992
      if(ju(j).eq.0)goto 18861                                             2992
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2992
18861 continue                                                             2993
18862 continue                                                             2993
      al0=al0/max(bta,1.0e-3)                                              2993
      al=alf*al0                                                           2994
18851 continue                                                             2995
18811 continue                                                             2995
      al2=al*omb                                                           2995
      al1=al*bta                                                           2995
      tlam=bta*(2.0*al-al0)                                                2996
18870 do 18871 k=1,ni                                                      2996
      if(ixx(k).eq.1)goto 18871                                            2996
      if(ju(k).eq.0)goto 18871                                             2997
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2998
18871 continue                                                             2999
18872 continue                                                             2999
10880 continue                                                             3000
18880 continue                                                             3000
18881 continue                                                             3000
      az0=az                                                               3001
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3002
18890 do 18891 j=1,ni                                                      3002
      if(ixx(j).eq.0)goto 18891                                            3002
      jb=ix(j)                                                             3002
      je=ix(j+1)-1                                                         3003
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3004
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   3006 
     *b(j)**2)/xs(j)**2
18891 continue                                                             3007
18892 continue                                                             3007
18900 continue                                                             3007
18901 continue                                                             3007
      nlp=nlp+1                                                            3008
      dlx=0.0                                                              3009
18910 do 18911 k=1,ni                                                      3009
      if(ixx(k).eq.0)goto 18911                                            3009
      jb=ix(k)                                                             3009
      je=ix(k+1)-1                                                         3009
      ak=a(k)                                                              3010
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3012 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3013
      if(au .gt. 0.0)goto 18931                                            3013
      a(k)=0.0                                                             3013
      goto 18941                                                           3014
18931 continue                                                             3015
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3016
18941 continue                                                             3017
18921 continue                                                             3017
      if(a(k).eq.ak)goto 18911                                             3018
      if(mm(k) .ne. 0)goto 18961                                           3018
      nin=nin+1                                                            3018
      if(nin.gt.nx)goto 18912                                              3019
      mm(k)=nin                                                            3019
      m(nin)=k                                                             3020
18961 continue                                                             3021
      d=a(k)-ak                                                            3021
      dlx=max(dlx,v(k)*d**2)                                               3021
      dv=d/xs(k)                                                           3022
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3023
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3024
      uu=uu-dv*xb(k)                                                       3024
      tt=tt-dv*xm(k)                                                       3025
18911 continue                                                             3026
18912 continue                                                             3026
      if(nin.gt.nx)goto 18902                                              3027
      if(intr .eq. 0)goto 18981                                            3027
      d=tt/ww-uu                                                           3028
      az=az+d                                                              3028
      dlx=max(dlx,ww*d**2)                                                 3028
      uu=uu+d                                                              3029
18981 continue                                                             3030
      if(dlx.lt.shr)goto 18902                                             3030
      if(nlp .le. maxit)goto 19001                                         3030
      jerr=-ilm                                                            3030
      return                                                               3030
19001 continue                                                             3031
19010 continue                                                             3031
19011 continue                                                             3031
      nlp=nlp+1                                                            3031
      dlx=0.0                                                              3032
19020 do 19021 l=1,nin                                                     3032
      k=m(l)                                                               3033
      jb=ix(k)                                                             3033
      je=ix(k+1)-1                                                         3033
      ak=a(k)                                                              3034
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3036 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3037
      if(au .gt. 0.0)goto 19041                                            3037
      a(k)=0.0                                                             3037
      goto 19051                                                           3038
19041 continue                                                             3039
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3040
19051 continue                                                             3041
19031 continue                                                             3041
      if(a(k).eq.ak)goto 19021                                             3041
      d=a(k)-ak                                                            3041
      dlx=max(dlx,v(k)*d**2)                                               3042
      dv=d/xs(k)                                                           3042
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3043
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3044
      uu=uu-dv*xb(k)                                                       3044
      tt=tt-dv*xm(k)                                                       3045
19021 continue                                                             3046
19022 continue                                                             3046
      if(intr .eq. 0)goto 19071                                            3046
      d=tt/ww-uu                                                           3046
      az=az+d                                                              3047
      dlx=max(dlx,ww*d**2)                                                 3047
      uu=uu+d                                                              3048
19071 continue                                                             3049
      if(dlx.lt.shr)goto 19012                                             3049
      if(nlp .le. maxit)goto 19091                                         3049
      jerr=-ilm                                                            3049
      return                                                               3049
19091 continue                                                             3050
      goto 19011                                                           3051
19012 continue                                                             3051
      goto 18901                                                           3052
18902 continue                                                             3052
      if(nin.gt.nx)goto 18882                                              3053
      euu=exp(sign(min(abs(uu),fmax),uu))                                  3054
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                3054
      ww=sum(w)                                                            3055
      wr=qy-w*(1.0-uu)                                                     3055
      tt=sum(wr)                                                           3056
      if(ww*(az-az0)**2 .ge. shr)goto 19111                                3056
      kx=0                                                                 3057
19120 do 19121 j=1,nin                                                     3057
      k=m(j)                                                               3058
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19121                            3058
      kx=1                                                                 3058
      goto 19122                                                           3059
19121 continue                                                             3060
19122 continue                                                             3060
      if(kx .ne. 0)goto 19141                                              3061
19150 do 19151 j=1,ni                                                      3061
      if(ixx(j).eq.1)goto 19151                                            3061
      if(ju(j).eq.0)goto 19151                                             3062
      jb=ix(j)                                                             3062
      je=ix(j+1)-1                                                         3063
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3064
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3066 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19171                                   3066
      ixx(j)=1                                                             3066
      kx=1                                                                 3066
19171 continue                                                             3067
19151 continue                                                             3068
19152 continue                                                             3068
      if(kx.eq.1) go to 10880                                              3069
      goto 18882                                                           3070
19141 continue                                                             3071
19111 continue                                                             3072
      goto 18881                                                           3073
18882 continue                                                             3073
      if(nin .le. nx)goto 19191                                            3073
      jerr=-10000-ilm                                                      3073
      goto 18802                                                           3073
19191 continue                                                             3074
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3074
      kin(ilm)=nin                                                         3075
      a0(ilm)=az                                                           3075
      alm(ilm)=al                                                          3075
      lmu=ilm                                                              3076
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        3077
      if(ilm.lt.mnl)goto 18801                                             3077
      if(flmin.ge.1.0)goto 18801                                           3078
      me=0                                                                 3078
19200 do 19201 j=1,nin                                                     3078
      if(ca(j,ilm).ne.0.0) me=me+1                                         3078
19201 continue                                                             3078
19202 continue                                                             3078
      if(me.gt.ne)goto 18802                                               3079
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18802              3080
      if(dev(ilm).gt.devmax)goto 18802                                     3081
18801 continue                                                             3082
18802 continue                                                             3082
      g=t+uu                                                               3083
12180 continue                                                             3083
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            3084
      return                                                               3085
      end                                                                  3086
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       3087
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           3088
      integer ix(*),jx(*)                                                  3089
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 19221                                     3092
      jerr=8888                                                            3092
      return                                                               3092
19221 continue                                                             3093
      allocate(w(1:no),stat=jerr)                                          3094
      allocate(f(1:no),stat=ierr)                                          3094
      jerr=jerr+ierr                                                       3095
      if(jerr.ne.0) return                                                 3096
      w=max(0.0,q)                                                         3096
      sw=sum(w)                                                            3096
      if(sw .gt. 0.0)goto 19241                                            3096
      jerr=9999                                                            3096
      go to 12180                                                          3096
19241 continue                                                             3097
      yb=dot_product(w,y)/sw                                               3097
      fmax=log(huge(y(1))*0.1)                                             3098
19250 do 19251 lam=1,nlam                                                  3098
      f=a0(lam)                                                            3099
19260 do 19261 j=1,ni                                                      3099
      if(a(j,lam).eq.0.0)goto 19261                                        3099
      jb=ix(j)                                                             3099
      je=ix(j+1)-1                                                         3100
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          3101
19261 continue                                                             3102
19262 continue                                                             3102
      f=f+g                                                                3103
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3104
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3105
19251 continue                                                             3106
19252 continue                                                             3106
12180 continue                                                             3106
      deallocate(w,f)                                                      3107
      return                                                               3108
      end                                                                  3109
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   3110 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          3111
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 3112
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 19281                                     3115
      jerr=8888                                                            3115
      return                                                               3115
19281 continue                                                             3116
      allocate(w(1:no),stat=jerr)                                          3117
      allocate(f(1:no),stat=ierr)                                          3117
      jerr=jerr+ierr                                                       3118
      if(jerr.ne.0) return                                                 3119
      w=max(0.0,q)                                                         3119
      sw=sum(w)                                                            3119
      if(sw .gt. 0.0)goto 19301                                            3119
      jerr=9999                                                            3119
      go to 12180                                                          3119
19301 continue                                                             3120
      yb=dot_product(w,y)/sw                                               3120
      fmax=log(huge(y(1))*0.1)                                             3121
19310 do 19311 lam=1,nlam                                                  3121
      f=a0(lam)                                                            3122
19320 do 19321 k=1,nin(lam)                                                3122
      j=ia(k)                                                              3122
      jb=ix(j)                                                             3122
      je=ix(j+1)-1                                                         3123
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         3124
19321 continue                                                             3125
19322 continue                                                             3125
      f=f+g                                                                3126
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3127
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3128
19311 continue                                                             3129
19312 continue                                                             3129
12180 continue                                                             3129
      deallocate(w,f)                                                      3130
      return                                                               3131
      end                                                                  3132
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3135 
     *in,ulam,thr,isd,jsd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      real x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)                   3136
      real ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,ni)             3137
      integer jd(*),ia(nx),nin(nlam)                                       3138
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 19341                                    3141
      jerr=10000                                                           3141
      return                                                               3141
19341 continue                                                             3142
      allocate(vq(1:ni),stat=jerr)                                         3142
      if(jerr.ne.0) return                                                 3143
      vq=max(0.0,vp)                                                       3143
      vq=vq*ni/sum(vq)                                                     3144
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   3146 
     *,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3147
      return                                                               3148
      end                                                                  3149
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3151 
     *in,ulam,thr,  isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      real vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)              3152
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3153
      integer jd(*),ia(nx),nin(nlam)                                       3154
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      real, dimension (:,:,:), allocatable :: clt                               
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      allocate(xm(1:ni),stat=ierr)                                         3160
      jerr=jerr+ierr                                                       3161
      allocate(xs(1:ni),stat=ierr)                                         3161
      jerr=jerr+ierr                                                       3162
      allocate(ym(1:nr),stat=ierr)                                         3162
      jerr=jerr+ierr                                                       3163
      allocate(ys(1:nr),stat=ierr)                                         3163
      jerr=jerr+ierr                                                       3164
      allocate(ju(1:ni),stat=ierr)                                         3164
      jerr=jerr+ierr                                                       3165
      allocate(xv(1:ni),stat=ierr)                                         3165
      jerr=jerr+ierr                                                       3166
      if(jerr.ne.0) return                                                 3167
      call chkvars(no,ni,x,ju)                                             3168
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3169
      if(maxval(ju) .gt. 0)goto 19361                                      3169
      jerr=7777                                                            3169
      return                                                               3169
19361 continue                                                             3170
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y   3171 
     *s0,jerr)
      if(jerr.ne.0) return                                                 3172
19370 do 19371 j=1,ni                                                      3172
19380 do 19381 k=1,nr                                                      3172
19390 do 19391 i=1,2                                                       3172
      clt(i,k,j)=cl(i,j)                                                   3172
19391 continue                                                             3172
19392 continue                                                             3172
19381 continue                                                             3172
19382 continue                                                             3172
19371 continue                                                             3173
19372 continue                                                             3173
      if(isd .le. 0)goto 19411                                             3173
19420 do 19421 j=1,ni                                                      3173
19430 do 19431 k=1,nr                                                      3173
19440 do 19441 i=1,2                                                       3173
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3173
19441 continue                                                             3173
19442 continue                                                             3173
19431 continue                                                             3173
19432 continue                                                             3173
19421 continue                                                             3173
19422 continue                                                             3173
19411 continue                                                             3174
      if(jsd .le. 0)goto 19461                                             3174
19470 do 19471 j=1,ni                                                      3174
19480 do 19481 k=1,nr                                                      3174
19490 do 19491 i=1,2                                                       3174
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3174
19491 continue                                                             3174
19492 continue                                                             3174
19481 continue                                                             3174
19482 continue                                                             3174
19471 continue                                                             3174
19472 continue                                                             3174
19461 continue                                                             3175
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   3177 
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3178
19500 do 19501 k=1,lmu                                                     3178
      nk=nin(k)                                                            3179
19510 do 19511 j=1,nr                                                      3180
19520 do 19521 l=1,nk                                                      3180
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3180
19521 continue                                                             3181
19522 continue                                                             3181
      if(intr .ne. 0)goto 19541                                            3181
      a0(j,k)=0.0                                                          3181
      goto 19551                                                           3182
19541 continue                                                             3182
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3182
19551 continue                                                             3183
19531 continue                                                             3183
19511 continue                                                             3184
19512 continue                                                             3184
19501 continue                                                             3185
19502 continue                                                             3185
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3186
      return                                                               3187
      end                                                                  3188
      subroutine multstandard1  (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym   3190 
     *,ys,xv,ys0,jerr)
      real x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)      3191
      integer ju(ni)                                                       3192
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                          3195
      if(jerr.ne.0) return                                                 3196
      w=w/sum(w)                                                           3196
      v=sqrt(w)                                                            3197
      if(intr .ne. 0)goto 19571                                            3198
19580 do 19581 j=1,ni                                                      3198
      if(ju(j).eq.0)goto 19581                                             3198
      xm(j)=0.0                                                            3198
      x(:,j)=v*x(:,j)                                                      3199
      z=dot_product(x(:,j),x(:,j))                                         3200
      if(isd .le. 0)goto 19601                                             3200
      xbq=dot_product(v,x(:,j))**2                                         3200
      vc=z-xbq                                                             3201
      xs(j)=sqrt(vc)                                                       3201
      x(:,j)=x(:,j)/xs(j)                                                  3201
      xv(j)=1.0+xbq/vc                                                     3202
      goto 19611                                                           3203
19601 continue                                                             3203
      xs(j)=1.0                                                            3203
      xv(j)=z                                                              3203
19611 continue                                                             3204
19591 continue                                                             3204
19581 continue                                                             3205
19582 continue                                                             3205
      ys0=0.0                                                              3206
19620 do 19621 j=1,nr                                                      3206
      ym(j)=0.0                                                            3206
      y(:,j)=v*y(:,j)                                                      3207
      z=dot_product(y(:,j),y(:,j))                                         3208
      if(jsd .le. 0)goto 19641                                             3208
      u=z-dot_product(v,y(:,j))**2                                         3208
      ys0=ys0+z/u                                                          3209
      ys(j)=sqrt(u)                                                        3209
      y(:,j)=y(:,j)/ys(j)                                                  3210
      goto 19651                                                           3211
19641 continue                                                             3211
      ys(j)=1.0                                                            3211
      ys0=ys0+z                                                            3211
19651 continue                                                             3212
19631 continue                                                             3212
19621 continue                                                             3213
19622 continue                                                             3213
      go to 10700                                                          3214
19571 continue                                                             3215
19660 do 19661 j=1,ni                                                      3215
      if(ju(j).eq.0)goto 19661                                             3216
      xm(j)=dot_product(w,x(:,j))                                          3216
      x(:,j)=v*(x(:,j)-xm(j))                                              3217
      xv(j)=dot_product(x(:,j),x(:,j))                                     3217
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3218
19661 continue                                                             3219
19662 continue                                                             3219
      if(isd .ne. 0)goto 19681                                             3219
      xs=1.0                                                               3219
      goto 19691                                                           3220
19681 continue                                                             3220
19700 do 19701 j=1,ni                                                      3220
      if(ju(j).eq.0)goto 19701                                             3220
      x(:,j)=x(:,j)/xs(j)                                                  3220
19701 continue                                                             3221
19702 continue                                                             3221
      xv=1.0                                                               3222
19691 continue                                                             3223
19671 continue                                                             3223
      ys0=0.0                                                              3224
19710 do 19711 j=1,nr                                                      3225
      ym(j)=dot_product(w,y(:,j))                                          3225
      y(:,j)=v*(y(:,j)-ym(j))                                              3226
      z=dot_product(y(:,j),y(:,j))                                         3227
      if(jsd .le. 0)goto 19731                                             3227
      ys(j)=sqrt(z)                                                        3227
      y(:,j)=y(:,j)/ys(j)                                                  3227
      goto 19741                                                           3228
19731 continue                                                             3228
      ys0=ys0+z                                                            3228
19741 continue                                                             3229
19721 continue                                                             3229
19711 continue                                                             3230
19712 continue                                                             3230
      if(jsd .ne. 0)goto 19761                                             3230
      ys=1.0                                                               3230
      goto 19771                                                           3230
19761 continue                                                             3230
      ys0=nr                                                               3230
19771 continue                                                             3231
19751 continue                                                             3231
10700 continue                                                             3231
      deallocate(v)                                                        3232
      return                                                               3233
      end                                                                  3234
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   3236 
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam)              3237
      real rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)                        3238
      integer ju(ni),ia(nx),kin(nlam)                                      3239
      real, dimension (:), allocatable :: g,gk,del,gj                           
      integer, dimension (:), allocatable :: mm,ix,isc                          
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3246
      allocate(gj(1:nr),stat=ierr)                                         3246
      jerr=jerr+ierr                                                       3247
      allocate(gk(1:nr),stat=ierr)                                         3247
      jerr=jerr+ierr                                                       3248
      allocate(del(1:nr),stat=ierr)                                        3248
      jerr=jerr+ierr                                                       3249
      allocate(mm(1:ni),stat=ierr)                                         3249
      jerr=jerr+ierr                                                       3250
      allocate(g(1:ni),stat=ierr)                                          3250
      jerr=jerr+ierr                                                       3251
      allocate(ix(1:ni),stat=ierr)                                         3251
      jerr=jerr+ierr                                                       3252
      allocate(isc(1:nr),stat=ierr)                                        3252
      jerr=jerr+ierr                                                       3253
      if(jerr.ne.0) return                                                 3254
      bta=beta                                                             3254
      omb=1.0-bta                                                          3254
      ix=0                                                                 3254
      thr=thri*ys0/nr                                                      3255
      if(flmin .ge. 1.0)goto 19791                                         3255
      eqs=max(eps,flmin)                                                   3255
      alf=eqs**(1.0/(nlam-1))                                              3255
19791 continue                                                             3256
      rsq=ys0                                                              3256
      a=0.0                                                                3256
      mm=0                                                                 3256
      nlp=0                                                                3256
      nin=nlp                                                              3256
      iz=0                                                                 3256
      mnl=min(mnlam,nlam)                                                  3256
      alm=0.0                                                              3257
19800 do 19801 j=1,ni                                                      3257
      if(ju(j).eq.0)goto 19801                                             3257
      g(j)=0.0                                                             3258
19810 do 19811 k=1,nr                                                      3258
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              3258
19811 continue                                                             3259
19812 continue                                                             3259
      g(j)=sqrt(g(j))                                                      3260
19801 continue                                                             3261
19802 continue                                                             3261
19820 do 19821 m=1,nlam                                                    3262
      if(flmin .lt. 1.0)goto 19841                                         3262
      alm=ulam(m)                                                          3262
      goto 19831                                                           3263
19841 if(m .le. 2)goto 19851                                               3263
      alm=alm*alf                                                          3263
      goto 19831                                                           3264
19851 if(m .ne. 1)goto 19861                                               3264
      alm=big                                                              3264
      goto 19871                                                           3265
19861 continue                                                             3265
      alm0=0.0                                                             3266
19880 do 19881 j=1,ni                                                      3266
      if(ju(j).eq.0)goto 19881                                             3267
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3268
19881 continue                                                             3269
19882 continue                                                             3269
      alm0=alm0/max(bta,1.0e-3)                                            3269
      alm=alf*alm0                                                         3270
19871 continue                                                             3271
19831 continue                                                             3271
      dem=alm*omb                                                          3271
      ab=alm*bta                                                           3271
      rsq0=rsq                                                             3271
      jz=1                                                                 3272
      tlam=bta*(2.0*alm-alm0)                                              3273
19890 do 19891 k=1,ni                                                      3273
      if(ix(k).eq.1)goto 19891                                             3273
      if(ju(k).eq.0)goto 19891                                             3274
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       3275
19891 continue                                                             3276
19892 continue                                                             3276
19900 continue                                                             3276
19901 continue                                                             3276
      if(iz*jz.ne.0) go to 10360                                           3277
10880 continue                                                             3277
      nlp=nlp+1                                                            3277
      dlx=0.0                                                              3278
19910 do 19911 k=1,ni                                                      3278
      if(ix(k).eq.0)goto 19911                                             3278
      gkn=0.0                                                              3279
19920 do 19921 j=1,nr                                                      3279
      gj(j)=dot_product(y(:,j),x(:,k))                                     3280
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3280
      gkn=gkn+gk(j)**2                                                     3282
19921 continue                                                             3282
19922 continue                                                             3282
      gkn=sqrt(gkn)                                                        3282
      u=1.0-ab*vp(k)/gkn                                                   3282
      del=a(:,k)                                                           3283
      if(u .gt. 0.0)goto 19941                                             3283
      a(:,k)=0.0                                                           3283
      goto 19951                                                           3284
19941 continue                                                             3284
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3285
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3287 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3288
19951 continue                                                             3289
19931 continue                                                             3289
      del=a(:,k)-del                                                       3289
      if(maxval(abs(del)).le.0.0)goto 19911                                3290
19960 do 19961 j=1,nr                                                      3290
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3291
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3291
      dlx=max(dlx,xv(k)*del(j)**2)                                         3292
19961 continue                                                             3293
19962 continue                                                             3293
      if(mm(k) .ne. 0)goto 19981                                           3293
      nin=nin+1                                                            3293
      if(nin.gt.nx)goto 19912                                              3294
      mm(k)=nin                                                            3294
      ia(nin)=k                                                            3295
19981 continue                                                             3296
19911 continue                                                             3297
19912 continue                                                             3297
      if(nin.gt.nx)goto 19902                                              3298
      if(dlx .ge. thr)goto 20001                                           3298
      ixx=0                                                                3299
20010 do 20011 k=1,ni                                                      3299
      if(ix(k).eq.1)goto 20011                                             3299
      if(ju(k).eq.0)goto 20011                                             3299
      g(k)=0.0                                                             3300
20020 do 20021 j=1,nr                                                      3300
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3300
20021 continue                                                             3301
20022 continue                                                             3301
      g(k)=sqrt(g(k))                                                      3302
      if(g(k) .le. ab*vp(k))goto 20041                                     3302
      ix(k)=1                                                              3302
      ixx=1                                                                3302
20041 continue                                                             3303
20011 continue                                                             3304
20012 continue                                                             3304
      if(ixx.eq.1) go to 10880                                             3305
      goto 19902                                                           3306
20001 continue                                                             3307
      if(nlp .le. maxit)goto 20061                                         3307
      jerr=-m                                                              3307
      return                                                               3307
20061 continue                                                             3308
10360 continue                                                             3308
      iz=1                                                                 3309
20070 continue                                                             3309
20071 continue                                                             3309
      nlp=nlp+1                                                            3309
      dlx=0.0                                                              3310
20080 do 20081 l=1,nin                                                     3310
      k=ia(l)                                                              3310
      gkn=0.0                                                              3311
20090 do 20091 j=1,nr                                                      3311
      gj(j)=dot_product(y(:,j),x(:,k))                                     3312
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3312
      gkn=gkn+gk(j)**2                                                     3314
20091 continue                                                             3314
20092 continue                                                             3314
      gkn=sqrt(gkn)                                                        3314
      u=1.0-ab*vp(k)/gkn                                                   3314
      del=a(:,k)                                                           3315
      if(u .gt. 0.0)goto 20111                                             3315
      a(:,k)=0.0                                                           3315
      goto 20121                                                           3316
20111 continue                                                             3316
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3317
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3319 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3320
20121 continue                                                             3321
20101 continue                                                             3321
      del=a(:,k)-del                                                       3321
      if(maxval(abs(del)).le.0.0)goto 20081                                3322
20130 do 20131 j=1,nr                                                      3322
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3323
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3323
      dlx=max(dlx,xv(k)*del(j)**2)                                         3324
20131 continue                                                             3325
20132 continue                                                             3325
20081 continue                                                             3326
20082 continue                                                             3326
      if(dlx.lt.thr)goto 20072                                             3326
      if(nlp .le. maxit)goto 20151                                         3326
      jerr=-m                                                              3326
      return                                                               3326
20151 continue                                                             3327
      goto 20071                                                           3328
20072 continue                                                             3328
      jz=0                                                                 3329
      goto 19901                                                           3330
19902 continue                                                             3330
      if(nin .le. nx)goto 20171                                            3330
      jerr=-10000-m                                                        3330
      goto 19822                                                           3330
20171 continue                                                             3331
      if(nin .le. 0)goto 20191                                             3331
20200 do 20201 j=1,nr                                                      3331
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3331
20201 continue                                                             3331
20202 continue                                                             3331
20191 continue                                                             3332
      kin(m)=nin                                                           3333
      rsqo(m)=1.0-rsq/ys0                                                  3333
      almo(m)=alm                                                          3333
      lmu=m                                                                3334
      if(m.lt.mnl)goto 19821                                               3334
      if(flmin.ge.1.0)goto 19821                                           3335
      me=0                                                                 3335
20210 do 20211 j=1,nin                                                     3335
      if(ao(j,1,m).ne.0.0) me=me+1                                         3335
20211 continue                                                             3335
20212 continue                                                             3335
      if(me.gt.ne)goto 19822                                               3336
      if(rsq0-rsq.lt.sml*rsq)goto 19822                                    3336
      if(rsqo(m).gt.rsqmax)goto 19822                                      3337
19821 continue                                                             3338
19822 continue                                                             3338
      deallocate(a,mm,g,ix,del,gj,gk)                                      3339
      return                                                               3340
      end                                                                  3341
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               3342
      real gk(nr),cl(2,nr),a(nr)                                           3342
      integer isc(nr)                                                      3343
      kerr=0                                                               3343
      al1p=1.0+al1/xv                                                      3343
      al2p=al2/xv                                                          3343
      isc=0                                                                3344
      gsq=gkn**2                                                           3344
      asq=dot_product(a,a)                                                 3344
      usq=0.0                                                              3345
20220 continue                                                             3345
20221 continue                                                             3345
      vmx=0.0                                                              3346
20230 do 20231 k=1,nr                                                      3346
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     3347
      if(v .le. vmx)goto 20251                                             3347
      vmx=v                                                                3347
      kn=k                                                                 3347
20251 continue                                                             3348
20231 continue                                                             3349
20232 continue                                                             3349
      if(vmx.le.0.0)goto 20222                                             3349
      if(isc(kn).ne.0)goto 20222                                           3350
      gsq=gsq-gk(kn)**2                                                    3350
      g=sqrt(gsq)/xv                                                       3351
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     3351
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     3352
      usq=usq+u**2                                                         3353
      if(usq .ne. 0.0)goto 20271                                           3353
      b=max(0.0,(g-al2p)/al1p)                                             3353
      goto 20281                                                           3354
20271 continue                                                             3354
      b0=sqrt(asq-a(kn)**2)                                                3355
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3355
      if(kerr.ne.0)goto 20222                                              3356
20281 continue                                                             3357
20261 continue                                                             3357
      asq=usq+b**2                                                         3357
      if(asq .gt. 0.0)goto 20301                                           3357
      a=0.0                                                                3357
      goto 20222                                                           3357
20301 continue                                                             3358
      a(kn)=u                                                              3358
      isc(kn)=1                                                            3358
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3359
20310 do 20311 j=1,nr                                                      3359
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3359
20311 continue                                                             3360
20312 continue                                                             3360
      goto 20221                                                           3361
20222 continue                                                             3361
      if(kerr.ne.0) jerr=kerr                                              3362
      return                                                               3363
      end                                                                  3364
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         3365
      real gk(nr),a(nr)                                                    3365
      integer isc(nr)                                                      3366
      kerr=0                                                               3366
      al1p=1.0+al1/xv                                                      3366
      al2p=al2/xv                                                          3366
      isc=0                                                                3367
      gsq=gkn**2                                                           3367
      asq=dot_product(a,a)                                                 3367
      usq=0.0                                                              3368
20320 continue                                                             3368
20321 continue                                                             3368
      vmx=0.0                                                              3369
20330 do 20331 k=1,nr                                                      3369
      v=max(a(k)-cl2,cl1-a(k))                                             3370
      if(v .le. vmx)goto 20351                                             3370
      vmx=v                                                                3370
      kn=k                                                                 3370
20351 continue                                                             3371
20331 continue                                                             3372
20332 continue                                                             3372
      if(vmx.le.0.0)goto 20322                                             3372
      if(isc(kn).ne.0)goto 20322                                           3373
      gsq=gsq-gk(kn)**2                                                    3373
      g=sqrt(gsq)/xv                                                       3374
      if(a(kn).lt.cl1) u=cl1                                               3374
      if(a(kn).gt.cl2) u=cl2                                               3375
      usq=usq+u**2                                                         3376
      if(usq .ne. 0.0)goto 20371                                           3376
      b=max(0.0,(g-al2p)/al1p)                                             3376
      goto 20381                                                           3377
20371 continue                                                             3377
      b0=sqrt(asq-a(kn)**2)                                                3378
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3378
      if(kerr.ne.0)goto 20322                                              3379
20381 continue                                                             3380
20361 continue                                                             3380
      asq=usq+b**2                                                         3380
      if(asq .gt. 0.0)goto 20401                                           3380
      a=0.0                                                                3380
      goto 20322                                                           3380
20401 continue                                                             3381
      a(kn)=u                                                              3381
      isc(kn)=1                                                            3381
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3382
20410 do 20411 j=1,nr                                                      3382
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3382
20411 continue                                                             3383
20412 continue                                                             3383
      goto 20321                                                           3384
20322 continue                                                             3384
      if(kerr.ne.0) jerr=kerr                                              3385
      return                                                               3386
      end                                                                  3387
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              3388
      data thr,mxit /1.0e-10,100/                                          3389
      b=b0                                                                 3389
      zsq=b**2+usq                                                         3389
      if(zsq .gt. 0.0)goto 20431                                           3389
      bnorm=0.0                                                            3389
      return                                                               3389
20431 continue                                                             3390
      z=sqrt(zsq)                                                          3390
      f=b*(al1p+al2p/z)-g                                                  3390
      jerr=0                                                               3391
20440 do 20441 it=1,mxit                                                   3391
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        3392
      zsq=b**2+usq                                                         3392
      if(zsq .gt. 0.0)goto 20461                                           3392
      bnorm=0.0                                                            3392
      return                                                               3392
20461 continue                                                             3393
      z=sqrt(zsq)                                                          3393
      f=b*(al1p+al2p/z)-g                                                  3394
      if(abs(f).le.thr)goto 20442                                          3394
      if(b .gt. 0.0)goto 20481                                             3394
      b=0.0                                                                3394
      goto 20442                                                           3394
20481 continue                                                             3395
20441 continue                                                             3396
20442 continue                                                             3396
      bnorm=b                                                              3396
      if(it.ge.mxit) jerr=90000                                            3397
      return                                                               3398
      entry chg_bnorm(arg,irg)                                             3398
      thr=arg                                                              3398
      mxit=irg                                                             3398
      return                                                               3399
      entry get_bnorm(arg,irg)                                             3399
      arg=thr                                                              3399
      irg=mxit                                                             3399
      return                                                               3400
      end                                                                  3401
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3402
      real a(nx,nr,lmu),b(ni,nr,lmu)                                       3402
      integer ia(nx),nin(lmu)                                              3403
20490 do 20491 lam=1,lmu                                                   3403
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3403
20491 continue                                                             3404
20492 continue                                                             3404
      return                                                               3405
      end                                                                  3406
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3407
      real ca(nx,nr),a(ni,nr)                                              3407
      integer ia(nx)                                                       3408
      a=0.0                                                                3409
      if(nin .le. 0)goto 20511                                             3409
20520 do 20521 j=1,nr                                                      3409
      a(ia(1:nin),j)=ca(1:nin,j)                                           3409
20521 continue                                                             3409
20522 continue                                                             3409
20511 continue                                                             3410
      return                                                               3411
      end                                                                  3412
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3413
      real a0(nr),ca(nx,nr),x(n,*),f(nr,n)                                 3413
      integer ia(nx)                                                       3414
20530 do 20531 i=1,n                                                       3414
      f(:,i)=a0                                                            3414
20531 continue                                                             3414
20532 continue                                                             3414
      if(nin.le.0) return                                                  3415
20540 do 20541 i=1,n                                                       3415
20550 do 20551 j=1,nr                                                      3415
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3415
20551 continue                                                             3415
20552 continue                                                             3415
20541 continue                                                             3416
20542 continue                                                             3416
      return                                                               3417
      end                                                                  3418
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   3421 
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      real x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)                  3422
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3423
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3424
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 20571                                    3427
      jerr=10000                                                           3427
      return                                                               3427
20571 continue                                                             3428
      allocate(vq(1:ni),stat=jerr)                                         3428
      if(jerr.ne.0) return                                                 3429
      vq=max(0.0,vp)                                                       3429
      vq=vq*ni/sum(vq)                                                     3430
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   3432 
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                       3433
      return                                                               3434
      end                                                                  3435
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   3437 
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      real x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)                  3438
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3439
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3440
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      real, dimension (:,:,:), allocatable :: clt                               
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      allocate(xm(1:ni),stat=ierr)                                         3446
      jerr=jerr+ierr                                                       3447
      allocate(xs(1:ni),stat=ierr)                                         3447
      jerr=jerr+ierr                                                       3448
      allocate(ym(1:nr),stat=ierr)                                         3448
      jerr=jerr+ierr                                                       3449
      allocate(ys(1:nr),stat=ierr)                                         3449
      jerr=jerr+ierr                                                       3450
      allocate(ju(1:ni),stat=ierr)                                         3450
      jerr=jerr+ierr                                                       3451
      allocate(xv(1:ni),stat=ierr)                                         3451
      jerr=jerr+ierr                                                       3452
      if(jerr.ne.0) return                                                 3453
      call spchkvars(no,ni,x,ix,ju)                                        3454
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3455
      if(maxval(ju) .gt. 0)goto 20591                                      3455
      jerr=7777                                                            3455
      return                                                               3455
20591 continue                                                             3456
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,   3458 
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                                 3459
20600 do 20601 j=1,ni                                                      3459
20610 do 20611 k=1,nr                                                      3459
20620 do 20621 i=1,2                                                       3459
      clt(i,k,j)=cl(i,j)                                                   3459
20621 continue                                                             3459
20622 continue                                                             3459
20611 continue                                                             3459
20612 continue                                                             3459
20601 continue                                                             3460
20602 continue                                                             3460
      if(isd .le. 0)goto 20641                                             3460
20650 do 20651 j=1,ni                                                      3460
20660 do 20661 k=1,nr                                                      3460
20670 do 20671 i=1,2                                                       3460
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3460
20671 continue                                                             3460
20672 continue                                                             3460
20661 continue                                                             3460
20662 continue                                                             3460
20651 continue                                                             3460
20652 continue                                                             3460
20641 continue                                                             3461
      if(jsd .le. 0)goto 20691                                             3461
20700 do 20701 j=1,ni                                                      3461
20710 do 20711 k=1,nr                                                      3461
20720 do 20721 i=1,2                                                       3461
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3461
20721 continue                                                             3461
20722 continue                                                             3461
20711 continue                                                             3461
20712 continue                                                             3461
20701 continue                                                             3461
20702 continue                                                             3461
20691 continue                                                             3462
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   3464 
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3465
20730 do 20731 k=1,lmu                                                     3465
      nk=nin(k)                                                            3466
20740 do 20741 j=1,nr                                                      3467
20750 do 20751 l=1,nk                                                      3467
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3467
20751 continue                                                             3468
20752 continue                                                             3468
      if(intr .ne. 0)goto 20771                                            3468
      a0(j,k)=0.0                                                          3468
      goto 20781                                                           3469
20771 continue                                                             3469
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3469
20781 continue                                                             3470
20761 continue                                                             3470
20741 continue                                                             3471
20742 continue                                                             3471
20731 continue                                                             3472
20732 continue                                                             3472
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3473
      return                                                               3474
      end                                                                  3475
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,     3477 
     *xm,xs,ym,ys,xv,ys0,jerr)
      real x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)          3478
      integer ix(*),jx(*),ju(ni)                                           3479
      w=w/sum(w)                                                           3480
      if(intr .ne. 0)goto 20801                                            3481
20810 do 20811 j=1,ni                                                      3481
      if(ju(j).eq.0)goto 20811                                             3481
      xm(j)=0.0                                                            3481
      jb=ix(j)                                                             3481
      je=ix(j+1)-1                                                         3482
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                              3483
      if(isd .le. 0)goto 20831                                             3483
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            3483
      vc=z-xbq                                                             3484
      xs(j)=sqrt(vc)                                                       3484
      xv(j)=1.0+xbq/vc                                                     3485
      goto 20841                                                           3486
20831 continue                                                             3486
      xs(j)=1.0                                                            3486
      xv(j)=z                                                              3486
20841 continue                                                             3487
20821 continue                                                             3487
20811 continue                                                             3488
20812 continue                                                             3488
      ys0=0.0                                                              3489
20850 do 20851 j=1,nr                                                      3489
      ym(j)=0.0                                                            3489
      z=dot_product(w,y(:,j)**2)                                           3490
      if(jsd .le. 0)goto 20871                                             3490
      u=z-dot_product(w,y(:,j))**2                                         3490
      ys0=ys0+z/u                                                          3491
      ys(j)=sqrt(u)                                                        3491
      y(:,j)=y(:,j)/ys(j)                                                  3492
      goto 20881                                                           3493
20871 continue                                                             3493
      ys(j)=1.0                                                            3493
      ys0=ys0+z                                                            3493
20881 continue                                                             3494
20861 continue                                                             3494
20851 continue                                                             3495
20852 continue                                                             3495
      return                                                               3496
20801 continue                                                             3497
20890 do 20891 j=1,ni                                                      3497
      if(ju(j).eq.0)goto 20891                                             3498
      jb=ix(j)                                                             3498
      je=ix(j+1)-1                                                         3498
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3499
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3500
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3501
20891 continue                                                             3502
20892 continue                                                             3502
      if(isd .ne. 0)goto 20911                                             3502
      xs=1.0                                                               3502
      goto 20921                                                           3502
20911 continue                                                             3502
      xv=1.0                                                               3502
20921 continue                                                             3503
20901 continue                                                             3503
      ys0=0.0                                                              3504
20930 do 20931 j=1,nr                                                      3505
      ym(j)=dot_product(w,y(:,j))                                          3505
      y(:,j)=y(:,j)-ym(j)                                                  3506
      z=dot_product(w,y(:,j)**2)                                           3507
      if(jsd .le. 0)goto 20951                                             3507
      ys(j)=sqrt(z)                                                        3507
      y(:,j)=y(:,j)/ys(j)                                                  3507
      goto 20961                                                           3508
20951 continue                                                             3508
      ys0=ys0+z                                                            3508
20961 continue                                                             3509
20941 continue                                                             3509
20931 continue                                                             3510
20932 continue                                                             3510
      if(jsd .ne. 0)goto 20981                                             3510
      ys=1.0                                                               3510
      goto 20991                                                           3510
20981 continue                                                             3510
      ys0=nr                                                               3510
20991 continue                                                             3511
20971 continue                                                             3511
      return                                                               3512
      end                                                                  3513
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   3515 
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      real y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)               3516
      real ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)       3517
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3518
      real, dimension (:), allocatable :: g,gj,gk,del,o                         
      integer, dimension (:), allocatable :: mm,iy,isc                          
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3525
      allocate(mm(1:ni),stat=ierr)                                         3525
      jerr=jerr+ierr                                                       3526
      allocate(g(1:ni),stat=ierr)                                          3526
      jerr=jerr+ierr                                                       3527
      allocate(gj(1:nr),stat=ierr)                                         3527
      jerr=jerr+ierr                                                       3528
      allocate(gk(1:nr),stat=ierr)                                         3528
      jerr=jerr+ierr                                                       3529
      allocate(del(1:nr),stat=ierr)                                        3529
      jerr=jerr+ierr                                                       3530
      allocate(o(1:nr),stat=ierr)                                          3530
      jerr=jerr+ierr                                                       3531
      allocate(iy(1:ni),stat=ierr)                                         3531
      jerr=jerr+ierr                                                       3532
      allocate(isc(1:nr),stat=ierr)                                        3532
      jerr=jerr+ierr                                                       3533
      if(jerr.ne.0) return                                                 3534
      bta=beta                                                             3534
      omb=1.0-bta                                                          3534
      alm=0.0                                                              3534
      iy=0                                                                 3534
      thr=thri*ys0/nr                                                      3535
      if(flmin .ge. 1.0)goto 21011                                         3535
      eqs=max(eps,flmin)                                                   3535
      alf=eqs**(1.0/(nlam-1))                                              3535
21011 continue                                                             3536
      rsq=ys0                                                              3536
      a=0.0                                                                3536
      mm=0                                                                 3536
      o=0.0                                                                3536
      nlp=0                                                                3536
      nin=nlp                                                              3536
      iz=0                                                                 3536
      mnl=min(mnlam,nlam)                                                  3537
21020 do 21021 j=1,ni                                                      3537
      if(ju(j).eq.0)goto 21021                                             3537
      jb=ix(j)                                                             3537
      je=ix(j+1)-1                                                         3537
      g(j)=0.0                                                             3538
21030 do 21031 k=1,nr                                                      3539
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3540 
     *)**2
21031 continue                                                             3541
21032 continue                                                             3541
      g(j)=sqrt(g(j))                                                      3542
21021 continue                                                             3543
21022 continue                                                             3543
21040 do 21041 m=1,nlam                                                    3543
      alm0=alm                                                             3544
      if(flmin .lt. 1.0)goto 21061                                         3544
      alm=ulam(m)                                                          3544
      goto 21051                                                           3545
21061 if(m .le. 2)goto 21071                                               3545
      alm=alm*alf                                                          3545
      goto 21051                                                           3546
21071 if(m .ne. 1)goto 21081                                               3546
      alm=big                                                              3546
      goto 21091                                                           3547
21081 continue                                                             3547
      alm0=0.0                                                             3548
21100 do 21101 j=1,ni                                                      3548
      if(ju(j).eq.0)goto 21101                                             3549
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3550
21101 continue                                                             3551
21102 continue                                                             3551
      alm0=alm0/max(bta,1.0e-3)                                            3551
      alm=alf*alm0                                                         3552
21091 continue                                                             3553
21051 continue                                                             3553
      dem=alm*omb                                                          3553
      ab=alm*bta                                                           3553
      rsq0=rsq                                                             3553
      jz=1                                                                 3554
      tlam=bta*(2.0*alm-alm0)                                              3555
21110 do 21111 k=1,ni                                                      3555
      if(iy(k).eq.1)goto 21111                                             3555
      if(ju(k).eq.0)goto 21111                                             3556
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3557
21111 continue                                                             3558
21112 continue                                                             3558
21120 continue                                                             3558
21121 continue                                                             3558
      if(iz*jz.ne.0) go to 10360                                           3559
10880 continue                                                             3559
      nlp=nlp+1                                                            3559
      dlx=0.0                                                              3560
21130 do 21131 k=1,ni                                                      3560
      if(iy(k).eq.0)goto 21131                                             3560
      jb=ix(k)                                                             3560
      je=ix(k+1)-1                                                         3560
      gkn=0.0                                                              3561
21140 do 21141 j=1,nr                                                      3562
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3563
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3563
      gkn=gkn+gk(j)**2                                                     3564
21141 continue                                                             3565
21142 continue                                                             3565
      gkn=sqrt(gkn)                                                        3565
      u=1.0-ab*vp(k)/gkn                                                   3565
      del=a(:,k)                                                           3566
      if(u .gt. 0.0)goto 21161                                             3566
      a(:,k)=0.0                                                           3566
      goto 21171                                                           3567
21161 continue                                                             3567
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3568
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3570 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3571
21171 continue                                                             3572
21151 continue                                                             3572
      del=a(:,k)-del                                                       3572
      if(maxval(abs(del)).le.0.0)goto 21131                                3573
      if(mm(k) .ne. 0)goto 21191                                           3573
      nin=nin+1                                                            3573
      if(nin.gt.nx)goto 21132                                              3574
      mm(k)=nin                                                            3574
      ia(nin)=k                                                            3575
21191 continue                                                             3576
21200 do 21201 j=1,nr                                                      3576
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3577
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3578
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3578
      dlx=max(xv(k)*del(j)**2,dlx)                                         3579
21201 continue                                                             3580
21202 continue                                                             3580
21131 continue                                                             3581
21132 continue                                                             3581
      if(nin.gt.nx)goto 21122                                              3582
      if(dlx .ge. thr)goto 21221                                           3582
      ixx=0                                                                3583
21230 do 21231 j=1,ni                                                      3583
      if(iy(j).eq.1)goto 21231                                             3583
      if(ju(j).eq.0)goto 21231                                             3584
      jb=ix(j)                                                             3584
      je=ix(j+1)-1                                                         3584
      g(j)=0.0                                                             3585
21240 do 21241 k=1,nr                                                      3585
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3587 
     *)/xs(j))**2
21241 continue                                                             3588
21242 continue                                                             3588
      g(j)=sqrt(g(j))                                                      3589
      if(g(j) .le. ab*vp(j))goto 21261                                     3589
      iy(j)=1                                                              3589
      ixx=1                                                                3589
21261 continue                                                             3590
21231 continue                                                             3591
21232 continue                                                             3591
      if(ixx.eq.1) go to 10880                                             3592
      goto 21122                                                           3593
21221 continue                                                             3594
      if(nlp .le. maxit)goto 21281                                         3594
      jerr=-m                                                              3594
      return                                                               3594
21281 continue                                                             3595
10360 continue                                                             3595
      iz=1                                                                 3596
21290 continue                                                             3596
21291 continue                                                             3596
      nlp=nlp+1                                                            3596
      dlx=0.0                                                              3597
21300 do 21301 l=1,nin                                                     3597
      k=ia(l)                                                              3597
      jb=ix(k)                                                             3597
      je=ix(k+1)-1                                                         3597
      gkn=0.0                                                              3598
21310 do 21311 j=1,nr                                                      3598
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3600 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3600
      gkn=gkn+gk(j)**2                                                     3601
21311 continue                                                             3602
21312 continue                                                             3602
      gkn=sqrt(gkn)                                                        3602
      u=1.0-ab*vp(k)/gkn                                                   3602
      del=a(:,k)                                                           3603
      if(u .gt. 0.0)goto 21331                                             3603
      a(:,k)=0.0                                                           3603
      goto 21341                                                           3604
21331 continue                                                             3604
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3605
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3607 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3608
21341 continue                                                             3609
21321 continue                                                             3609
      del=a(:,k)-del                                                       3609
      if(maxval(abs(del)).le.0.0)goto 21301                                3610
21350 do 21351 j=1,nr                                                      3610
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3611
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3612
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3612
      dlx=max(xv(k)*del(j)**2,dlx)                                         3613
21351 continue                                                             3614
21352 continue                                                             3614
21301 continue                                                             3615
21302 continue                                                             3615
      if(dlx.lt.thr)goto 21292                                             3615
      if(nlp .le. maxit)goto 21371                                         3615
      jerr=-m                                                              3615
      return                                                               3615
21371 continue                                                             3616
      goto 21291                                                           3617
21292 continue                                                             3617
      jz=0                                                                 3618
      goto 21121                                                           3619
21122 continue                                                             3619
      if(nin .le. nx)goto 21391                                            3619
      jerr=-10000-m                                                        3619
      goto 21042                                                           3619
21391 continue                                                             3620
      if(nin .le. 0)goto 21411                                             3620
21420 do 21421 j=1,nr                                                      3620
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3620
21421 continue                                                             3620
21422 continue                                                             3620
21411 continue                                                             3621
      kin(m)=nin                                                           3622
      rsqo(m)=1.0-rsq/ys0                                                  3622
      almo(m)=alm                                                          3622
      lmu=m                                                                3623
      if(m.lt.mnl)goto 21041                                               3623
      if(flmin.ge.1.0)goto 21041                                           3624
      me=0                                                                 3624
21430 do 21431 j=1,nin                                                     3624
      if(ao(j,1,m).ne.0.0) me=me+1                                         3624
21431 continue                                                             3624
21432 continue                                                             3624
      if(me.gt.ne)goto 21042                                               3625
      if(rsq0-rsq.lt.sml*rsq)goto 21042                                    3625
      if(rsqo(m).gt.rsqmax)goto 21042                                      3626
21041 continue                                                             3627
21042 continue                                                             3627
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3628
      return                                                               3629
      end                                                                  3630
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   3632 
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),cl(2,ni)     3633
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(ni)            3634
      integer ju(ni),m(nx),kin(nlam)                                       3635
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del                    
      integer, dimension (:), allocatable :: mm,is,ixx,isc                      
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr;                         
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3644
      exmn=-exmx                                                           3645
      allocate(mm(1:ni),stat=ierr)                                         3645
      jerr=jerr+ierr                                                       3646
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3646
      jerr=jerr+ierr                                                       3647
      allocate(sxp(1:no),stat=ierr)                                        3647
      jerr=jerr+ierr                                                       3648
      allocate(sxpl(1:no),stat=ierr)                                       3648
      jerr=jerr+ierr                                                       3649
      allocate(ga(1:ni),stat=ierr)                                         3649
      jerr=jerr+ierr                                                       3650
      allocate(ixx(1:ni),stat=ierr)                                        3650
      jerr=jerr+ierr                                                       3651
      allocate(gk(1:nc),stat=ierr)                                         3651
      jerr=jerr+ierr                                                       3652
      allocate(del(1:nc),stat=ierr)                                        3652
      jerr=jerr+ierr                                                       3653
      allocate(isc(1:nc),stat=ierr)                                        3653
      jerr=jerr+ierr                                                       3654
      if(jerr.ne.0) return                                                 3655
      pmax=1.0-pmin                                                        3655
      emin=pmin/pmax                                                       3655
      emax=1.0/emin                                                        3656
      bta=parm                                                             3656
      omb=1.0-bta                                                          3656
      dev1=0.0                                                             3656
      dev0=0.0                                                             3657
21440 do 21441 ic=1,nc                                                     3657
      q0=dot_product(w,y(:,ic))                                            3658
      if(q0 .gt. pmin)goto 21461                                           3658
      jerr =8000+ic                                                        3658
      return                                                               3658
21461 continue                                                             3659
      if(q0 .lt. pmax)goto 21481                                           3659
      jerr =9000+ic                                                        3659
      return                                                               3659
21481 continue                                                             3660
      if(intr .ne. 0)goto 21501                                            3660
      q0=1.0/nc                                                            3660
      b(0,ic)=0.0                                                          3660
      goto 21511                                                           3661
21501 continue                                                             3661
      b(0,ic)=log(q0)                                                      3661
      dev1=dev1-q0*b(0,ic)                                                 3661
21511 continue                                                             3662
21491 continue                                                             3662
      b(1:ni,ic)=0.0                                                       3663
21441 continue                                                             3664
21442 continue                                                             3664
      if(intr.eq.0) dev1=log(float(nc))                                    3664
      ixx=0                                                                3664
      al=0.0                                                               3665
      if(nonzero(no*nc,g) .ne. 0)goto 21531                                3666
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3666
      sxp=0.0                                                              3667
21540 do 21541 ic=1,nc                                                     3667
      q(:,ic)=exp(b(0,ic))                                                 3667
      sxp=sxp+q(:,ic)                                                      3667
21541 continue                                                             3668
21542 continue                                                             3668
      goto 21551                                                           3669
21531 continue                                                             3669
21560 do 21561 i=1,no                                                      3669
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3669
21561 continue                                                             3669
21562 continue                                                             3669
      sxp=0.0                                                              3670
      if(intr .ne. 0)goto 21581                                            3670
      b(0,:)=0.0                                                           3670
      goto 21591                                                           3671
21581 continue                                                             3671
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3671
      if(jerr.ne.0) return                                                 3671
21591 continue                                                             3672
21571 continue                                                             3672
      dev1=0.0                                                             3673
21600 do 21601 ic=1,nc                                                     3673
      q(:,ic)=b(0,ic)+g(:,ic)                                              3674
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3675
      q(:,ic)=exp(q(:,ic))                                                 3675
      sxp=sxp+q(:,ic)                                                      3676
21601 continue                                                             3677
21602 continue                                                             3677
      sxpl=w*log(sxp)                                                      3677
21610 do 21611 ic=1,nc                                                     3677
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3677
21611 continue                                                             3678
21612 continue                                                             3678
21551 continue                                                             3679
21521 continue                                                             3679
21620 do 21621 ic=1,nc                                                     3679
21630 do 21631 i=1,no                                                      3679
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3679
21631 continue                                                             3679
21632 continue                                                             3679
21621 continue                                                             3680
21622 continue                                                             3680
      dev0=dev0+dev1                                                       3681
      if(flmin .ge. 1.0)goto 21651                                         3681
      eqs=max(eps,flmin)                                                   3681
      alf=eqs**(1.0/(nlam-1))                                              3681
21651 continue                                                             3682
      m=0                                                                  3682
      mm=0                                                                 3682
      nin=0                                                                3682
      nlp=0                                                                3682
      mnl=min(mnlam,nlam)                                                  3682
      bs=0.0                                                               3682
      shr=shri*dev0                                                        3683
      ga=0.0                                                               3684
21660 do 21661 ic=1,nc                                                     3684
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3685
21670 do 21671 j=1,ni                                                      3685
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            3685
21671 continue                                                             3686
21672 continue                                                             3686
21661 continue                                                             3687
21662 continue                                                             3687
      ga=sqrt(ga)                                                          3688
21680 do 21681 ilm=1,nlam                                                  3688
      al0=al                                                               3689
      if(flmin .lt. 1.0)goto 21701                                         3689
      al=ulam(ilm)                                                         3689
      goto 21691                                                           3690
21701 if(ilm .le. 2)goto 21711                                             3690
      al=al*alf                                                            3690
      goto 21691                                                           3691
21711 if(ilm .ne. 1)goto 21721                                             3691
      al=big                                                               3691
      goto 21731                                                           3692
21721 continue                                                             3692
      al0=0.0                                                              3693
21740 do 21741 j=1,ni                                                      3693
      if(ju(j).eq.0)goto 21741                                             3693
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3693
21741 continue                                                             3694
21742 continue                                                             3694
      al0=al0/max(bta,1.0e-3)                                              3694
      al=alf*al0                                                           3695
21731 continue                                                             3696
21691 continue                                                             3696
      al2=al*omb                                                           3696
      al1=al*bta                                                           3696
      tlam=bta*(2.0*al-al0)                                                3697
21750 do 21751 k=1,ni                                                      3697
      if(ixx(k).eq.1)goto 21751                                            3697
      if(ju(k).eq.0)goto 21751                                             3698
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3699
21751 continue                                                             3700
21752 continue                                                             3700
10880 continue                                                             3701
21760 continue                                                             3701
21761 continue                                                             3701
      ix=0                                                                 3701
      jx=ix                                                                3701
      kx=jx                                                                3701
      t=0.0                                                                3702
21770 do 21771 ic=1,nc                                                     3702
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3702
21771 continue                                                             3703
21772 continue                                                             3703
      if(t .ge. eps)goto 21791                                             3703
      kx=1                                                                 3703
      goto 21762                                                           3703
21791 continue                                                             3703
      t=2.0*t                                                              3703
      alt=al1/t                                                            3703
      al2t=al2/t                                                           3704
21800 do 21801 ic=1,nc                                                     3705
      bs(0,ic)=b(0,ic)                                                     3705
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3706
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3707
      d=0.0                                                                3707
      if(intr.ne.0) d=sum(r(:,ic))                                         3708
      if(d .eq. 0.0)goto 21821                                             3709
      b(0,ic)=b(0,ic)+d                                                    3709
      r(:,ic)=r(:,ic)-d*w                                                  3709
      dlx=max(dlx,d**2)                                                    3710
21821 continue                                                             3711
21801 continue                                                             3712
21802 continue                                                             3712
21830 continue                                                             3712
21831 continue                                                             3712
      nlp=nlp+nc                                                           3712
      dlx=0.0                                                              3713
21840 do 21841 k=1,ni                                                      3713
      if(ixx(k).eq.0)goto 21841                                            3713
      gkn=0.0                                                              3714
21850 do 21851 ic=1,nc                                                     3714
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3715
      gkn=gkn+gk(ic)**2                                                    3716
21851 continue                                                             3717
21852 continue                                                             3717
      gkn=sqrt(gkn)                                                        3717
      u=1.0-alt*vp(k)/gkn                                                  3717
      del=b(k,:)                                                           3718
      if(u .gt. 0.0)goto 21871                                             3718
      b(k,:)=0.0                                                           3718
      goto 21881                                                           3719
21871 continue                                                             3719
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3720
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   3722 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3723
21881 continue                                                             3724
21861 continue                                                             3724
      del=b(k,:)-del                                                       3724
      if(maxval(abs(del)).le.0.0)goto 21841                                3725
21890 do 21891 ic=1,nc                                                     3725
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3726
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3727
21891 continue                                                             3728
21892 continue                                                             3728
      if(mm(k) .ne. 0)goto 21911                                           3728
      nin=nin+1                                                            3729
      if(nin .le. nx)goto 21931                                            3729
      jx=1                                                                 3729
      goto 21842                                                           3729
21931 continue                                                             3730
      mm(k)=nin                                                            3730
      m(nin)=k                                                             3731
21911 continue                                                             3732
21841 continue                                                             3733
21842 continue                                                             3733
      if(jx.gt.0)goto 21832                                                3733
      if(dlx.lt.shr)goto 21832                                             3734
      if(nlp .le. maxit)goto 21951                                         3734
      jerr=-ilm                                                            3734
      return                                                               3734
21951 continue                                                             3735
21960 continue                                                             3735
21961 continue                                                             3735
      nlp=nlp+nc                                                           3735
      dlx=0.0                                                              3736
21970 do 21971 l=1,nin                                                     3736
      k=m(l)                                                               3736
      gkn=0.0                                                              3737
21980 do 21981 ic=1,nc                                                     3737
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3738
      gkn=gkn+gk(ic)**2                                                    3739
21981 continue                                                             3740
21982 continue                                                             3740
      gkn=sqrt(gkn)                                                        3740
      u=1.0-alt*vp(k)/gkn                                                  3740
      del=b(k,:)                                                           3741
      if(u .gt. 0.0)goto 22001                                             3741
      b(k,:)=0.0                                                           3741
      goto 22011                                                           3742
22001 continue                                                             3742
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3743
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   3745 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3746
22011 continue                                                             3747
21991 continue                                                             3747
      del=b(k,:)-del                                                       3747
      if(maxval(abs(del)).le.0.0)goto 21971                                3748
22020 do 22021 ic=1,nc                                                     3748
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3749
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3750
22021 continue                                                             3751
22022 continue                                                             3751
21971 continue                                                             3752
21972 continue                                                             3752
      if(dlx.lt.shr)goto 21962                                             3752
      if(nlp .le. maxit)goto 22041                                         3752
      jerr=-ilm                                                            3752
      return                                                               3752
22041 continue                                                             3754
      goto 21961                                                           3755
21962 continue                                                             3755
      goto 21831                                                           3756
21832 continue                                                             3756
      if(jx.gt.0)goto 21762                                                3757
22050 do 22051 ic=1,nc                                                     3758
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                3759
      if(ix .ne. 0)goto 22071                                              3760
22080 do 22081 j=1,nin                                                     3760
      k=m(j)                                                               3761
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22101                   3761
      ix=1                                                                 3761
      goto 22082                                                           3761
22101 continue                                                             3763
22081 continue                                                             3764
22082 continue                                                             3764
22071 continue                                                             3765
22110 do 22111 i=1,no                                                      3765
      fi=b(0,ic)+g(i,ic)                                                   3767
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         3768
      fi=min(max(exmn,fi),exmx)                                            3768
      sxp(i)=sxp(i)-q(i,ic)                                                3769
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    3770
      sxp(i)=sxp(i)+q(i,ic)                                                3771
22111 continue                                                             3772
22112 continue                                                             3772
22051 continue                                                             3773
22052 continue                                                             3773
      s=-sum(b(0,:))/nc                                                    3773
      b(0,:)=b(0,:)+s                                                      3774
      if(jx.gt.0)goto 21762                                                3775
      if(ix .ne. 0)goto 22131                                              3776
22140 do 22141 k=1,ni                                                      3776
      if(ixx(k).eq.1)goto 22141                                            3776
      if(ju(k).eq.0)goto 22141                                             3776
      ga(k)=0.0                                                            3776
22141 continue                                                             3777
22142 continue                                                             3777
22150 do 22151 ic=1,nc                                                     3777
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3778
22160 do 22161 k=1,ni                                                      3778
      if(ixx(k).eq.1)goto 22161                                            3778
      if(ju(k).eq.0)goto 22161                                             3779
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           3780
22161 continue                                                             3781
22162 continue                                                             3781
22151 continue                                                             3782
22152 continue                                                             3782
      ga=sqrt(ga)                                                          3783
22170 do 22171 k=1,ni                                                      3783
      if(ixx(k).eq.1)goto 22171                                            3783
      if(ju(k).eq.0)goto 22171                                             3784
      if(ga(k) .le. al1*vp(k))goto 22191                                   3784
      ixx(k)=1                                                             3784
      ix=1                                                                 3784
22191 continue                                                             3785
22171 continue                                                             3786
22172 continue                                                             3786
      if(ix.eq.1) go to 10880                                              3787
      goto 21762                                                           3788
22131 continue                                                             3789
      goto 21761                                                           3790
21762 continue                                                             3790
      if(kx .le. 0)goto 22211                                              3790
      jerr=-20000-ilm                                                      3790
      goto 21682                                                           3790
22211 continue                                                             3791
      if(jx .le. 0)goto 22231                                              3791
      jerr=-10000-ilm                                                      3791
      goto 21682                                                           3791
22231 continue                                                             3791
      devi=0.0                                                             3792
22240 do 22241 ic=1,nc                                                     3793
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3793
      a0(ic,ilm)=b(0,ic)                                                   3794
22250 do 22251 i=1,no                                                      3794
      if(y(i,ic).le.0.0)goto 22251                                         3795
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3796
22251 continue                                                             3797
22252 continue                                                             3797
22241 continue                                                             3798
22242 continue                                                             3798
      kin(ilm)=nin                                                         3798
      alm(ilm)=al                                                          3798
      lmu=ilm                                                              3799
      dev(ilm)=(dev1-devi)/dev0                                            3800
      if(ilm.lt.mnl)goto 21681                                             3800
      if(flmin.ge.1.0)goto 21681                                           3801
      me=0                                                                 3801
22260 do 22261 j=1,nin                                                     3801
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3801
22261 continue                                                             3801
22262 continue                                                             3801
      if(me.gt.ne)goto 21682                                               3802
      if(dev(ilm).gt.devmax)goto 21682                                     3802
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21682                             3803
21681 continue                                                             3804
21682 continue                                                             3804
      g=log(q)                                                             3804
22270 do 22271 i=1,no                                                      3804
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3804
22271 continue                                                             3805
22272 continue                                                             3805
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    3806
      return                                                               3807
      end                                                                  3808
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   3810 
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni),   3811 
     *xv(ni)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          3812
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3813
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del,sc,svr             
      integer, dimension (:), allocatable :: mm,is,iy,isc                       
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3822
      exmn=-exmx                                                           3823
      allocate(mm(1:ni),stat=ierr)                                         3823
      jerr=jerr+ierr                                                       3824
      allocate(ga(1:ni),stat=ierr)                                         3824
      jerr=jerr+ierr                                                       3825
      allocate(gk(1:nc),stat=ierr)                                         3825
      jerr=jerr+ierr                                                       3826
      allocate(del(1:nc),stat=ierr)                                        3826
      jerr=jerr+ierr                                                       3827
      allocate(iy(1:ni),stat=ierr)                                         3827
      jerr=jerr+ierr                                                       3828
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3828
      jerr=jerr+ierr                                                       3829
      allocate(sxp(1:no),stat=ierr)                                        3829
      jerr=jerr+ierr                                                       3830
      allocate(sxpl(1:no),stat=ierr)                                       3830
      jerr=jerr+ierr                                                       3831
      allocate(svr(1:nc),stat=ierr)                                        3831
      jerr=jerr+ierr                                                       3832
      allocate(sc(1:no),stat=ierr)                                         3832
      jerr=jerr+ierr                                                       3833
      allocate(isc(1:nc),stat=ierr)                                        3833
      jerr=jerr+ierr                                                       3834
      if(jerr.ne.0) return                                                 3835
      pmax=1.0-pmin                                                        3835
      emin=pmin/pmax                                                       3835
      emax=1.0/emin                                                        3836
      bta=parm                                                             3836
      omb=1.0-bta                                                          3836
      dev1=0.0                                                             3836
      dev0=0.0                                                             3837
22280 do 22281 ic=1,nc                                                     3837
      q0=dot_product(w,y(:,ic))                                            3838
      if(q0 .gt. pmin)goto 22301                                           3838
      jerr =8000+ic                                                        3838
      return                                                               3838
22301 continue                                                             3839
      if(q0 .lt. pmax)goto 22321                                           3839
      jerr =9000+ic                                                        3839
      return                                                               3839
22321 continue                                                             3840
      b(1:ni,ic)=0.0                                                       3841
      if(intr .ne. 0)goto 22341                                            3841
      q0=1.0/nc                                                            3841
      b(0,ic)=0.0                                                          3841
      goto 22351                                                           3842
22341 continue                                                             3842
      b(0,ic)=log(q0)                                                      3842
      dev1=dev1-q0*b(0,ic)                                                 3842
22351 continue                                                             3843
22331 continue                                                             3843
22281 continue                                                             3844
22282 continue                                                             3844
      if(intr.eq.0) dev1=log(float(nc))                                    3844
      iy=0                                                                 3844
      al=0.0                                                               3845
      if(nonzero(no*nc,g) .ne. 0)goto 22371                                3846
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3846
      sxp=0.0                                                              3847
22380 do 22381 ic=1,nc                                                     3847
      q(:,ic)=exp(b(0,ic))                                                 3847
      sxp=sxp+q(:,ic)                                                      3847
22381 continue                                                             3848
22382 continue                                                             3848
      goto 22391                                                           3849
22371 continue                                                             3849
22400 do 22401 i=1,no                                                      3849
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3849
22401 continue                                                             3849
22402 continue                                                             3849
      sxp=0.0                                                              3850
      if(intr .ne. 0)goto 22421                                            3850
      b(0,:)=0.0                                                           3850
      goto 22431                                                           3851
22421 continue                                                             3851
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3851
      if(jerr.ne.0) return                                                 3851
22431 continue                                                             3852
22411 continue                                                             3852
      dev1=0.0                                                             3853
22440 do 22441 ic=1,nc                                                     3853
      q(:,ic)=b(0,ic)+g(:,ic)                                              3854
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3855
      q(:,ic)=exp(q(:,ic))                                                 3855
      sxp=sxp+q(:,ic)                                                      3856
22441 continue                                                             3857
22442 continue                                                             3857
      sxpl=w*log(sxp)                                                      3857
22450 do 22451 ic=1,nc                                                     3857
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3857
22451 continue                                                             3858
22452 continue                                                             3858
22391 continue                                                             3859
22361 continue                                                             3859
22460 do 22461 ic=1,nc                                                     3859
22470 do 22471 i=1,no                                                      3859
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3859
22471 continue                                                             3859
22472 continue                                                             3859
22461 continue                                                             3860
22462 continue                                                             3860
      dev0=dev0+dev1                                                       3861
      if(flmin .ge. 1.0)goto 22491                                         3861
      eqs=max(eps,flmin)                                                   3861
      alf=eqs**(1.0/(nlam-1))                                              3861
22491 continue                                                             3862
      m=0                                                                  3862
      mm=0                                                                 3862
      nin=0                                                                3862
      nlp=0                                                                3862
      mnl=min(mnlam,nlam)                                                  3862
      bs=0.0                                                               3863
      shr=shri*dev0                                                        3863
      ga=0.0                                                               3864
22500 do 22501 ic=1,nc                                                     3864
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3864
      svr(ic)=sum(r(:,ic))                                                 3865
22510 do 22511 j=1,ni                                                      3865
      if(ju(j).eq.0)goto 22511                                             3866
      jb=ix(j)                                                             3866
      je=ix(j+1)-1                                                         3867
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3868
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3869
22511 continue                                                             3870
22512 continue                                                             3870
22501 continue                                                             3871
22502 continue                                                             3871
      ga=sqrt(ga)                                                          3872
22520 do 22521 ilm=1,nlam                                                  3872
      al0=al                                                               3873
      if(flmin .lt. 1.0)goto 22541                                         3873
      al=ulam(ilm)                                                         3873
      goto 22531                                                           3874
22541 if(ilm .le. 2)goto 22551                                             3874
      al=al*alf                                                            3874
      goto 22531                                                           3875
22551 if(ilm .ne. 1)goto 22561                                             3875
      al=big                                                               3875
      goto 22571                                                           3876
22561 continue                                                             3876
      al0=0.0                                                              3877
22580 do 22581 j=1,ni                                                      3877
      if(ju(j).eq.0)goto 22581                                             3877
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3877
22581 continue                                                             3878
22582 continue                                                             3878
      al0=al0/max(bta,1.0e-3)                                              3878
      al=alf*al0                                                           3879
22571 continue                                                             3880
22531 continue                                                             3880
      al2=al*omb                                                           3880
      al1=al*bta                                                           3880
      tlam=bta*(2.0*al-al0)                                                3881
22590 do 22591 k=1,ni                                                      3881
      if(iy(k).eq.1)goto 22591                                             3881
      if(ju(k).eq.0)goto 22591                                             3882
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      3883
22591 continue                                                             3884
22592 continue                                                             3884
10880 continue                                                             3885
22600 continue                                                             3885
22601 continue                                                             3885
      ixx=0                                                                3885
      jxx=ixx                                                              3885
      kxx=jxx                                                              3885
      t=0.0                                                                3886
22610 do 22611 ic=1,nc                                                     3886
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3886
22611 continue                                                             3887
22612 continue                                                             3887
      if(t .ge. eps)goto 22631                                             3887
      kxx=1                                                                3887
      goto 22602                                                           3887
22631 continue                                                             3887
      t=2.0*t                                                              3887
      alt=al1/t                                                            3887
      al2t=al2/t                                                           3888
22640 do 22641 ic=1,nc                                                     3888
      bs(0,ic)=b(0,ic)                                                     3888
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3889
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3889
      svr(ic)=sum(r(:,ic))                                                 3890
      if(intr .eq. 0)goto 22661                                            3890
      b(0,ic)=b(0,ic)+svr(ic)                                              3890
      r(:,ic)=r(:,ic)-svr(ic)*w                                            3891
      dlx=max(dlx,svr(ic)**2)                                              3892
22661 continue                                                             3893
22641 continue                                                             3894
22642 continue                                                             3894
22670 continue                                                             3894
22671 continue                                                             3894
      nlp=nlp+nc                                                           3894
      dlx=0.0                                                              3895
22680 do 22681 k=1,ni                                                      3895
      if(iy(k).eq.0)goto 22681                                             3896
      jb=ix(k)                                                             3896
      je=ix(k+1)-1                                                         3896
      del=b(k,:)                                                           3896
      gkn=0.0                                                              3897
22690 do 22691 ic=1,nc                                                     3898
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        3899
      gk(ic)=u+del(ic)*xv(k)                                               3899
      gkn=gkn+gk(ic)**2                                                    3900
22691 continue                                                             3901
22692 continue                                                             3901
      gkn=sqrt(gkn)                                                        3901
      u=1.0-alt*vp(k)/gkn                                                  3902
      if(u .gt. 0.0)goto 22711                                             3902
      b(k,:)=0.0                                                           3902
      goto 22721                                                           3903
22711 continue                                                             3904
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3905
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   3907 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3908
22721 continue                                                             3909
22701 continue                                                             3909
      del=b(k,:)-del                                                       3909
      if(maxval(abs(del)).le.0.0)goto 22681                                3910
22730 do 22731 ic=1,nc                                                     3910
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3911
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3913 
     *b(k))/xs(k)
22731 continue                                                             3914
22732 continue                                                             3914
      if(mm(k) .ne. 0)goto 22751                                           3914
      nin=nin+1                                                            3915
      if(nin .le. nx)goto 22771                                            3915
      jxx=1                                                                3915
      goto 22682                                                           3915
22771 continue                                                             3916
      mm(k)=nin                                                            3916
      m(nin)=k                                                             3917
22751 continue                                                             3918
22681 continue                                                             3919
22682 continue                                                             3919
      if(jxx.gt.0)goto 22672                                               3920
      if(dlx.lt.shr)goto 22672                                             3920
      if(nlp .le. maxit)goto 22791                                         3920
      jerr=-ilm                                                            3920
      return                                                               3920
22791 continue                                                             3921
22800 continue                                                             3921
22801 continue                                                             3921
      nlp=nlp+nc                                                           3921
      dlx=0.0                                                              3922
22810 do 22811 l=1,nin                                                     3922
      k=m(l)                                                               3922
      jb=ix(k)                                                             3922
      je=ix(k+1)-1                                                         3922
      del=b(k,:)                                                           3922
      gkn=0.0                                                              3923
22820 do 22821 ic=1,nc                                                     3924
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      3926
      gk(ic)=u+del(ic)*xv(k)                                               3926
      gkn=gkn+gk(ic)**2                                                    3927
22821 continue                                                             3928
22822 continue                                                             3928
      gkn=sqrt(gkn)                                                        3928
      u=1.0-alt*vp(k)/gkn                                                  3929
      if(u .gt. 0.0)goto 22841                                             3929
      b(k,:)=0.0                                                           3929
      goto 22851                                                           3930
22841 continue                                                             3931
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3932
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   3934 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3935
22851 continue                                                             3936
22831 continue                                                             3936
      del=b(k,:)-del                                                       3936
      if(maxval(abs(del)).le.0.0)goto 22811                                3937
22860 do 22861 ic=1,nc                                                     3937
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3938
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3940 
     *b(k))/xs(k)
22861 continue                                                             3941
22862 continue                                                             3941
22811 continue                                                             3942
22812 continue                                                             3942
      if(dlx.lt.shr)goto 22802                                             3942
      if(nlp .le. maxit)goto 22881                                         3942
      jerr=-ilm                                                            3942
      return                                                               3942
22881 continue                                                             3944
      goto 22801                                                           3945
22802 continue                                                             3945
      goto 22671                                                           3946
22672 continue                                                             3946
      if(jxx.gt.0)goto 22602                                               3947
22890 do 22891 ic=1,nc                                                     3948
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               3949
      if(ixx .ne. 0)goto 22911                                             3950
22920 do 22921 j=1,nin                                                     3950
      k=m(j)                                                               3951
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22941                   3951
      ixx=1                                                                3951
      goto 22922                                                           3951
22941 continue                                                             3953
22921 continue                                                             3954
22922 continue                                                             3954
22911 continue                                                             3955
      sc=b(0,ic)+g(:,ic)                                                   3955
      b0=0.0                                                               3956
22950 do 22951 j=1,nin                                                     3956
      l=m(j)                                                               3956
      jb=ix(l)                                                             3956
      je=ix(l+1)-1                                                         3957
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   3958
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            3959
22951 continue                                                             3960
22952 continue                                                             3960
      sc=min(max(exmn,sc+b0),exmx)                                         3961
      sxp=sxp-q(:,ic)                                                      3962
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          3963
      sxp=sxp+q(:,ic)                                                      3964
22891 continue                                                             3965
22892 continue                                                             3965
      s=sum(b(0,:))/nc                                                     3965
      b(0,:)=b(0,:)-s                                                      3966
      if(jxx.gt.0)goto 22602                                               3967
      if(ixx .ne. 0)goto 22971                                             3968
22980 do 22981 j=1,ni                                                      3968
      if(iy(j).eq.1)goto 22981                                             3968
      if(ju(j).eq.0)goto 22981                                             3968
      ga(j)=0.0                                                            3968
22981 continue                                                             3969
22982 continue                                                             3969
22990 do 22991 ic=1,nc                                                     3969
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3970
23000 do 23001 j=1,ni                                                      3970
      if(iy(j).eq.1)goto 23001                                             3970
      if(ju(j).eq.0)goto 23001                                             3971
      jb=ix(j)                                                             3971
      je=ix(j+1)-1                                                         3972
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3973
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3974
23001 continue                                                             3975
23002 continue                                                             3975
22991 continue                                                             3976
22992 continue                                                             3976
      ga=sqrt(ga)                                                          3977
23010 do 23011 k=1,ni                                                      3977
      if(iy(k).eq.1)goto 23011                                             3977
      if(ju(k).eq.0)goto 23011                                             3978
      if(ga(k) .le. al1*vp(k))goto 23031                                   3978
      iy(k)=1                                                              3978
      ixx=1                                                                3978
23031 continue                                                             3979
23011 continue                                                             3980
23012 continue                                                             3980
      if(ixx.eq.1) go to 10880                                             3981
      goto 22602                                                           3982
22971 continue                                                             3983
      goto 22601                                                           3984
22602 continue                                                             3984
      if(kxx .le. 0)goto 23051                                             3984
      jerr=-20000-ilm                                                      3984
      goto 22522                                                           3984
23051 continue                                                             3985
      if(jxx .le. 0)goto 23071                                             3985
      jerr=-10000-ilm                                                      3985
      goto 22522                                                           3985
23071 continue                                                             3985
      devi=0.0                                                             3986
23080 do 23081 ic=1,nc                                                     3987
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3987
      a0(ic,ilm)=b(0,ic)                                                   3988
23090 do 23091 i=1,no                                                      3988
      if(y(i,ic).le.0.0)goto 23091                                         3989
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3990
23091 continue                                                             3991
23092 continue                                                             3991
23081 continue                                                             3992
23082 continue                                                             3992
      kin(ilm)=nin                                                         3992
      alm(ilm)=al                                                          3992
      lmu=ilm                                                              3993
      dev(ilm)=(dev1-devi)/dev0                                            3994
      if(ilm.lt.mnl)goto 22521                                             3994
      if(flmin.ge.1.0)goto 22521                                           3995
      me=0                                                                 3995
23100 do 23101 j=1,nin                                                     3995
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3995
23101 continue                                                             3995
23102 continue                                                             3995
      if(me.gt.ne)goto 22522                                               3996
      if(dev(ilm).gt.devmax)goto 22522                                     3996
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22522                             3997
22521 continue                                                             3998
22522 continue                                                             3998
      g=log(q)                                                             3998
23110 do 23111 i=1,no                                                      3998
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3998
23111 continue                                                             3999
23112 continue                                                             3999
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  4000
      return                                                               4001
      end                                                                  4002
      subroutine psort7 (v,a,ii,jj)                                             
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      real v                                                                    
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       

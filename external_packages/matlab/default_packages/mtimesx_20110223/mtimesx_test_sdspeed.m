% Test routine for mtimesx, op(single) * op(double) speed vs MATLAB
%******************************************************************************
%
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
%
%  Function:    mtimesx_test_sdspeed
%  Filename:    mtimesx_test_sdspeed.m
%  Programmer:  James Tursa
%  Version:     1.0
%  Date:        September 27, 2009
%  Copyright:   (c) 2009 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are
%  met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%  POSSIBILITY OF SUCH DAMAGE.
%
%  Syntax (arguments in brackets [ ] are optional):
%
%    T = mtimesx_test_ddspeed( [N [,D]] )
%
%  Inputs:
%
%    N = Number of runs to make for each individual test. The test result will
%        be the median of N runs. N must be even. If N is odd, it will be
%        automatically increased to the next even number. The default is 10,
%        which can take *hours* to run. Best to run this program overnight.
%    D = The string 'details'. If present, this will cause all of the
%        individual intermediate run results to print as they happen.
%
%  Output:
%
%    T = A character array containing a summary of the results.
%
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function ttable = mtimesx_test_sdspeed(nn,details)
                                                                                                                                                                                                                                                                                                            
global mtimesx_ttable
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('****************************************************************************');
disp('*                                                                          *');
disp('*  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  *');
disp('*                                                                          *');
disp('*  This test program can take several *hours* to complete, particularly    *');
disp('*  when using the default number of runs as 10. It is strongly suggested   *');
disp('*  to close all applications and run this program overnight to get the     *');
disp('*  best possible result with minimal impacts to your computer usage.       *');
disp('*                                                                          *');
disp('*  The program will be done when you see the message:  DONE !              *');
disp('*                                                                          *');
disp('****************************************************************************');
disp(' ');
try
    input('Press Enter to start test, or Ctrl-C to exit ','s');
catch
    ttable = '';
    return
end
                                                                                                                                                                                                                                                                                                            
start_time = datenum(clock);
                                                                                                                                                                                                                                                                                                            
if nargin >= 1
    n = nn;
else
    n = 10;
end
if nargin < 2
    details = false;
else
    if( isempty(details) )  % code to get rid of the lint message
        details = true;
    else
        details = true;
    end
end
                                                                                                                                                                                                                                                                                                            
RC = '                                Real*Real  Real*Cplx  Cplx*Real  Cplx*Cplx';
                                                                                                                                                                                                                                                                                                            
compver = [computer ', ' version ', mtimesx mode ' mtimesx ', median of ' num2str(n) ' runs'];
k = length(compver);
                                                                                                                                                                                                                                                                                                            
mtimesx_ttable = char([]);
mtimesx_ttable(100,74) = ' ';
mtimesx_ttable(1,1:k) = compver;
mtimesx_ttable(2,:) = RC;
for r=3:170
mtimesx_ttable(r,:) = '                                       --         --         --         --';
end
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(compver);
disp('Test program for function mtimesx:')
disp('----------------------------------');
                                                                                                                                                                                                                                                                                                            
rsave = 2;
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeNN('Scalar  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1);
maxtimeNN('Vector  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeNN('Scalar  * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeNN('Array   * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1);
maxtimeNN('Vector  i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500);
maxtimeNN('Vector  o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeNN('Vector  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeNN('Matrix  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeNN('Matrix  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNN('Scalar  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNN('Vector  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeNN('Scalar  * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNN('Array   * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeNN('Vector  i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeNN('Vector  o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNN('Vector  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeNN('Matrix  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNN('Matrix  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeNN('Scalar  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1);
maxtimeNN('Vector  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeNN('Scalar  * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1);
maxtimeNN('Array   * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1);
maxtimeNN('Vector  i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500);
maxtimeNN('Vector  o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeNN('Vector  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeNN('Matrix  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeNN('Matrix  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNN('Scalar  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNN('Vector  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeNN('Scalar  * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNN('Array   * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeNN('Vector  i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeNN('Vector  o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNN('Vector  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeNN('Matrix  * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNN('Matrix  * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)  * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeNT('Scalar  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeNT('Vector  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeNT('Array   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000);
maxtimeNT('Vector  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1);
maxtimeNT('Vector  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeNT('Vector  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeNT('Matrix  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeNT('Matrix  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNT('Scalar  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNT('Vector  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNT('Array   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeNT('Vector  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeNT('Vector  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNT('Vector  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeNT('Matrix  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNT('Matrix  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)  * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeNT('Scalar  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeNT('Vector  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = rand(1,1);
maxtimeNT('Array   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000);
maxtimeNT('Vector  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1);
maxtimeNT('Vector  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeNT('Vector  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeNT('Matrix  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeNT('Matrix  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNT('Scalar  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNT('Vector  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNT('Array   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeNT('Vector  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeNT('Vector  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNT('Vector  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeNT('Matrix  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNT('Matrix  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)  * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeNC('Scalar  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeNC('Vector  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeNC('Array   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000);
maxtimeNC('Vector  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1);
maxtimeNC('Vector  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeNC('Vector  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeNC('Matrix  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeNC('Matrix  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNC('Scalar  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNC('Vector  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNC('Array   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeNC('Vector  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeNC('Vector  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNC('Vector  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeNC('Matrix  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNC('Matrix  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)  * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeNC('Scalar  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeNC('Vector  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = rand(1,1);
maxtimeNC('Array   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000);
maxtimeNC('Vector  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1);
maxtimeNC('Vector  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeNC('Vector  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeNC('Matrix  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeNC('Matrix  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNC('Scalar  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNC('Vector  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNC('Array   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeNC('Vector  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeNC('Vector  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNC('Vector  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeNC('Matrix  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNC('Matrix  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real) * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeNG('Scalar  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1);
maxtimeNG('Vector  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeNG('Scalar  * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeNG('Array   * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1);
maxtimeNG('Vector  i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500);
maxtimeNG('Vector  o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeNG('Vector  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeNG('Matrix  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeNG('Matrix  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real) * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNG('Scalar  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNG('Vector  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeNG('Scalar  * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeNG('Array   * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeNG('Vector  i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeNG('Vector  o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNG('Vector  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeNG('Matrix  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNG('Matrix  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeNG('Scalar  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1);
maxtimeNG('Vector  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeNG('Scalar  * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1);
maxtimeNG('Array   * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1);
maxtimeNG('Vector  i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500);
maxtimeNG('Vector  o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeNG('Vector  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeNG('Matrix  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeNG('Matrix  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeNG('Scalar  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNG('Vector  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeNG('Scalar  * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeNG('Array   * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeNG('Vector  i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeNG('Vector  o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNG('Vector  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeNG('Matrix  * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeNG('Matrix  * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real).'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeTN('Scalar.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeTN('Vector.'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeTN('Scalar.'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1);
maxtimeTN('Vector.'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500);
maxtimeTN('Vector.'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeTN('Vector.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeTN('Matrix.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeTN('Matrix.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real).'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTN('Scalar.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeTN('Vector.'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeTN('Scalar.'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeTN('Vector.'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeTN('Vector.'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTN('Vector.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeTN('Matrix.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTN('Matrix.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeTN('Scalar.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeTN('Vector.'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeTN('Scalar.'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1);
maxtimeTN('Vector.'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500);
maxtimeTN('Vector.'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeTN('Vector.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeTN('Matrix.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeTN('Matrix.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTN('Scalar.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeTN('Vector.'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeTN('Scalar.'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeTN('Vector.'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeTN('Vector.'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTN('Vector.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeTN('Matrix.'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTN('Matrix.'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real).'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeTT('Scalar.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeTT('Vector.'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000);
maxtimeTT('Vector.'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1);
maxtimeTT('Vector.'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeTT('Vector.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeTT('Matrix.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeTT('Matrix.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real).'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTT('Scalar.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeTT('Vector.'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeTT('Vector.'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeTT('Vector.'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTT('Vector.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeTT('Matrix.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTT('Matrix.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeTT('Scalar.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeTT('Vector.'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000);
maxtimeTT('Vector.'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1);
maxtimeTT('Vector.'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeTT('Vector.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeTT('Matrix.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeTT('Matrix.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTT('Scalar.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeTT('Vector.'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeTT('Vector.'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeTT('Vector.'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTT('Vector.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeTT('Matrix.'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTT('Matrix.'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real).'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeTC('Scalar.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeTC('Vector.'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000);
maxtimeTC('Vector.'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1);
maxtimeTC('Vector.'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeTC('Vector.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeTC('Matrix.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeTC('Matrix.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real).'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTC('Scalar.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeTC('Vector.'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeTC('Vector.'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeTC('Vector.'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTC('Vector.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeTC('Matrix.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTC('Matrix.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeTC('Scalar.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeTC('Vector.'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000);
maxtimeTC('Vector.'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1);
maxtimeTC('Vector.'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeTC('Vector.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeTC('Matrix.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeTC('Matrix.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTC('Scalar.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeTC('Vector.'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeTC('Vector.'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeTC('Vector.'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTC('Vector.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeTC('Matrix.'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTC('Matrix.'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real).'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeTG('Scalar.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeTG('Vector.'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeTG('Scalar.'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1);
maxtimeTG('Vector.'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500);
maxtimeTG('Vector.'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeTG('Vector.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeTG('Matrix.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeTG('Matrix.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real).'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTG('Scalar.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeTG('Vector.'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeTG('Scalar.'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeTG('Vector.'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeTG('Vector.'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTG('Vector.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeTG('Matrix.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTG('Matrix.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeTG('Scalar.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeTG('Vector.'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeTG('Scalar.'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1);
maxtimeTG('Vector.'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500);
maxtimeTG('Vector.'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeTG('Vector.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeTG('Matrix.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeTG('Matrix.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeTG('Scalar.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeTG('Vector.'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeTG('Scalar.'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeTG('Vector.'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeTG('Vector.'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTG('Vector.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeTG('Matrix.'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeTG('Matrix.'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeCN('Scalar'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeCN('Vector'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeCN('Scalar'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1);
maxtimeCN('Vector'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500);
maxtimeCN('Vector'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeCN('Vector'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeCN('Matrix'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeCN('Matrix'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCN('Scalar'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeCN('Vector'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeCN('Scalar'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeCN('Vector'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeCN('Vector'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCN('Vector'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeCN('Matrix'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCN('Matrix'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeCN('Scalar'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeCN('Vector'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeCN('Scalar'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1);
maxtimeCN('Vector'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500);
maxtimeCN('Vector'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeCN('Vector'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeCN('Matrix'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeCN('Matrix'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCN('Scalar'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeCN('Vector'' * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeCN('Scalar'' * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeCN('Vector'' i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeCN('Vector'' o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCN('Vector'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeCN('Matrix'' * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCN('Matrix'' * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeCT('Scalar'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeCT('Vector'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000);
maxtimeCT('Vector'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1);
maxtimeCT('Vector'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeCT('Vector'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeCT('Matrix'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeCT('Matrix'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCT('Scalar'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeCT('Vector'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeCT('Vector'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeCT('Vector'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCT('Vector'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeCT('Matrix'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCT('Matrix'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeCT('Scalar'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeCT('Vector'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000);
maxtimeCT('Vector'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1);
maxtimeCT('Vector'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeCT('Vector'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeCT('Matrix'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeCT('Matrix'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCT('Scalar'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeCT('Vector'' * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeCT('Vector'' i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeCT('Vector'' o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCT('Vector'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeCT('Matrix'' * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCT('Matrix'' * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeCC('Scalar'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeCC('Vector'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000);
maxtimeCC('Vector'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1);
maxtimeCC('Vector'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeCC('Vector'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeCC('Matrix'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeCC('Matrix'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCC('Scalar'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeCC('Vector'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeCC('Vector'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeCC('Vector'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCC('Vector'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeCC('Matrix'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCC('Matrix'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeCC('Scalar'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeCC('Vector'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000);
maxtimeCC('Vector'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1);
maxtimeCC('Vector'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeCC('Vector'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeCC('Matrix'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeCC('Matrix'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCC('Scalar'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeCC('Vector'' * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeCC('Vector'' i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeCC('Vector'' o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCC('Vector'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeCC('Matrix'' * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCC('Matrix'' * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('(real)'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeCG('Scalar'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeCG('Vector'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeCG('Scalar'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1);
maxtimeCG('Vector'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500);
maxtimeCG('Vector'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000);
maxtimeCG('Vector'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeCG('Matrix'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeCG('Matrix'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(real)'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCG('Scalar'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeCG('Vector'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeCG('Scalar'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeCG('Vector'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeCG('Vector'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCG('Vector'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeCG('Matrix'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCG('Matrix'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeCG('Scalar'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeCG('Vector'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeCG('Scalar'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1);
maxtimeCG('Vector'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500);
maxtimeCG('Vector'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000);
maxtimeCG('Vector'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeCG('Matrix'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeCG('Matrix'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeCG('Scalar'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeCG('Vector'' * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeCG('Scalar'' * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeCG('Vector'' i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeCG('Vector'' o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,1) + rand(2000,1)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCG('Vector'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeCG('Matrix'' * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeCG('Matrix'' * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('conj(real)   * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeGN('conj(Scalar) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1);
maxtimeGN('conj(Vector) * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeGN('conj(Scalar) * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeGN('conj(Array)  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1);
maxtimeGN('conj(Vector) i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500);
maxtimeGN('conj(Vector) o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeGN('conj(Vector) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeGN('conj(Matrix) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeGN('conj(Matrix) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('conj(real)   * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGN('conj(Scalar) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGN('conj(Vector) * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeGN('conj(Scalar) * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGN('conj(Array)  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeGN('conj(Vector) i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeGN('conj(Vector) o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGN('conj(Vector) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeGN('conj(Matrix) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGN('conj(Matrix) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeGN('conj(Scalar) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1);
maxtimeGN('conj(Vector) * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeGN('conj(Scalar) * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1);
maxtimeGN('conj(Array)  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1);
maxtimeGN('conj(Vector) i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500);
maxtimeGN('conj(Vector) o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeGN('conj(Vector) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeGN('conj(Matrix) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeGN('conj(Matrix) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGN('conj(Scalar) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGN('conj(Vector) * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeGN('conj(Scalar) * Array  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGN('conj(Array)  * Scalar ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeGN('conj(Vector) i Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeGN('conj(Vector) o Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGN('conj(Vector) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeGN('conj(Matrix) * Vector ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGN('conj(Matrix) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('conj(real)  * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeGT('conj(Scalar)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeGT('conj(Vector)  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeGT('conj(Array)   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000);
maxtimeGT('conj(Vector)  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1);
maxtimeGT('conj(Vector)  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeGT('conj(Vector)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeGT('conj(Matrix)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeGT('conj(Matrix)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('conj(real)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGT('conj(Scalar)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGT('conj(Vector)  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGT('conj(Array)   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeGT('conj(Vector)  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeGT('conj(Vector)  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGT('conj(Vector)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeGT('conj(Matrix)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGT('conj(Matrix)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)  * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeGT('conj(Scalar)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeGT('conj(Vector)  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = rand(1,1);
maxtimeGT('conj(Array)   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000);
maxtimeGT('conj(Vector)  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1);
maxtimeGT('conj(Vector)  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeGT('conj(Vector)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeGT('conj(Matrix)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeGT('conj(Matrix)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGT('conj(Scalar)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGT('conj(Vector)  * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGT('conj(Array)   * Scalar.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeGT('conj(Vector)  i Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeGT('conj(Vector)  o Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGT('conj(Vector)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeGT('conj(Matrix)  * Vector.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGT('conj(Matrix)  * Matrix.'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('conj(real)    * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeGC('conj(Scalar)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1);
maxtimeGC('conj(Vector)  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeGC('conj(Array)   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000);
maxtimeGC('conj(Vector)  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1);
maxtimeGC('conj(Vector)  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeGC('conj(Vector)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000);
maxtimeGC('conj(Matrix)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeGC('conj(Matrix)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('conj(real)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGC('conj(Scalar)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGC('conj(Vector)  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGC('conj(Array)   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeGC('conj(Vector)  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeGC('conj(Vector)  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGC('conj(Vector)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeGC('conj(Matrix)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGC('conj(Matrix)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)  * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeGC('conj(Scalar)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1);
maxtimeGC('conj(Vector)  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = rand(1,1);
maxtimeGC('conj(Array)   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000);
maxtimeGC('conj(Vector)  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1);
maxtimeGC('conj(Vector)  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeGC('conj(Vector)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000);
maxtimeGC('conj(Matrix)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeGC('conj(Matrix)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGC('conj(Scalar)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000000,1) + rand(1000000,1)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGC('conj(Vector)  * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGC('conj(Array)   * Scalar'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(1,10000000) + rand(1,10000000)*1i;
maxtimeGC('conj(Vector)  i Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(2500,1) + rand(2500,1)*1i;
maxtimeGC('conj(Vector)  o Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGC('conj(Vector)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(1,2000) + rand(1,2000)*1i;
maxtimeGC('conj(Matrix)  * Vector'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGC('conj(Matrix)  * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs']);
disp(' ');
disp('conj(real)   * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000);
maxtimeGG('conj(Scalar) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1);
maxtimeGG('conj(Vector) * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400);
maxtimeGG('conj(Scalar) * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1);
maxtimeGG('conj(Array)  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1);
maxtimeGG('conj(Vector) i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500);
maxtimeGG('conj(Vector) o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000);
maxtimeGG('conj(Vector) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1);
maxtimeGG('conj(Matrix) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000);
maxtimeGG('conj(Matrix) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('conj(real)   * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGG('conj(Scalar) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGG('conj(Vector) * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeGG('conj(Scalar) * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400));
B = rand(1,1) + rand(1,1)*1i;
maxtimeGG('conj(Array)  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeGG('conj(Vector) i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeGG('conj(Vector) o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGG('conj(Vector) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeGG('conj(Matrix) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000));
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGG('conj(Matrix) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000);
maxtimeGG('conj(Scalar) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1);
maxtimeGG('conj(Vector) * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400);
maxtimeGG('conj(Scalar) * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1);
maxtimeGG('conj(Array)  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1);
maxtimeGG('conj(Vector) i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500);
maxtimeGG('conj(Vector) o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000);
maxtimeGG('conj(Vector) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1);
maxtimeGG('conj(Matrix) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000);
maxtimeGG('conj(Matrix) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(1,1000000) + rand(1,1000000)*1i;
maxtimeGG('conj(Scalar) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000000) + rand(1,1000000)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGG('conj(Vector) * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = rand(10,20,30,400) + rand(10,20,30,400)*1i;
maxtimeGG('conj(Scalar) * conj(Array)  ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,400) + rand(10,20,30,400)*1i);
B = rand(1,1) + rand(1,1)*1i;
maxtimeGG('conj(Array)  * conj(Scalar) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = rand(10000000,1) + rand(10000000,1)*1i;
maxtimeGG('conj(Vector) i conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = rand(1,2500) + rand(1,2500)*1i;
maxtimeGG('conj(Vector) o conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2000) + rand(1,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGG('conj(Vector) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,1) + rand(2000,1)*1i;
maxtimeGG('conj(Matrix) * conj(Vector) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000,2000) + rand(2000,2000)*1i);
B = rand(2000,2000) + rand(2000,2000)*1i;
maxtimeGG('conj(Matrix) * conj(Matrix) ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs ... symmetric cases op(A) * op(A)']);
disp(' ');
disp('real');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymCN('Matrix''      * Same      ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymNC('Matrix       * Same''     ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymTN('Matrix.''     * Same      ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymNT('Matrix       * Same.''    ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymGC('conj(Matrix) * Same''     ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymCG('Matrix''      * conj(Same)',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymGT('conj(Matrix) * Same.''    ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxtimesymTG('Matrix.''     * conj(Same)',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('complex');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymCN('Matrix''      * Same      ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymNC('Matrix       * Same''     ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymTN('Matrix.''     * Same      ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymNT('Matrix       * Same.''    ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymGC('conj(Matrix) * Same''     ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymCG('Matrix''      * conj(Same)',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymGT('conj(Matrix) * Same.''    ',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxtimesymTG('Matrix.''     * conj(Same)',A,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(['Timing Tests ... median of ' num2str(n) ' runs ... special scalar cases']);
disp(' ');
disp('(scalar) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_ttable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = rand(2500);
maxtimeNN('( 1+0i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = rand(2500);
maxtimeNN('( 1+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = rand(2500);
maxtimeNN('( 1-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = rand(2500);
maxtimeNN('( 1+2i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = rand(2500);
maxtimeNN('(-1+0i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = rand(2500);
maxtimeNN('(-1+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = rand(2500);
maxtimeNN('(-1-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = rand(2500);
maxtimeNN('(-1+2i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = rand(2500);
maxtimeNN('( 2+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = rand(2500);
maxtimeNN('( 2-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(scalar) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 1+0i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 1+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 1-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 1+2i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('(-1+0i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('(-1+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('(-1-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('(-1+2i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 2+1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNN('( 2-1i) * Matrix ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(scalar) * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
%r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 1+0i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 1+1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 1-1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 1+2i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('(-1+0i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('(-1+1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('(-1-1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('(-1+2i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 2+1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = rand(2500) + rand(2500)*1i;
maxtimeNC('( 2-1i) * Matrix'' ',A,B,n,details,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(' --- DONE ! ---');
disp(' ');
disp(['Summary of Timing Tests, ' num2str(n) ' runs, + = percent faster, - = percent slower:']);
disp(' ');
mtimesx_ttable(1,1:k) = compver;
disp(mtimesx_ttable);
disp(' ');
                                                                                                                                                                                                                                                                                                            
ttable = mtimesx_ttable;
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeNN(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*B;
    mtoc(k) = toc;
    tic;
    mtimesx(A,B);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeNT(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*B.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,B,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeNC(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*B';
    mtoc(k) = toc;
    tic;
    mtimesx(A,B,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeNG(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*conj(B);
    mtoc(k) = toc;
    tic;
    mtimesx(A,B,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeTN(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*B;
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',B);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeTT(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*B.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',B,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeTC(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*B';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',B,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeTG(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*conj(B);
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',B,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeCN(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*B;
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',B);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeCT(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*B.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',B,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeCC(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*B';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',B,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeCG(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*conj(B);
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',B,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeGN(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*B;
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',B);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeGT(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*B.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',B,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeGC(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*B';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',B,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeGG(T,A,B,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*conj(B);
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',B,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
    B(1,1) = 2*B(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimeout(T,A,B,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymCN(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*A;
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',A);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymNC(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*A';
    mtoc(k) = toc;
    tic;
    mtimesx(A,A,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymTN(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*A;
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',A);
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymNT(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A*A.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,A,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymCG(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A'*conj(A);
    mtoc(k) = toc;
    tic;
    mtimesx(A,'C',A,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymGC(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*A';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',A,'C');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymTG(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    A.'*conj(A);
    mtoc(k) = toc;
    tic;
    mtimesx(A,'T',A,'G');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymGT(T,A,n,details,r)
pp(n) = 0;
mtoc(n) = 0;
xtoc(n) = 0;
for k=1:n
    tic;
    conj(A)*A.';
    mtoc(k) = toc;
    tic;
    mtimesx(A,'G',A,'T');
    xtoc(k) = toc;
    pp(k) = (100 * (xtoc(k) - mtoc(k)) / min(mtoc(k),xtoc(k)));
    A(1,1) = 2*A(1,1); % prevent JIT accelerator from interfering with timing
end
if( details )
    disp('MATLAB mtimes times:');
    disp(mtoc);
    disp('mtimesx times:')
    disp(xtoc);
    disp('mtimesx percent faster times (+ = faster, - = slower)');
    disp(-pp);
end
p = median(pp);
ap = abs(p);
sp = sprintf('%6.1f',ap);
if( ap < 5 )
    c = '(not significant)';
else
    c = '';
end
if( p < 0 )
    a = [' <' repmat('-',[1,floor((ap+5)/10)])];
    disp([T ' mtimesx is ' sp '% faster than MATLAB mtimes' a c]);
else
    disp([T ' mtimesx is ' sp '% slower than MATLAB mtimes  ' c]);
end
                                                                                                                                                                                                                                                                                                            
maxtimesymout(T,A,p,r);
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimeout(T,A,B,p,r)
global mtimesx_ttable
mtimesx_ttable(r,1:length(T)) = T;
if( isreal(A) && isreal(B) )
    lt = length(T);
    b = repmat(' ',1,30-lt);
    x = [T b sprintf('%10.0f%%',-p)];
    mtimesx_ttable(r,1:length(x)) = x;
elseif( isreal(A) && ~isreal(B) )
    x = sprintf('%10.0f%%',-p);
    mtimesx_ttable(r,42:41+length(x)) = x;
elseif( ~isreal(A) && isreal(B) )
    x = sprintf('%10.0f%%',-p);
    mtimesx_ttable(r,53:52+length(x)) = x;
else
    x = sprintf('%10.0f%%',-p);
    mtimesx_ttable(r,64:63+length(x)) = x;
end
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxtimesymout(T,A,p,r)
global mtimesx_ttable
if( isreal(A) )
    lt = length(T);
    b = repmat(' ',1,30-lt);
    x = [T b sprintf('%10.0f%%',-p)];
    mtimesx_ttable(r,1:length(x)) = x;
else
    x = sprintf('%10.0f%%',-p);
    mtimesx_ttable(r,1:length(T)) = T;
    mtimesx_ttable(r,64:63+length(x)) = x;
end
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function running_time(d)
h = 24*d;
hh = floor(h);
m = 60*(h - hh);
mm = floor(m);
s = 60*(m - mm);
ss = floor(s);
disp(' ');
rt = sprintf('Running time hh:mm:ss = %2.0f:%2.0f:%2.0f',hh,mm,ss);
if( rt(28) == ' ' )
    rt(28) = '0';
end
if( rt(31) == ' ' )
    rt(31) = '0';
end
disp(rt);
disp(' ');
return
end

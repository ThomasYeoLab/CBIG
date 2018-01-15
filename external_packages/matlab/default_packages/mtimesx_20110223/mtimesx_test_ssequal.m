% Test routine for mtimesx, op(single) * op(single) equality vs MATLAB
%******************************************************************************
%
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
%
%  Function:    mtimesx_test_ssequal
%  Filename:    mtimesx_test_ssequal.m
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
%  Syntax:
%
%    T = mtimesx_test_ssequal
%
%  Output:
%
%    T = A character array containing a summary of the results.
%
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function dtable = mtimesx_test_ssequal
                                                                                                                                                                                                                                                                                                            
global mtimesx_dtable
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('****************************************************************************');
disp('*                                                                          *');
disp('*  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  *');
disp('*                                                                          *');
disp('*  This test program can take an hour or so to complete. It is suggested   *');
disp('*  that you close all applications and run this program during your lunch  *');
disp('*  break or overnight to minimize impacts to your computer usage.          *');
disp('*                                                                          *');
disp('*  The program will be done when you see the message:  DONE !              *');
disp('*                                                                          *');
disp('****************************************************************************');
disp(' ');
try
    input('Press Enter to start test, or Ctrl-C to exit ','s');
catch
    dtable = '';
    return
end
                                                                                                                                                                                                                                                                                                            
start_time = datenum(clock);
                                                                                                                                                                                                                                                                                                            
compver = [computer ', ' version ', mtimesx mode ' mtimesx];
k = length(compver);
RC = '                                Real*Real  Real*Cplx  Cplx*Real  Cplx*Cplx';
                                                                                                                                                                                                                                                                                                            
mtimesx_dtable = char([]);
mtimesx_dtable(157,74) = ' ';
mtimesx_dtable(1,1:k) = compver;
mtimesx_dtable(2,:) = RC;
for r=3:157
mtimesx_dtable(r,:) = '                                       --         --         --         --';
end
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(compver);
disp('Test program for function mtimesx:')
disp('----------------------------------');
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
rsave = 2;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
%if( false ) % debug jump
                                                                                                                                                                                                                                                                                                            
if( isequal([]*[],mtimesx([],[])) )
    disp('Empty   * Empty                EQUAL');
else
    disp('Empty   * Empty                NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffNN('Scalar  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1));
maxdiffNN('Vector  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffNN('Scalar  * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffNN('Array   * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1));
maxdiffNN('Vector  i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500));
maxdiffNN('Vector  o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffNN('Vector  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffNN('Matrix  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffNN('Matrix  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNN('Scalar  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNN('Vector  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffNN('Scalar  * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNN('Array   * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffNN('Vector  i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffNN('Vector  o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNN('Vector  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffNN('Matrix  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNN('Matrix  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffNN('Scalar  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1));
maxdiffNN('Vector  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffNN('Scalar  * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffNN('Array   * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1));
maxdiffNN('Vector  i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500));
maxdiffNN('Vector  o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffNN('Vector  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffNN('Matrix  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffNN('Matrix  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNN('Scalar  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNN('Vector  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffNN('Scalar  * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNN('Array   * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffNN('Vector  i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffNN('Vector  o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNN('Vector  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffNN('Matrix  * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNN('Matrix  * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real) * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]*[].',mtimesx([],[],'T')) )
    disp('Empty   *  Empty.''             EQUAL');
else
    disp('Empty   *  Empty.''             NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffNT('Scalar  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffNT('Vector  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffNT('Array   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000));
maxdiffNT('Vector  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1));
maxdiffNT('Vector  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffNT('Vector  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffNT('Matrix  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffNT('Matrix  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNT('Scalar  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNT('Vector  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNT('Array   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffNT('Vector  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffNT('Vector  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNT('Vector  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffNT('Matrix  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNT('Matrix  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffNT('Scalar  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffNT('Vector  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffNT('Array   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000));
maxdiffNT('Vector  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1));
maxdiffNT('Vector  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffNT('Vector  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffNT('Matrix  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffNT('Matrix  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNT('Scalar  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNT('Vector  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNT('Array   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffNT('Vector  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffNT('Vector  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNT('Vector  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffNT('Matrix  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNT('Matrix  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real) * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]*[]',mtimesx([],[],'C')) )
    disp('Empty   *  Empty''              EQUAL');
else
    disp('Empty   *  Empty''              NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffNC('Scalar  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffNC('Vector  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffNC('Array   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000));
maxdiffNC('Vector  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1));
maxdiffNC('Vector  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffNC('Vector  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffNC('Matrix  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffNC('Matrix  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNC('Scalar  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNC('Vector  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNC('Array   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffNC('Vector  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffNC('Vector  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNC('Vector  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffNC('Matrix  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNC('Matrix  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffNC('Scalar  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffNC('Vector  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffNC('Array   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000));
maxdiffNC('Vector  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1));
maxdiffNC('Vector  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffNC('Vector  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffNC('Matrix  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffNC('Matrix  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNC('Scalar  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNC('Vector  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNC('Array   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffNC('Vector  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffNC('Vector  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNC('Vector  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffNC('Matrix  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNC('Matrix  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real) * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
%if( false ) % debug jump
                                                                                                                                                                                                                                                                                                            
if( isequal([]*conj([]),mtimesx([],[],'G')) )
    disp('Empty   * conj(Empty)          EQUAL');
else
    disp('Empty   * conj(Empty)          NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffNG('Scalar  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1));
maxdiffNG('Vector  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffNG('Scalar  * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffNG('Array   * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1));
maxdiffNG('Vector  i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500));
maxdiffNG('Vector  o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffNG('Vector  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffNG('Matrix  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffNG('Matrix  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real) * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNG('Scalar  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNG('Vector  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffNG('Scalar  * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNG('Array   * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffNG('Vector  i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffNG('Vector  o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNG('Vector  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffNG('Matrix  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNG('Matrix  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * conj((real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffNG('Scalar  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1));
maxdiffNG('Vector  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffNG('Scalar  * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffNG('Array   * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1));
maxdiffNG('Vector  i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500));
maxdiffNG('Vector  o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffNG('Vector  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffNG('Matrix  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffNG('Matrix  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex) * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffNG('Scalar  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNG('Vector  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffNG('Scalar  * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffNG('Array   * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffNG('Vector  i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffNG('Vector  o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNG('Vector  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffNG('Matrix  * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffNG('Matrix  * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real).'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*[],mtimesx([],'C',[])) )
    disp('Empty.''  * Empty               EQUAL');
else
    disp('Empty.''  * Empty               NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffTN('Scalar.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffTN('Vector.'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffTN('Scalar.'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1));
maxdiffTN('Vector.'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500));
maxdiffTN('Vector.'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffTN('Vector.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffTN('Matrix.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffTN('Matrix.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real).'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTN('Scalar.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTN('Vector.'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffTN('Scalar.'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffTN('Vector.'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffTN('Vector.'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTN('Vector.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffTN('Matrix.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTN('Matrix.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffTN('Scalar.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffTN('Vector.'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffTN('Scalar.'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1));
maxdiffTN('Vector.'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500));
maxdiffTN('Vector.'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffTN('Vector.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffTN('Matrix.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffTN('Matrix.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTN('Scalar.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTN('Vector.'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffTN('Scalar.'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffTN('Vector.'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffTN('Vector.'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTN('Vector.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffTN('Matrix.'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTN('Matrix.'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real).'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([].'*[]',mtimesx([],'T',[],'C')) )
    disp('Empty.''  *  Empty.''            EQUAL');
else
    disp('Empty.''  *  Empty.''            NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffTT('Scalar.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffTT('Vector.'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000));
maxdiffTT('Vector.'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1));
maxdiffTT('Vector.'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffTT('Vector.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffTT('Matrix.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffTT('Matrix.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real).''  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTT('Scalar.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTT('Vector.'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffTT('Vector.'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffTT('Vector.'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTT('Vector.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffTT('Matrix.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTT('Matrix.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffTT('Scalar.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffTT('Vector.'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000));
maxdiffTT('Vector.'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1));
maxdiffTT('Vector.'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffTT('Vector.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffTT('Matrix.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffTT('Matrix.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTT('Scalar.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTT('Vector.'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffTT('Vector.'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffTT('Vector.'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTT('Vector.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffTT('Matrix.'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTT('Matrix.'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real).'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([].'*[]',mtimesx([],'T',[],'C')) )
    disp('Empty.''  *  Empty''             EQUAL');
else
    disp('Empty.''  *  Empty''             NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffTC('Scalar.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffTC('Vector.'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000));
maxdiffTC('Vector.'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1));
maxdiffTC('Vector.'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffTC('Vector.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffTC('Matrix.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffTC('Matrix.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real).''  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTC('Scalar.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTC('Vector.'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffTC('Vector.'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffTC('Vector.'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTC('Vector.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffTC('Matrix.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTC('Matrix.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffTC('Scalar.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffTC('Vector.'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000));
maxdiffTC('Vector.'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1));
maxdiffTC('Vector.'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffTC('Vector.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffTC('Matrix.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffTC('Matrix.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTC('Scalar.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTC('Vector.'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffTC('Vector.'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffTC('Vector.'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTC('Vector.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffTC('Matrix.'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTC('Matrix.'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real).'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*conj([]),mtimesx([],'C',[],'G')) )
    disp('Empty.''  * conj(Empty)         EQUAL');
else
    disp('Empty.''  * conj(Empty)         NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffTG('Scalar.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffTG('Vector.'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffTG('Scalar.'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1));
maxdiffTG('Vector.'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500));
maxdiffTG('Vector.'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffTG('Vector.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffTG('Matrix.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffTG('Matrix.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real).'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTG('Scalar.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTG('Vector.'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffTG('Scalar.'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffTG('Vector.'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffTG('Vector.'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTG('Vector.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffTG('Matrix.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTG('Matrix.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffTG('Scalar.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffTG('Vector.'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffTG('Scalar.'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1));
maxdiffTG('Vector.'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500));
maxdiffTG('Vector.'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffTG('Vector.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffTG('Matrix.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffTG('Matrix.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex).'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffTG('Scalar.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffTG('Vector.'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffTG('Scalar.'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffTG('Vector.'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffTG('Vector.'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTG('Vector.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffTG('Matrix.'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffTG('Matrix.'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real)'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*[],mtimesx([],'C',[])) )
    disp('Empty''  * Empty                EQUAL');
else
    disp('Empty''  * Empty                NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffCN('Scalar'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffCN('Vector'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffCN('Scalar'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1));
maxdiffCN('Vector'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500));
maxdiffCN('Vector'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffCN('Vector'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffCN('Matrix'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffCN('Matrix'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCN('Scalar'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCN('Vector'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffCN('Scalar'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffCN('Vector'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffCN('Vector'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCN('Vector'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffCN('Matrix'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCN('Matrix'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffCN('Scalar'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffCN('Vector'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffCN('Scalar'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1));
maxdiffCN('Vector'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500));
maxdiffCN('Vector'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffCN('Vector'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffCN('Matrix'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffCN('Matrix'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCN('Scalar'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCN('Vector'' * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffCN('Scalar'' * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffCN('Vector'' i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffCN('Vector'' o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCN('Vector'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffCN('Matrix'' * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCN('Matrix'' * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real)'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*[]',mtimesx([],'C',[],'C')) )
    disp('Empty''  *  Empty.''             EQUAL');
else
    disp('Empty''  *  Empty.''             NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffCT('Scalar'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffCT('Vector'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000));
maxdiffCT('Vector'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1));
maxdiffCT('Vector'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffCT('Vector'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffCT('Matrix'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffCT('Matrix'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)''  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCT('Scalar'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCT('Vector'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffCT('Vector'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffCT('Vector'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCT('Vector'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffCT('Matrix'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCT('Matrix'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffCT('Scalar'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffCT('Vector'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000));
maxdiffCT('Vector'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1));
maxdiffCT('Vector'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffCT('Vector'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffCT('Matrix'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffCT('Matrix'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCT('Scalar'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCT('Vector'' * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffCT('Vector'' i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffCT('Vector'' o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCT('Vector'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffCT('Matrix'' * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCT('Matrix'' * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real)'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*[]',mtimesx([],'C',[],'C')) )
    disp('Empty''  *  Empty''              EQUAL');
else
    disp('Empty''  *  Empty''              NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffCC('Scalar'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffCC('Vector'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000));
maxdiffCC('Vector'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1));
maxdiffCC('Vector'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffCC('Vector'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffCC('Matrix'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffCC('Matrix'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)''  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCC('Scalar'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCC('Vector'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffCC('Vector'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffCC('Vector'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCC('Vector'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffCC('Matrix'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCC('Matrix'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffCC('Scalar'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffCC('Vector'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000));
maxdiffCC('Vector'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1));
maxdiffCC('Vector'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffCC('Vector'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffCC('Matrix'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffCC('Matrix'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCC('Scalar'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCC('Vector'' * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffCC('Vector'' i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffCC('Vector'' o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCC('Vector'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffCC('Matrix'' * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCC('Matrix'' * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('(real)'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal([]'*conj([]),mtimesx([],'C',[],'G')) )
    disp('Empty''  * conj(Empty)          EQUAL');
else
    disp('Empty''  * conj(Empty)          NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffCG('Scalar'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffCG('Vector'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffCG('Scalar'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1));
maxdiffCG('Vector'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500));
maxdiffCG('Vector'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000));
maxdiffCG('Vector'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffCG('Matrix'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffCG('Matrix'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(real)'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCG('Scalar'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCG('Vector'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffCG('Scalar'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffCG('Vector'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffCG('Vector'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCG('Vector'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffCG('Matrix'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCG('Matrix'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffCG('Scalar'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffCG('Vector'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffCG('Scalar'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1));
maxdiffCG('Vector'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500));
maxdiffCG('Vector'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000));
maxdiffCG('Vector'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffCG('Matrix'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffCG('Matrix'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('(complex)'' * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffCG('Scalar'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffCG('Vector'' * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffCG('Scalar'' * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000000,1) + rand(10000000,1)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffCG('Vector'' i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,2500) + rand(1,2500)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffCG('Vector'' o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1) + rand(1000,1)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCG('Vector'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffCG('Matrix'' * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffCG('Matrix'' * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('conj(real)   * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal(conj([])*[],mtimesx([],'G',[])) )
    disp('conj(Empty)  * Empty           EQUAL');
else
    disp('conj(Empty)  * Empty           NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffGN('conj(Scalar) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1));
maxdiffGN('conj(Vector) * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffGN('conj(Scalar) * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffGN('conj(Array)  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1));
maxdiffGN('conj(Vector) i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500));
maxdiffGN('conj(Vector) o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffGN('conj(Vector) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffGN('conj(Matrix) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffGN('conj(Matrix) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(real)  * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGN('conj(Scalar) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGN('conj(Vector) * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffGN('conj(Scalar) * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGN('conj(Array)  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffGN('conj(Vector) i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffGN('conj(Vector) o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGN('conj(Vector) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffGN('conj(Matrix) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGN('conj(Matrix) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffGN('conj(Scalar) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1));
maxdiffGN('conj(Vector) * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffGN('conj(Scalar) * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffGN('conj(Array)  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1));
maxdiffGN('conj(Vector) i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500));
maxdiffGN('conj(Vector) o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffGN('conj(Vector) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffGN('conj(Matrix) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffGN('conj(Matrix) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGN('conj(Scalar) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGN('conj(Vector) * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffGN('conj(Scalar) * Array  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGN('conj(Array)  * Scalar ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffGN('conj(Vector) i Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffGN('conj(Vector) o Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGN('conj(Vector) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffGN('conj(Matrix) * Vector ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGN('conj(Matrix) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('conj(real) * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal(conj([])*[].',mtimesx([],'G',[],'T')) )
    disp('conj(Empty)   *  Empty.''       EQUAL');
else
    disp('conj(Empty)   *  Empty.''       NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffGT('conj(Scalar)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffGT('conj(Vector)  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffGT('conj(Array)   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000));
maxdiffGT('conj(Vector)  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1));
maxdiffGT('conj(Vector)  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffGT('conj(Vector)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffGT('conj(Matrix)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffGT('conj(Matrix)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(real)  * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGT('conj(Scalar)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGT('conj(Vector)  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGT('conj(Array)   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffGT('conj(Vector)  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffGT('conj(Vector)  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGT('conj(Vector)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffGT('conj(Matrix)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGT('conj(Matrix)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex) * (real).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffGT('conj(Scalar)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffGT('conj(Vector)  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffGT('conj(Array)   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000));
maxdiffGT('conj(Vector)  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1));
maxdiffGT('conj(Vector)  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffGT('conj(Vector)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffGT('conj(Matrix)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffGT('conj(Matrix)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex) * (complex).''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGT('conj(Scalar)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGT('conj(Vector)  * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGT('conj(Array)   * Scalar.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffGT('conj(Vector)  i Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffGT('conj(Vector)  o Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGT('conj(Vector)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffGT('conj(Matrix)  * Vector.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGT('conj(Matrix)  * Matrix.'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('conj(real) * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal(conj([])*[]',mtimesx([],'G',[],'C')) )
    disp('conj(Empty)   *  Empty''        EQUAL');
else
    disp('conj(Empty)   *  Empty''        NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10000,1));
maxdiffGC('conj(Scalar)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1));
maxdiffGC('conj(Vector)  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffGC('conj(Array)   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000));
maxdiffGC('conj(Vector)  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1));
maxdiffGC('conj(Vector)  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffGC('conj(Vector)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000));
maxdiffGC('conj(Matrix)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffGC('conj(Matrix)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(real)  * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGC('conj(Scalar)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGC('conj(Vector)  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGC('conj(Array)   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffGC('conj(Vector)  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffGC('conj(Vector)  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGC('conj(Vector)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffGC('conj(Matrix)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGC('conj(Matrix)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex) * (real)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffGC('conj(Scalar)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1));
maxdiffGC('conj(Vector)  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffGC('conj(Array)   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000));
maxdiffGC('conj(Vector)  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1));
maxdiffGC('conj(Vector)  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffGC('conj(Vector)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000));
maxdiffGC('conj(Matrix)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffGC('conj(Matrix)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex) * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGC('conj(Scalar)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10000,1)+ rand(10000,1)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGC('conj(Vector)  * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGC('conj(Array)   * Scalar'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(1,10000000) + rand(1,10000000)*1i);
maxdiffGC('conj(Vector)  i Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(2500,1) + rand(2500,1)*1i);
maxdiffGC('conj(Vector)  o Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGC('conj(Vector)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1,1000) + rand(1,1000)*1i);
maxdiffGC('conj(Matrix)  * Vector'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGC('conj(Matrix)  * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('Numerical Comparison Tests ...');
disp(' ');
disp('conj(real)   * conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
if( isequal(conj([])*conj([]),mtimesx([],'G',[],'G')) )
    disp('conj(Empty)  * conj(Empty)     EQUAL');
else
    disp('conj(Empty)  * conj(Empty)     NOT EQUAL <---');
end
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000));
maxdiffGG('conj(Scalar) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1));
maxdiffGG('conj(Vector) * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40));
maxdiffGG('conj(Scalar) * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1));
maxdiffGG('conj(Array)  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1));
maxdiffGG('conj(Vector) i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500));
maxdiffGG('conj(Vector) o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000));
maxdiffGG('conj(Vector) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1));
maxdiffGG('conj(Matrix) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000));
maxdiffGG('conj(Matrix) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(real)  * conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGG('conj(Scalar) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGG('conj(Vector) * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1));
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffGG('conj(Scalar) * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40));
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGG('conj(Array)  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000));
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffGG('conj(Vector) i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1));
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffGG('conj(Vector) o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGG('conj(Vector) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffGG('conj(Matrix) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000));
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGG('conj(Matrix) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* conj(real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000));
maxdiffGG('conj(Scalar) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1));
maxdiffGG('conj(Vector) * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40));
maxdiffGG('conj(Scalar) * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1));
maxdiffGG('conj(Array)  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1));
maxdiffGG('conj(Vector) i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500));
maxdiffGG('conj(Vector) o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000));
maxdiffGG('conj(Vector) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1));
maxdiffGG('conj(Matrix) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000));
maxdiffGG('conj(Matrix) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
disp(' ');
disp('conj(complex)* conj(complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(1,10000) + rand(1,10000)*1i);
maxdiffGG('conj(Scalar) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000)+ rand(1,10000)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGG('conj(Vector) * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1) + rand(1,1)*1i);
B = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
maxdiffGG('conj(Scalar) * conj(Array)  ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(10,20,30,40) + rand(10,20,30,40)*1i);
B = single(rand(1,1) + rand(1,1)*1i);
maxdiffGG('conj(Array)  * conj(Scalar) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,10000000) + rand(1,10000000)*1i);
B = single(rand(10000000,1) + rand(10000000,1)*1i);
maxdiffGG('conj(Vector) i conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2500,1) + rand(2500,1)*1i);
B = single(rand(1,2500) + rand(1,2500)*1i);
maxdiffGG('conj(Vector) o conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1,1000) + rand(1,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGG('conj(Vector) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1) + rand(1000,1)*1i);
maxdiffGG('conj(Matrix) * conj(Vector) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(1000,1000) + rand(1000,1000)*1i);
B = single(rand(1000,1000) + rand(1000,1000)*1i);
maxdiffGG('conj(Matrix) * conj(Matrix) ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp('----------------------------------');
disp(' ');
disp('Numerical Comparison Tests ... symmetric cases op(A) * op(A)');
disp(' ');
disp('real');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = RC;
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymCN('Matrix''      * Same ',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymNC('Matrix       * Same''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymTN('Matrix.''     * Same ',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymNT('Matrix       * Same.''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymGC('conj(Matrix) * Same''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymCG('Matrix''      * conj(Same)',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymGT('conj(Matrix) * Same.'' ',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000));
maxdiffsymTG('Matrix.''     * conj(Same)',A,r);
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('complex');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymCN('Matrix''      * Same ',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymNC('Matrix       * Same''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymTN('Matrix.''     * Same ',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymNT('Matrix       * Same.''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymGC('conj(Matrix) * Same''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymCG('Matrix''      * conj(Same)',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymGT('conj(Matrix) * Same.''',A,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(rand(2000) + rand(2000)*1i);
maxdiffsymTG('Matrix.''     * conj(Same)',A,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
%end % debug jump
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('Numerical Comparison Tests ... special scalar cases');
disp(' ');
disp('(scalar) * (real)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = r + 1;
mtimesx_dtable(r,:) = '                                Real*Real  Real*Cplx  Cplx*Real  Cplx*Cplx';
                                                                                                                                                                                                                                                                                                            
rsave = r;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = single(rand(2500));
maxdiffNN('( 1+0i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = single(rand(2500));
maxdiffNN('( 1+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = single(rand(2500));
maxdiffNN('( 1-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = single(rand(2500));
maxdiffNN('( 1+2i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = single(rand(2500));
maxdiffNN('(-1+0i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = single(rand(2500));
maxdiffNN('(-1+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = single(rand(2500));
maxdiffNN('(-1-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = single(rand(2500));
maxdiffNN('(-1+2i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = single(rand(2500));
maxdiffNN('( 2+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = single(rand(2500));
maxdiffNN('( 2-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(scalar) * (complex)');
disp(' ');
                                                                                                                                                                                                                                                                                                            
r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 1+0i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 1+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 1-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 1+2i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('(-1+0i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('(-1+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('(-1-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('(-1+2i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 2+1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNN('( 2-1i) * Matrix ',A,B,r);
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp('(scalar) * (complex)''');
disp(' ');
                                                                                                                                                                                                                                                                                                            
%r = rsave;
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 1+0i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 1+1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 1-1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(1 + 2i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 1+2i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('(-1+0i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('(-1+1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('(-1-1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(-1 + 2i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('(-1+2i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 + 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 2+1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
r = r + 1;
A = single(2 - 1i);
B = single(rand(2500) + rand(2500)*1i);
maxdiffNC('( 2-1i) * Matrix'' ',A,B,r);
                                                                                                                                                                                                                                                                                                            
running_time(datenum(clock) - start_time);
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
disp(' ');
disp(' --- DONE ! ---');
disp(' ');
disp('Summary of Numerical Comparison Tests, max relative element difference:');
disp(' ');
mtimesx_dtable(1,1:k) = compver;
disp(mtimesx_dtable);
disp(' ');
                                                                                                                                                                                                                                                                                                            
dtable = mtimesx_dtable;
                                                                                                                                                                                                                                                                                                            
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffNN(T,A,B,r)
Cm = A*B;
Cx = mtimesx(A,B);
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffCN(T,A,B,r)
Cm = A'*B;
Cx = mtimesx(A,'C',B);
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffTN(T,A,B,r)
Cm = A.'*B;
Cx = mtimesx(A,'T',B);
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffGN(T,A,B,r)
Cm = conj(A)*B;
Cx = mtimesx(A,'G',B);
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffNC(T,A,B,r)
Cm = A*B';
Cx = mtimesx(A,B,'C');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffCC(T,A,B,r)
Cm = A'*B';
Cx = mtimesx(A,'C',B,'C');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffTC(T,A,B,r)
Cm = A.'*B';
Cx = mtimesx(A,'T',B,'C');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffGC(T,A,B,r)
Cm = conj(A)*B';
Cx = mtimesx(A,'G',B,'C');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffNT(T,A,B,r)
Cm = A*B.';
Cx = mtimesx(A,B,'T');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffCT(T,A,B,r)
Cm = A'*B.';
Cx = mtimesx(A,'C',B,'T');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffTT(T,A,B,r)
Cm = A.'*B.';
Cx = mtimesx(A,'T',B,'T');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffGT(T,A,B,r)
Cm = conj(A)*B.';
Cx = mtimesx(A,'G',B,'T');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffNG(T,A,B,r)
Cm = A*conj(B);
Cx = mtimesx(A,B,'G');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffCG(T,A,B,r)
Cm = A'*conj(B);
Cx = mtimesx(A,'C',B,'G');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffTG(T,A,B,r)
Cm = A.'*conj(B);
Cx = mtimesx(A,'T',B,'G');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffGG(T,A,B,r)
Cm = conj(A)*conj(B);
Cx = mtimesx(A,'G',B,'G');
maxdiffout(T,A,B,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymCN(T,A,r)
Cm = A'*A;
Cx = mtimesx(A,'C',A);
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymNC(T,A,r)
Cm = A*A';
Cx = mtimesx(A,A,'C');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymTN(T,A,r)
Cm = A.'*A;
Cx = mtimesx(A,'T',A);
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymNT(T,A,r)
Cm = A*A.';
Cx = mtimesx(A,A,'T');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymTG(T,A,r)
Cm = A.'*conj(A);
Cx = mtimesx(A,'T',A,'G');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymGT(T,A,r)
Cm = conj(A)*A.';
Cx = mtimesx(A,'G',A,'T');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymCG(T,A,r)
Cm = A'*conj(A);
Cx = mtimesx(A,'C',A,'G');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymGC(T,A,r)
Cm = conj(A)*A';
Cx = mtimesx(A,'G',A,'C');
maxdiffsymout(T,A,Cm,Cx,r);
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffout(T,A,B,Cm,Cx,r)
global mtimesx_dtable
lt = length(T);
b = repmat(' ',1,30-lt);
if( isequal(Cm,Cx) )
    disp([T b ' EQUAL']);
    d = 0;
else
    Cm = Cm(:);
    Cx = Cx(:);
    if( isreal(Cm) && isreal(Cx) )
        rx = Cx ~= Cm;
        d = max(abs((Cx(rx)-Cm(rx))./Cm(rx)));
    else
        Cmr = real(Cm);
        Cmi = imag(Cm);
        Cxr = real(Cx);
        Cxi = imag(Cx);
        rx = Cxr ~= Cmr;
        ix = Cxi ~= Cmi;
        dr = max(abs((Cxr(rx)-Cmr(rx))./max(abs(Cmr(rx)),abs(Cmr(rx)))));
        di = max(abs((Cxi(ix)-Cmi(ix))./max(abs(Cmi(ix)),abs(Cxi(ix)))));
        if( isempty(dr) )
            d = di;
        elseif( isempty(di) )
            d = dr;
        else
            d = max(dr,di);
        end
    end
    disp([T b ' NOT EQUAL <--- Max relative difference: ' num2str(d)]);
end
mtimesx_dtable(r,1:length(T)) = T;
if( isreal(A) && isreal(B) )
    if( d == 0 )
        x = [T b '          0'];
    else
        x = [T b sprintf('%11.2e',d)];
    end
    mtimesx_dtable(r,1:length(x)) = x;
elseif( isreal(A) && ~isreal(B) )
    if( d == 0 )
        x = '          0';
    else
        x = sprintf('%11.2e',d);
    end
    mtimesx_dtable(r,42:41+length(x)) = x;
elseif( ~isreal(A) && isreal(B) )
    if( d == 0 )
        x = '          0';
    else
        x = sprintf('%11.2e',d);
    end
    mtimesx_dtable(r,53:52+length(x)) = x;
else
    if( d == 0 )
        x = '          0';
    else
        x = sprintf('%11.2e',d);
    end
    mtimesx_dtable(r,64:63+length(x)) = x;
end
                                                                                                                                                                                                                                                                                                            
return
end
                                                                                                                                                                                                                                                                                                            
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
                                                                                                                                                                                                                                                                                                            
function maxdiffsymout(T,A,Cm,Cx,r)
global mtimesx_dtable
lt = length(T);
b = repmat(' ',1,30-lt);
if( isequal(Cm,Cx) )
    disp([T b ' EQUAL']);
    d = 0;
else
    Cm = Cm(:);
    Cx = Cx(:);
    if( isreal(Cm) && isreal(Cx) )
        rx = Cx ~= Cm;
        d = max(abs((Cx(rx)-Cm(rx))./Cm(rx)));
    else
        Cmr = real(Cm);
        Cmi = imag(Cm);
        Cxr = real(Cx);
        Cxi = imag(Cx);
        rx = Cxr ~= Cmr;
        ix = Cxi ~= Cmi;
        dr = max(abs((Cxr(rx)-Cmr(rx))./max(abs(Cmr(rx)),abs(Cmr(rx)))));
        di = max(abs((Cxi(ix)-Cmi(ix))./max(abs(Cmi(ix)),abs(Cxi(ix)))));
        if( isempty(dr) )
            d = di;
        elseif( isempty(di) )
            d = dr;
        else
            d = max(dr,di);
        end
    end
    disp([T b ' NOT EQUAL <--- Max relative difference: ' num2str(d)]);
end
if( isreal(A) )
    if( d == 0 )
        x = [T b '          0'];
    else
        x = [T b sprintf('%11.2e',d)];
    end
    mtimesx_dtable(r,1:length(x)) = x;
else
    if( d == 0 )
        x = '          0';
    else
        x = sprintf('%11.2e',d);
    end
    mtimesx_dtable(r,1:length(T)) = T;
    mtimesx_dtable(r,64:63+length(x)) = x;
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

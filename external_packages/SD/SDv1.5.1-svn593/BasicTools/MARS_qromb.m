% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================
function value = MARS_qromb(functor, a, b, varargin)

EPS = 1.0e-3;
DIFF = 1.0e-3;
JMAX = 25;
JMAXP = 25+1;
K = 3;

if(a == b)
   value = 0;
   return;
end

if(isempty(varargin))

    s = zeros(JMAXP, 1);
    h = zeros(JMAXP+1, 1);

    h(1) = 1;
    for j=1:JMAX

        if(j==1)
            s(j) = MARS_dumbtrapzd(functor, a, b, j, 0);
        else
            dummy = MARS_dumbtrapzd(functor, a, b, j, s(j-1));
            s(j) = dummy;
        end

        if(j >= K)
            [ss, dss] = MARS_polint(h(j-K+1:j), s(j-K+1:j), K, 0);
            if (abs(dss)<= EPS*abs(ss) || abs(dss)<EPS*DIFF)
                value = ss;
                return;
            end
        end

        h(j+1) = 0.25*h(j);

    end
    error('Too many steps in qromb');
    

else

    s = zeros(JMAXP, 1);
    h = zeros(JMAXP+1, 1);

    h(1) = 1;
    for j=1:JMAX

        if(j==1)
            s(j) = MARS_dumbtrapzd(functor, a, b, j, 0, varargin{1:end});
        else
            dummy = MARS_dumbtrapzd(functor, a, b, j, s(j-1), varargin{1:end});
            s(j) = dummy;
        end

        if(j >= K)
            [ss, dss] = MARS_polint(h(j-K+1:j), s(j-K+1:j), K, 0);
            if (abs(dss)<= EPS*abs(ss) || abs(dss)<EPS*DIFF)
                value = ss;
                return;
            end
        end

        h(j+1) = 0.25*h(j);

    end
    error('Too many steps in qromb');


end
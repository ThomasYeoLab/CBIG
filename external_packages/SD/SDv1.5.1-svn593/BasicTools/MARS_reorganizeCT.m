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
function MARS_ct = MARS_reorganizeCT(ct, structures_of_interest)

% MARS_ct = MARS_reorganizeCT(ct, structures_of_interest)


if(~isempty(structures_of_interest))
    
    MARS_ct.orig_tab = ct.orig_tab;
    MARS_ct.numEntries = length(structures_of_interest);
    
    MARS_ct.struct_names = structures_of_interest;
    MARS_ct.table = zeros(MARS_ct.numEntries, size(ct.table, 2));
    
    
    for i = 1:MARS_ct.numEntries
        index = find(strcmp(ct.struct_names, structures_of_interest{i}) == 1);
        
        if(length(index) ~= 1)
           error('More than one structures have same names!!'); 
        end
        
        MARS_ct.table(i,:) = ct.table(index, :);    
    end
else
    count = 0;
    MARS_ct.orig_tab = ct.orig_tab;

    for i = 1:ct.numEntries
        if(~isempty(ct.struct_names{i}))
            count = count + 1;
        end
    end

    MARS_ct.numEntries = count;
    MARS_ct.struct_names = cell(count, 1);
    MARS_ct.table = zeros(count, size(ct.table, 2));

    count = 0;
    for i = 1:ct.numEntries
        if(~isempty(ct.struct_names{i}))
            count = count + 1;
            MARS_ct.struct_names{count} = ct.struct_names{i};
            MARS_ct.table(count, :) = ct.table(i,:);
        end
    end
end
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
function [MARS_label, MARS_ct, index] = MARS_reorganizeLabels(label, ct, vertNbors)

MARS_ct = ct;
MARS_ct.table = [MARS_ct.table (1:ct.numEntries)'];


%first check for labels that are not defined within the ct!!
temp = label;
index = find(temp == -1);
if(isempty(index))
    for i = 1:ct.numEntries
        temp(label == ct.table(i,5)) = -1;
    end
    index = find(temp~=-1);
    
    if(~isempty(index))
        warning(['There exists ' num2str(length(index)) ' points not found in the colortable!']); 
        disp(['index is ' num2str(index')]); 
        disp(['value is ' num2str(label(index)')]);
        
        if(length(index) > 15)
           error('More than 5 unlabeled points. Manual Intervention needed!');
        else
           disp('Filling with neighbors...');
        end
        
        count = 0;
        for i = 1:length(index)
            vno = index(i);
            j = 1;
            while(vertNbors(j, vno) ~= 0)
                
                v_neighbor = vertNbors(j, vno);
                if(isempty(find(index == v_neighbor)))
                    %neighbor is labeled!
                    label(vno) = label(v_neighbor);
                    count = count + 1;
                    break;
                    
                end
                j = j + 1;
            end
        end
        if(count~=length(index))
            error('Not all the unlabeled points are filled!!');
        else
            disp('All labels filled');
        end
    end
    
else
   error('There exists a point with label -1!!'); 
end




MARS_label = label;
for i = 1:ct.numEntries
    MARS_label(MARS_label == ct.table(i,5)) = i;
end

MARS_label = int32(MARS_label);
MARS_ct.numEntries = int32(MARS_ct.numEntries);
MARS_ct.table = int32(MARS_ct.table);



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
function [unfold_vertices, folded_vertices_list] = MARS_unfoldMesh(MARS_atlas, sbjWarp, folding_energy, STEP_SIZE, MAX_STEPS, lambda, folded_vertices_list, MAX_UNFOLD_ITER)

% unfold_vertices = MARS_unfoldMesh(MARS_atlas, sbjWarp, folding_energy,
% STEP_SIZE, MAX_STEPS, lambda)

if(nargin < 8)
    MAX_UNFOLD_ITER = 1000;
end
temp_sbjWarp = sbjWarp;

disp('********** UNFOLD *************');

num_step_memory = 5;
step_mem_vec = zeros(1, num_step_memory);
step_adaptation = 1.5;
MAX_STEP_SIZE_FACTOR = 1.5;
MIN_STEP_SIZE_FACTOR = 0.1;
STEP_SIZE_FACTOR = STEP_SIZE/sqrt(min(MARS_atlas.vertexDistSq2Nbors(1,:)));

%time1 = toc;
disp(['j: 0, Length of folded list is : ' num2str(length(folded_vertices_list))]);
prev_folded_list = folded_vertices_list;
for j = 1:MAX_UNFOLD_ITER
   
    %First compute gradient of folds. 
    [energy, grad, folded_vertices_list] = MARS_computeFoldingGradFast(sbjWarp.curr_vertices, MARS_atlas.MARS_uniformMesh_SpatialPrior, folded_vertices_list);    
      
    grad = grad./repmat(sqrt(sum(grad.^2,1) + eps), [3 1]);    
    
    if(j > 1000)
        if(isempty(setxor(prev_folded_list, folded_vertices_list)))
            %might be stuck in an infinite loop
            grad(:, folded_vertices_list) = grad(:, folded_vertices_list) + randn(3, length(folded_vertices_list))*0.1;        
            grad = grad./repmat(sqrt(sum(grad.^2,1) + eps), [3 1]); 
        end
    end
    prev_folded_list = folded_vertices_list;
    
    grad = -lambda * grad;
    max_grad = sqrt(max(sum(grad.^2, 1)));       
    if(max_grad == 0)
       unfold_vertices = temp_sbjWarp.curr_vertices;
       break;
    end 
    
    
%     if(length(folded_vertices_list) > 1000 && j == 1)
%         length(folded_vertices_list)
%         keyboard;
%     end
    %Perform line search
    for i = 1:MAX_STEPS 
    
        curr_step_size = 2^(-i+1)*STEP_SIZE / max_grad;
        temp_sbjWarp.curr_vertices(:, folded_vertices_list) = MARS_warpPointbyTangentVec(sbjWarp.curr_vertices(:,folded_vertices_list), curr_step_size*grad(:,folded_vertices_list));

        [new_folding_energy, folded_vertices_list] = MARS_computeFoldingEnergyFast(temp_sbjWarp.curr_vertices, MARS_atlas.MARS_uniformMesh_SpatialPrior, folded_vertices_list);       
        
        new_folding_energy = -lambda * new_folding_energy;
               
        
        if(new_folding_energy > folding_energy)
            %disp(['Unfold iter. ' num2str(j) ', Line search iter. ' num2str(i) ': curr_step: ' num2str(curr_step_size) ', old_ener: ' num2str(folding_energy) ' < new_energy: ' num2str(new_folding_energy)]);
            break;
        else
            %disp(['Unfold iter. ' num2str(j) ', Line search iter. ' num2str(i) ': curr_step: ' num2str(curr_step_size) ', old_ener: ' num2str(folding_energy) ' > new_energy: ' num2str(new_folding_energy)]);
        end
    end
    
    %%%%%diagnostics %%%%%%%%%%
    if(find(j == [1:9 10:10:90 1e2:1e2:9e2 1e3:1e3:10e4]))
        disp(['j: ' num2str(j) ', Length of folded list is : ' num2str(length(folded_vertices_list))]);
    end
    
    
    step_mem_vec(mod(j, num_step_memory)+1) = 1/i*STEP_SIZE_FACTOR; %fill in the step_mem_vec in a round robin manner
    if(j >= num_step_memory) %When there's enough memory change adaptive step.
        avg_step = mean(step_mem_vec); %average step taken in the past num_step_memory steps
        STEP_SIZE_FACTOR = max(min(avg_step * step_adaptation, MAX_STEP_SIZE_FACTOR), MIN_STEP_SIZE_FACTOR); % change it accordingly...
        STEP_SIZE = STEP_SIZE_FACTOR * sqrt(min(MARS_atlas.vertexDistSq2Nbors(1,:)));
        %disp(['STEP_TAKEN: ' num2str(i) ', STEP_SIZE_FACTOR: ' num2str(STEP_SIZE_FACTOR) ', STEP_SIZE:' num2str(STEP_SIZE)]);
    end
            
    %Whether line search succeeds or not, change the folding energy and
    %vertices
    folding_energy = new_folding_energy;
    sbjWarp.curr_vertices = temp_sbjWarp.curr_vertices;
    
    %If unfolded...
    if(folding_energy == 0.0)
        break;
    end
    
end

unfold_vertices = temp_sbjWarp.curr_vertices;
%time2 = toc;

disp(['j: ' num2str(j) ', Length of folded list is : ' num2str(length(folded_vertices_list))]);
% if(folding_energy < 0)
%    keyboard;
%    error('Unfolding failed!!'); 
% end

%disp(['Time Taken: ' num2str(time2-time1) ', takes iteration: ' num2str(j)]);
disp('******** END UNFOLD ***********');






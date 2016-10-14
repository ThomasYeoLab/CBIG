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
function plot_ty_boundary_dists_diff(left_hemi_file1, right_hemi_file1, left_hemi_file2, right_hemi_file2)

% plot_ty_boundary_dists_diff('DistMat.sphere.mat', 'DistMat.sphere.rh.mat', 'lh.DistMat.sphere.DT.TLOO2.DIRECT.0.mat ', 'rh.DistMat.sphere.DT.TLOO2.DIRECT.0.mat')

dm_lh1 = ComputeStatisticsOnDistMat(left_hemi_file1) ;
dm_rh1 = ComputeStatisticsOnDistMat(right_hemi_file1) ;

dm_lh2 = ComputeStatisticsOnDistMat(left_hemi_file2) ;
dm_rh2 = ComputeStatisticsOnDistMat(right_hemi_file2) ;

dm_lh = dm_lh1 - dm_lh2;
dm_rh = dm_rh1 - dm_rh2;

means_lh = dm_lh([1:4 8],:);
for i=1:3
  tmp = dm_lh(i+4,:);
  ind = find(tmp ~= 0);
  means46_lh(i,:) = tmp(ind) ;
  clear('tmp') ;
end

means_rh = dm_rh([1:4 8],:);
for i=1:3
  tmp = dm_rh(i+4,:);
  ind = find(tmp ~= 0);
  means46_rh(i,:) = tmp(ind) ;
  clear('tmp') ;
end

means_lh = means_lh' ;
means_rh = means_rh' ;
means46_lh = means46_lh' ;
means46_rh = means46_rh' ;


%h = boxplot(squeeze(means(:,1,:))', 'labels',mlabels,'whisker',6);
%hold on ;
%h2 = boxplot(squeeze(means46(:,1,:))', 'positions', nlabels+1:nlabels+2, ...
%             'labels', ['BA4a', 'BA4p', 'BA6']);
xt = str2mat('V1', 'V2', 'BA44', 'BA45', 'BA2', 'BA4a', 'BA4p', 'BA6');
% v1 v2 44 45 2 4a 4p 6  --> v1 4a 4p 2 v2 6 44 45
order = [1 6 7 5 2 8 3 4] ;

sensory = [1 2 5 6 7 ] ;
primary_sensory = [1  5 6 7 ] ;
higher = [3 4] ;

for h=1:2

  if (h == 1)
    hemi = 'left' ;
    d = boxprep(means_lh, means46_lh) ;
  else
    hemi = 'right' ;
    d = boxprep(means_rh, means46_rh) ;
  end
  for l=1:size(d,2)
    if (l <=5)
      med(h,l) = median(d(:,l)) ;
    else
      med(h,l) = median(d(1:8,l)) ;
    end
  end
  figure('name', sprintf('%s hemi BAs', hemi)) ;
  
  hndl = boxplot(d(:,order), 'labels', xt(order,:), 'notch', 'on') ;
  fs = 14;
  set(gca, 'fontsize', fs, 'fontweight', 'bold') ;
  xlabel('Brodmann Area', 'fontsize', fs, 'fontweight', 'bold') ;
  ylabel('average mm between boundaries (diff)', 'fontsize', fs, 'fontweight', 'bold');
  hold on ;
  for i=1:length(order)
    plot(i, d(:,order(i)), 'go', 'linewidth', 2) ;
  end
  hold off ;
  for i=1:length(hndl(:))
    if (isnan(hndl(i)) == 0)
      set(hndl(i), 'markersize', 10, 'linewidth', 2)
    end
  end
  set(gca, 'ylim', [-2 4.5]) ;
  set(gca, 'YTick', -2:2:4);
  title(sprintf('predictability of BAs (%s hemi)',hemi),'fontsize', fs, 'fontweight', 'bold');
end


if 0
smean = zeros(nlabels+3, 2) ;
smean(1:nlabels,:) = mean(sigmas,3) ;
smean(nlabels+1:nlabels+2,:) = mean(sigmas4,3) ;
smean(nlabels+3,:) = mean(sigmas6,3) ;

sstd(1:nlabels,:) = std(sigmas,0, 3) / sqrt(nsubjects) ;
sstd(nlabels+1:nlabels+2,:) = std(sigmas4,0,3) / sqrt(nsubjects4) ;
sstd(nlabels+3,:) = std(sigmas6,0,3) / sqrt(nsubjects6) ;


order = [1 7 2 6 3 8 4 5];
bh = bar(smean(order,:)) ;
xt = [mlabels;'BA4a'; 'BA4p'; 'BA6 '];
set(gca, 'xticklabel', xt(order,:)) ;
fs = 14 ;
set(gca, 'fontsize', fs, 'fontweight', 'bold') ;
ylabel('mm blurring (\sigma)', 'fontsize', fs, 'fontweight', 'bold') ;
xlabel('Brodmann Area', 'fontsize', fs, 'fontweight', 'bold') ;
set(gca, 'xlim', [.5 nlabels+3.5]) ;

bw = .15 ;
for l=1:size(smean,1)
  ind = order(l) ;
  ln = line([l-bw l-bw], [smean(ind,1) smean(ind,1)+sstd(ind,1)]) ;
  set(ln, 'linewidth', 4, 'color', [0 0 0]);
  ln = line([l+bw l+bw], [smean(ind,2) smean(ind,2)+sstd(ind,2)]) ;
  set(ln, 'linewidth', 4, 'color', [0 0 0]);
end
title('Variability of Brodmann areas (blue=lh,red=rh)','fontsize', fs, 'fontweight', 'bold') ;
end
mean_sensory = mean(med(:,sensory)')
mean_higher = mean(med(:,higher)')

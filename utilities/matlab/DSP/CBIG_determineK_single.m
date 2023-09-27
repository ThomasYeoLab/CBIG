function stabilityvals = CBIG_determineK_single(x,res)

% stabilityvals = CBIG_determineK_single(x,res)
%
% Compute stability of parcellation computed from re-sampled data. The
% stability is defined by the cost of Hungarian matching between two label
% assignments. The first label assignment is computing the probabilistic
% cluster assignments on the first half using the clustering parameters
% from the SAME FIRST half. The second one is computing the probabilistic
% cluster assignments on the first half using the clustering parameters
% from the second half.
%
% x:   input correlation profiles data
% res: clustering results based on re-sampled data
%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


B = size(res,2);
Kmax = size(res,1); 


for k = 2:Kmax
    
    for i=1:B
        
        %[matching cost2] = Hungarian(CBIG_dist2(res{k,i,1}.mtc,res{k,i,2}.mtc));
        
        % Computing the probabilistic cluster assignments on the first half using 
        % the clustering parameters from the second half

        distctrs1 = res{k,i,2}.lambda*x(res{k,i,1}.inds,:)*(res{k,i,2}.mtc');
        maxdists1 = max(distctrs1,[],2);
        
        rr1 = exp(distctrs1 - maxdists1*ones(1,k));
        sumrr1 = sum(rr1,2);
        prob2 = rr1 ./ sumrr1(:,ones(k,1));

        % Computing the probabilistic cluster assignments on the first half using 
        % the clustering parameters from the SAME FIRST half
                
        distctrs1 = res{k,i,1}.lambda*x(res{k,i,1}.inds,:)*(res{k,i,1}.mtc');
        maxdists1 = max(distctrs1,[],2);
        
        r1 = exp(distctrs1 - maxdists1*ones(1,k));
        sumr1 = sum(r1,2);
        prob1 = r1 ./ sumr1(:,ones(k,1));

        % Turning the probabilistic clustering assignments into hard cluster memberships

        maxr1 = max(prob1,[],2);
        assign1 = ~(prob1-maxr1*ones(1,size(prob1,2)));

        maxr2 = max(prob2,[],2);
        assign2 = ~(prob2-maxr2*ones(1,size(prob2,2)));

        % Do the matching 

        [matching cost1] = Hungarian(-assign1*assign2');

        % Creating random cluster assignments

        r1r = rand(size(r1));
        rr1r=rand(size(rr1));
        
        maxr1r = max(r1r,[],2);
        arrign1_random = ~(r1r-maxr1r*ones(1,size(r1,2)));
        
        maxrr1r = max(rr1r,[],2);
        arrign2_random = ~(rr1r-maxrr1r*ones(1,size(r1,2)));

        [matching cost2] = Hungarian(-assign1_random*assign2_random');

        % Computing stability values
        
        %stab(k,i) = mean(sum(assign1.*assign2,2));
        %rstab(k,i)= mean(sum(arrign1_random.*arrign2_random,2));
        
        stab(k,i) = 1 + cost1 ./ size(assign1,1);
        rstab(k,i) = 1 + cost2 ./ size(assign1,1);

    end

  
end

stabilityvals{1}=stab;
stabilityvals{2}=rstab;

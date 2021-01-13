function labels = watershed_algorithm_all_par_cifti(edgemetrics, ...
    minimametrics, stepnum, fracmaxh, neighbors, minh, maxh)

%Label initial markers with unique value
labels = zeros(size(minimametrics), 'single');

if ~exist('minh', 'var') || isempty(minh)
    minh = min(edgemetrics(:));
end
if ~exist('maxh', 'var') || isempty(maxh)
    maxh = max(edgemetrics(:));
end

stoph = maxh*fracmaxh;
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;

%system('mkdir countingDir')
%system('touch countingDir/countfile')
%poolname = parpool(8);

numlabels = size(labels,2);
divisions = 4;
labelsperdivision = floor(numlabels / divisions);
for j = 1:divisions
    
    if j <= divisions
        labelstodo = ((j-1)*labelsperdivision + 1):(j*labelsperdivision);
    else
        labelstodo = ((j-1)*labelsperdivision + 1):numlabels;
    end
    %parfor l = labelstodo%size(labels,2)
    for l = labelstodo%size(labels,2)
        %system(['echo ' num2str(l) ' >> countingDir/countfile']);
        %[ign count] = system('cat countingDir/countfile | wc -l');
        %disp(['Calculating edges on metric #:' num2str(l) ', Count = ' num2str(count)])
        label = labels(:,l);
        edgemetric = edgemetrics(:,l);
        
        labelpos = find(minimametrics(:,l)==1);
        randval = randn(length(labelpos),1);
        [~, randind] = sort(randval);
        temp = 1:length(labelpos);
        labelnums = temp(randind);
        label(labelpos) = labelnums;
        watershed_zones = zeros(size(label));
        
        for i = 1:length(hiter)
            %disp(['Number of iterations will be ' num2str(length(hiter)) ', Iteration = ' num2str(i)])
            %maskpos = find(edgemetric<sortedge(i)); % Take values in metric less than current iteration
            maskmetrics = edgemetric<hiter(i); % Take values in metric less than current iteration
            maskmetrics = maskmetrics & ~label>0 & ~watershed_zones;
            
            maskpos = find(sum(maskmetrics,2)>0);
            randnum = randperm(length(maskpos));
            maskpos = maskpos(randnum);
            
            for m = 1:length(maskpos) %For all nodes at this threshold
                nodeneigh = neighbors(maskpos(m),2:end);
                maskinthismetric = maskmetrics(maskpos(m),:);
                %nodeneigh = neighbors(sortedgepos(i),2:7);
                
                nodeneigh(isnan(nodeneigh)) = [];
                nodeneighlab = label(nodeneigh,:);
                
                %Find minimum value other than 0 among neighbors
                minfindnodeneighlab = nodeneighlab;
                minfindnodeneighlab(nodeneighlab==0) = 100000;
                minnodeneighlab = min(minfindnodeneighlab,[],1);
                
                %Find maximum value other than 0 among neighbors
                maxfindnodeneighlab = nodeneighlab;
                maxfindnodeneighlab(nodeneighlab==0) = -100000;
                maxnodeneighlab = max(maxfindnodeneighlab,[],1);
                
                %If min and max differ (i.e. two or more neighbor water), watershed
                %zone
                watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                watershed_zones(maskpos(m),watershed_nodes) = 1;
                label(maskpos(m),watershed_nodes) = 0;
                
                %If min and max the same but different from 0, add to neighbor
                %water
                next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water); 
            end
        end
        labels(:,l) = label;
    end
end
%system('rm -r countingDir')
%delete(parpool)
%delete(poolname);
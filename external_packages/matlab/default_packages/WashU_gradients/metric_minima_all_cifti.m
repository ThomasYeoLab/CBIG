function minimametric = metric_minima_all_cifti(metric,neighdist,neighbors)

if ~exist('neighbors')    
    bufsize=16384;
    % Read in node neighbor file generated from caret -surface-topology-neighbors
    [neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
        neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
        textread(['/data/cn/data1/scripts/CIFTI_RELATED/Resources/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
        'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
    neighbors = neighbors+1;
end

if(neighdist == 3 && exist('./input/NborDist3_fs_LR_32k.mat', 'file'))
    % load precomputed neigbours to save time
    load('input/NborDist3_fs_LR_32k.mat', 'NborDist3')
    precomFlag = 1;
else
    precomFlag = 0;
end

minimametric = zeros(size(metric), 'single');
for i = 1:size(metric, 1)
    if(precomFlag == 0)
        nodeneigh = i;
        newneigh = nodeneigh;
        curneigh = newneigh;
        
        for n = 1:neighdist
            for t = 1:length(curneigh)
                try
                    newneigh = [newneigh neighbors(curneigh(t),2:end)];
                catch
                    1;
                end
                newneigh(isnan(newneigh)) = [];
            end
            curneigh = setdiff(newneigh,nodeneigh);
            nodeneigh = union(nodeneigh,newneigh,'stable');
        end
    else
        nodeneigh = nonzeros(NborDist3(i,:));
    end
    
    nodeneighval_all = metric(nodeneigh(2:end),:);
    origval_all = repmat(metric(nodeneigh(1),:),[length(nodeneigh)-1 1]);
    minval = sum(origval_all<nodeneighval_all,1);
    minimametric(i,:) = single(minval==(length(nodeneigh)-1));
end
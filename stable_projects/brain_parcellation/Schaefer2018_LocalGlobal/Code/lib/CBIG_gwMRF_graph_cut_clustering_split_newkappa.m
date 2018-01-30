function [ likelihood, results] = CBIG_gwMRF_graph_cut_clustering_split_newkappa(x,prams,hemisphere)

    % This function performs the actual clustering
    %
    % Input
    %   - x = input timeseries
    %   - prams = a struct containing various input parameters
    %   - hemisphere = a string indicating the hemisphere 'lh' or 'rh'
    %
    % Ouput
    %   - likelihood = likelihood of the found clustering result
    %   - results = a struct containing the labels and various additional results
    %
    % Example
    %   - [ likelihood, results] = CBIG_gwMRF_graph_cut_clustering_split_newkappa(vol,prams,'lh')
    %
    %Written by Alexander Schaefer and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

    if(prams.data_normalization==1)
        x=bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));% Normalize the data to the unit sphere
    end

    avg_mesh = CBIG_ReadNCAvgMesh(hemisphere, prams.fsaverage, 'inflated', 'cortex');
    cortex_vertices=max(size(find(avg_mesh.MARS_label==2)));
    
    if(~prams.potts) 
        if (strcmp(hemisphere,'lh'))
            [Neighborhood]=CBIG_build_sparse_gradient(avg_mesh,prams.lh_grad_file,prams);
        elseif(strcmp(hemisphere,'rh'))
            [Neighborhood]=CBIG_build_sparse_gradient(avg_mesh,prams.rh_grad_file,prams);
        end
        fprintf(prams.fileID,'with gradient \n');
    else%% Potts Model
        [Neighborhood]=CBIG_build_sparse_neighborhood(avg_mesh);
        fprintf(prams.fileID,'with potts \n');
    end
    
    
    %build border matrix which will be the gradient prior
    idx_cortex_vertices=(avg_mesh.MARS_label==2);

    if (strcmp(hemisphere,'lh'))
        load(prams.lh_grad_file);
    elseif(strcmp(hemisphere,'rh'))
        load(prams.rh_grad_file);
    end

    grad_matrix=border_matrix(:,idx_cortex_vertices);

    %%%
    rng(prams.seed,'twister'); % fixate the random initialization for each prams.seed
    %%%%
    fprintf(prams.fileID,'seed number %i \n',prams.seed);
    %%%

    [likeli]=CBIG_initialize_clusters(prams.cluster,prams.dim,x,grad_matrix);% random initialization over the cortex, except vertices with gradient>0.05
    prams.graphCutIterations=10;
    gamma=zeros([1,prams.cluster])+prams.start_gamma;
    %%%%%%%%%%%%%%%%
    
    [label,Enew,D,S]=CBIG_compute_labels(cortex_vertices,likeli,Neighborhood,prams);% compute labels the first time
    [label]=CBIG_assign_empty_clusters(label,prams.cluster,grad_matrix,gamma,prams);% if a cluster is of size zero then we randomly reassign this cluster to a vertex
   
    
    [ likelihood, results,label] = CBIG_gwMRF_graph_cut_clustering_split_standard(x,prams,hemisphere,label,Neighborhood,gamma,1,grad_matrix);% 0.01 correspondes to 1% termination criteria
    initial_label=label;
    for i=0:2:4
        for j=1:1000
            gamma_head=CBIG_UpdateGamma(gamma,label,prams,hemisphere);
            if((j>1)&&(isequal(gamma,gamma_head))) % Even after at least one optimization there was no change in gamma
                break
            else
                gamma=gamma_head;
                i=i
                prams.graphCutIterations=10^(i+2);
                [ likelihood, results,label] = CBIG_gwMRF_graph_cut_clustering_split_standard(x,prams,hemisphere,label,Neighborhood,gamma,10^-i,grad_matrix);% 0.01 correspondes to 1% termination criteria
                gamma=results.gamma; % in case some cluster was empty
            end
        end
    end
    results.initial_full_label=zeros(1,max(size(avg_mesh.vertices)));
    results.initial_full_label(avg_mesh.MARS_label==2)=initial_label(1:cortex_vertices);
    
end

function [ gamma_head ] = CBIG_UpdateGamma(gamma,label,prams,hemisphere)
    avg_mesh = CBIG_ReadNCAvgMesh(hemisphere, prams.fsaverage, 'inflated', 'cortex');
    
    cortex_vertices=max(size(find(avg_mesh.MARS_label==2)));
    full_label=zeros(1,max(size(avg_mesh.vertices)));
    full_label(avg_mesh.MARS_label==2)=label(1:cortex_vertices);
    
    [lh_ci,rh_ci,lh_sizes,rh_sizes]=CBIG_gwMRF_generate_components(avg_mesh,avg_mesh,full_label,full_label);%%we are only interested in one hemisphere /however overhead for both should be small /constant :)
    lh_ci=lh_ci(avg_mesh.MARS_label==2);
    for i=1:prams.cluster
        idx=find(label==i);
        uniquevals=unique(lh_ci(idx));
        if(length(uniquevals)>1)
            binVector(i)=1;
        else
            binVector(i)=0;
        end
    end
    gamma_head=gamma;
    gamma_head(binVector==1 & gamma==0)=1000;
    gamma_head(binVector==1 & gamma~=0)=prams.reduce_speed*gamma_head(binVector==1 & gamma~=0);
    gamma_head
    number_of_components=length(unique(lh_ci))
end

function [ likelihood, results,label] = CBIG_gwMRF_graph_cut_clustering_split_standard(x,prams,hemisphere,label,Neighborhood,gamma,termination,grad_matrix)
    %UNTITLED Summary of this function goes here
    %   x - fMRI data
    %   prams - struct containing all the necessary parameters. also alex_set_prams.m
    
    Enew=1000000000000000000;
    avg_mesh = CBIG_ReadNCAvgMesh(hemisphere, prams.fsaverage, 'inflated', 'cortex');
    cortex_vertices=max(size(find(avg_mesh.MARS_label==2)));
    for j=1:(prams.iterations)
       Eold=Enew;
        if(prams.kappa_vector==1)
            [likeli,max_max_likeli,kappa,mu]=CBIG_compute_likelihood_vector_of_kappa(x,label',prams.cluster,prams.dim); % compute current likelihood with a one kappa for each cluster
        else
            [likeli,max_max_likeli,kappa,mu]=CBIG_compute_likelihood(x,label',prams.cluster,prams.dim); % compute current likelihood with one kappa for all cluster
        end
       [likeli_pos,max_max_likeli_local]=CBIG_get_position(prams,label,hemisphere,gamma); % compute spatial likelihood
       [label,Enew,E_current_D,E_current_S]=CBIG_compute_labels_spatial(cortex_vertices,likeli,Neighborhood,prams,likeli_pos);
       Enew=(Enew-max_max_likeli*cortex_vertices);
       E_current_D=(E_current_D-max_max_likeli*cortex_vertices);
 
       
       fprintf(prams.fileID,'improvement after %i iterations of %f percent \n' ,j,(Eold/Enew - 1)*100); % times 100 to make percentage, write in external text file
       fprintf(prams.fileID,'Smoothcost: %f, DataCost: %f, Energy %f \n' ,E_current_S,E_current_D,Enew);
       Eold=Eold
       Enew=Enew
       abs(Eold/Enew-1)
       if (length(unique(label))<prams.cluster)%check if there are empty cluster
            fprintf(prams.fileID,'empty cluster \n');
            [label,gamma]=CBIG_assign_empty_clusters(label,prams.cluster,grad_matrix,gamma,prams);% if a cluster is of size zero then we randomly reassign this cluster to a vertex
       elseif ((abs(Eold/Enew-1) *100) < termination)
           fprintf(prams.fileID,'hit termination of %f  \n' ,termination);
           termination=termination
           break
       else
           fprintf(prams.fileID,'missed termination of %f  \n' ,termination);
       end
       j
    end
    
    %%%% Prepare Output
    likeli=likeli+max_max_likeli;
    likeli_pos=likeli_pos+max_max_likeli_local;
    likelihood=mean(likeli([label']));
    idx=sub2ind(size(likeli),[label'],[1:max(size([label']))]); % pointwise output
    final_likeli=likeli(idx);

    %%%%%%%%%%%%%%%%%%%
    results.D=E_current_D;
    results.S=E_current_S;
    results.UnormalizedE=Enew;
    results.final_likeli=final_likeli;
    results.kappa=kappa;
    results.E= CBIG_ComputeNormalizedEnergy(kappa,mu,x,results.S,label); % to account for differences in concentration, closer to MAP
    results.gamma=gamma;
    results.mu=mu;
    results.likeli_pos=likeli_pos(idx);
    
    results.full_label=zeros(1,max(size(avg_mesh.vertices)));
    results.full_label(avg_mesh.MARS_label==2)=label(1:cortex_vertices);
end

function [ NormalizedEnergy] = CBIG_ComputeNormalizedEnergy(kappa,mu,y,smoothcost,label)
        d=size(y,2);
        n = size(y,1);
        t1=1:1:n;
        t2=label';
        weighted_data=mu'*y';
        clear y;
        datacost=bsxfun(@times,kappa',weighted_data);
        datacost=datacost(sub2ind(size(datacost),t2,t1));%indicate to which cluster each vertex belongs
        NormalizedEnergy=-sum(datacost)+smoothcost;
        
end
%-------------------------------------------------------------------------
function [likeli,initial_assigned,max_max_likeli]=CBIG_initialize_clusters(k,d,x,grad)% using randomly selected vertices for initialization, which have gradient less than 0.05
        grad=mean(grad,1);%go from edges to vertices,(only cortex vertices)
        low_grad_idx=find(grad<0.05); % we only randomly assign to vertices with low gradient values
        indices=datasample(1:length(low_grad_idx),k,'Replace',false);
        initial_assigned=low_grad_idx(indices);
        %miu=bsxfun(@times,x(initial_assigned,:),1./sqrt(sum((x(initial_assigned,:)).^2)));
        norm=1./sqrt(sum((x(initial_assigned,:)).^2,2));
        for i=1:length(initial_assigned)
            miu(i,:)=x(initial_assigned(i),:).*norm(i);
        end
        if((d>2000)&&(d<10000))
            kappa=1800;
        elseif(d>10000)
            kappa=12500;
        else
            kappa = 500;
        end
        likeli=CBIG_Cdln(kappa,d) +kappa*miu*x';
        %likeli=kappa*miu'*x'; % for the initialization we ignore the normalization term
        max_max_likeli=max(max(likeli));%global setoff
        likeli=likeli-max_max_likeli;%
end


%-----------------------------------------------------------------------
function [label,gamma]=CBIG_assign_empty_clusters(label,k,input_grad,gamma,prams)    %%%assigned empty cluster to random vertices
    grad=mean(input_grad,1);%go from edges to vertices,(only cortex vertices)
    low_grad_idx=find(grad<0.05); % we only randomly assign to vertices with low gradient values
    empty=[];
    for i=1:k
        idx=find(label==i);
        if(min(size(idx))==0) % if empty
            empty=[empty,i];
        end
    end
    assigned_vertices=[];
    if(min(size(empty)>0))%%assign new labels,gamma
        assigned_vertices=low_grad_idx(datasample(1:max(size(low_grad_idx)),max(size(empty)),'Replace',false));%%%random reassignement on low gradient vertices
        label(assigned_vertices)=empty;
    end
    gamma(empty)=prams.start_gamma; % cluster was reinitialized, clear gamma
end

%-----------------------------------------------------------------------
function [likeli,max_max_likeli]=CBIG_get_position(prams,label,hemisphere,gamma)% computes the loglikehood
    avg_mesh = CBIG_ReadNCAvgMesh(hemisphere, prams.fsaverage, 'sphere', 'cortex');
    data=avg_mesh.vertices(:,avg_mesh.MARS_label==2)';
    data=bsxfun(@rdivide,data,sqrt(sum(data.^2,2)));
    [likeli,max_max_likeli]=CBIG_compute_likelihood_given_concentration(data,label',prams.cluster,3,gamma);
end

%-----------------------------------------------------------------------
function [likeli,max_max_likeli]=CBIG_compute_likelihood_given_concentration(x,label,k,d,gamma)% computes the loglikehood
    n = size(x,1);
    t1=1:1:n;
    t2=label;
    binary_matrix=zeros(size(x,1),k);
    binary_matrix(sub2ind(size(binary_matrix),t1,t2))=1;%% we use a nxk binary maytrix to indicate to which cluster each vertex belongs
    nu=bsxfun(@times,(x'*binary_matrix),1./sqrt(sum((x'*binary_matrix).^2)));
    gammas(gamma==0)=0;
    gammas(gamma~=0)=CBIG_Cdln(single(gamma(gamma~=0)),d);
    likeli=bsxfun(@plus,gammas,bsxfun(@times,gamma,(nu'*x')'));%compute likelihood c(gamma) + gamma*nu'*x'
    max_max_likeli=max(max(likeli));
    likeli=likeli-max_max_likeli;%global setoff to get mincut and avoid maxcut (np hard)
    
end

%-----------------------------------------------------------------------
function [likeli,max_max_likeli,kappa,miu]=CBIG_compute_likelihood(x,label,k,d)% computes the loglikehood
    n = size(x,1);
    t1=1:1:n;
    t2=label;
    binary_matrix=zeros(size(x,1),k);
    binary_matrix(sub2ind(size(binary_matrix),t1,t2))=1;%% we use a nxk binary maytrix to indicate to which cluster each vertex belongs
    %miu=bsxfun(@times,(x'*binary_matrix),1./sqrt(sum((x'*binary_matrix).^2)));
    test=x'*binary_matrix;
    norm=1./sqrt(sum((test).^2));
    for i=1:k
        miu(i,:)=test(:,i).*norm(i)';
    end
    r_bar=sum(sum(binary_matrix.*(x*miu)))/n;
    kappa=CBIG_invAd(d,double(r_bar));
    %kappa=(r_bar*d-power(r_bar,3))/(1-r_bar*r_bar);
    likeli=CBIG_Cdln(kappa,d) +kappa*miu'*x';
    max_max_likeli=max(max(likeli));%global setoff
    likeli=likeli-max_max_likeli;%
    
end

%-----------------------------------------------------------------------
function [likeli,max_max_likeli,kappa,miu]=CBIG_compute_likelihood_vector_of_kappa(x,label,k,d)% computes the loglikehood
    n = size(x,1);
    t1=1:1:n;

    for i=1:k
        idx=find(label==i);
        indexed_data=x(label==i,:);
        miu=(mean(indexed_data,1))/norm(mean(indexed_data,1));
        %miu=bsxfun(@times,(mean(x(idx,:),1)),1./sqrt(sum((mean(x(idx,:),1)).^2)));
        r_bar=sum((indexed_data*miu'))/length(idx);
        if(length(idx)==1)
            if(r_bar>=0.99) % this happens when there is only vertex with the particular label, would otherwise inflate the likelihood 
                r_bar=0.99;
            end
        end
        kappa(i)=CBIG_invAd(d,double(r_bar));
    end
    t2=label;
    binary_matrix=zeros(size(x,1),k);
    binary_matrix(sub2ind(size(binary_matrix),t1,t2))=1;%% we use a nxk binary maytrix to indicate to which cluster each vertex belongs
    %miu=bsxfun(@times,(x'*binary_matrix),1./sqrt(sum((x'*binary_matrix).^2)));
    test=x'*binary_matrix;
    normdat=1./sqrt(sum((test).^2));
    for i=1:k
        miu(i,:)=test(:,i).*normdat(i)';
    end
    miu=miu';% transpose to be consistent with the rest of the code, in case we want to go back to using bsxfun
    likeli=bsxfun(@plus,Cdln(kappa,d)',bsxfun(@times,kappa',miu'*x'));
    max_max_likeli=max(max(likeli));%global setoff
    likeli=likeli-max_max_likeli;%
    
end

%-----------------------------------------------------------------------
function [out exitflag] =CBIG_invAd(D,rbar)

    outu = (D-1)*rbar/(1-rbar^2) + D/(D-1)*rbar;

    [i] = besseli(D/2-1,outu);

    if ((i == Inf)||(isnan(i))||(i==0))
        out = outu - D/(D-1)*rbar/2;
        exitflag = Inf;
    else
        [outNew fval exitflag]  = fzero(@(argum) CBIG_Ad(argum,D)-rbar,outu);
        if exitflag == 1
            out = outNew;
        else
            out = outu - D/(D-1)*rbar/2;
        end
    end
end
%-----------------------------------------------------------------------
function out =CBIG_Ad(in,D)
    % out = CBIG_Ad(in,D)
    out = besseli(D/2,in) ./ besseli(D/2-1,in);
end
%-----------------------------------------------------------------------
function out = CBIG_Cdln(k,d)

    % Computes the logarithm of the partition function of vonMises-Fisher as
    % a function of kappa

    sizek = size(k);
    k = k(:);

    out = (d/2-1).*log(k)-log(besseli((d/2-1)*ones(size(k)),k));
    if(d<1000)
        k0 = 10;
    elseif(d>1000&&d<2000)
        k0 = 500;
    elseif(d>=2000&&d<5000)
        k0=1600;
    elseif(d>=30000&&d<40000)
        k0=12500;
    elseif(d>=40000&&d<50000)
        k0=13700;
    elseif(d>=150000)
        k0=51200;    
    end

    fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
    nGrids = 1000;
    if isinf(fk0)% more general function works at least up to d=350k
        k0=0.331*d;
        fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
    end
    maskof = find(max((k>k0),(isinf(out))));
    nkof = length(maskof);

    % The kappa values higher than the overflow

    if (nkof > 0 )

        kof = k(maskof);

        ofintv = (kof - k0)/nGrids;
        tempcnt = (1:nGrids) - 0.5;
        ks = k0 + repmat(tempcnt,nkof,1).*repmat(ofintv,1,nGrids);
        adsum = sum( 1./((0.5*(d-1)./ks) + sqrt(1+(0.5*(d-1)./ks).^2)) ,2);

        % with correction:
        % 1./((0.5*(d-1)./x) + sqrt(1+(0.5*(d-1)./x).^2*(d/2/(d-1)^2)))

        out(maskof) =  fk0 - ofintv .* adsum;

    end

    out = reshape(out,sizek);
end



function [Neighborhood]=CBIG_build_sparse_neighborhood(avg_mesh)
    idx_cortex_vertices=find(avg_mesh.MARS_label==2);
    vertices=max(size(avg_mesh.vertexNbors));
    r = reshape(repmat(13:vertices,6,1),1,6*(vertices-12));
    Neighborhood=sparse(double(r),double(reshape(avg_mesh.vertexNbors(1:6,13:end),[1,6*(vertices-12)])),double(ones(size(r))));
    for i=1:12
        Neighborhood(i,avg_mesh.vertexNbors(1:5,i))=1; %% account for the first 12 vertices having only 5 neighbors, might be slow
    end
    Neighborhood=Neighborhood(idx_cortex_vertices,idx_cortex_vertices); % remove medial wall
end

function [Neighborhood]=CBIG_build_sparse_gradient(avg_mesh,gradient,prams)
    idx_cortex_vertices=find(avg_mesh.MARS_label==2);
    vertices=max(size(avg_mesh.vertexNbors));
    load(gradient)
    r = reshape(repmat(13:vertices,6,1),1,6*(vertices-12));
    Neighborhood=sparse(double(r),double(reshape(avg_mesh.vertexNbors(1:6,13:end),[1,6*(vertices-12)])),reshape(StableE(border_matrix(:,13:end),prams.exponential),(size(r)))); %we build a sparse matrix,stableE will ad epsilon if values are 0, be aware if you use something else
    for i=1:12
        Neighborhood(i,avg_mesh.vertexNbors(1:5,i))=StableE(border_matrix(1:5,i),prams.exponential); %% account for the first 12 vertices having only 5 neighbors, might be slow
    end
    Neighborhood=Neighborhood(idx_cortex_vertices,idx_cortex_vertices); % remove medial wall
end

function [label,Enew,D,S]=CBIG_compute_labels(cortex_vertices,probs,Neighborhood,prams)

    h = GCO_Create(cortex_vertices,prams.cluster);
    data_cost=single(-(1/cortex_vertices)*((probs(:,1:cortex_vertices)+eps))+eps);

    GCO_SetDataCost(h,data_cost);% minimize negativ loglikehood

    Smooth_cost=ones(prams.cluster,prams.cluster,'single');
    Smooth_cost(1:(prams.cluster+1):end)=0;% sets the diagonal to 0
    if(~prams.potts)
        GCO_SetNeighbors(h,Neighborhood);
    end
    GCO_SetSmoothCost(h,Smooth_cost*single(prams.smoothcost*(1/cortex_vertices))); % increase the smooth cost, to get smoother results
    GCO_Expansion(h,prams.graphCutIterations);% maximum number of iterations is 1000, this might be changed but should be a finite number
    label=GCO_GetLabeling(h);
    [Enew,D,S] = GCO_ComputeEnergy(h);

    GCO_Delete(h);
    Enew=Enew*cortex_vertices;% scale the cost directly back
    D=D*cortex_vertices;% scale the cost directly back
    S=S*cortex_vertices;% scale the cost directly back
end

function [label,Enew,D,S]=CBIG_compute_labels_spatial(cortex_vertices,probs,Neighborhood,prams,likeli_pos)

    h = GCO_Create(cortex_vertices,prams.cluster);
    normalize=cortex_vertices;% do the normalization dynamically
    data_cost=single(-(1/normalize)*((probs+eps))+eps);
    pos_cost=single(-(1/normalize)*((likeli_pos+eps))+eps);
    GCO_SetDataCost(h,data_cost+pos_cost');
    
    
    Smooth_cost=ones(prams.cluster,prams.cluster,'single');
    Smooth_cost(1:(prams.cluster+1):end)=0;% sets the diagonal to 0
    if(~prams.potts)
        GCO_SetNeighbors(h,Neighborhood);
    end
    GCO_SetSmoothCost(h,Smooth_cost*single(prams.smoothcost*(1/normalize))); % increase the smooth cost, to get smoother results
    GCO_Expansion(h,prams.graphCutIterations);% maximum number of iterations should be a finite number
    label=GCO_GetLabeling(h);
    [Enew,D,S] = GCO_ComputeEnergy(h);

    GCO_Delete(h);
    Enew=Enew*normalize;% scale the cost directly back
    D=D*normalize;% scale the cost directly back
    S=S*normalize;% scale the cost directly back
end

function x = CBIG_StableE(x,k)
%%% this function realizes: exp(-k * gradient) - exp(-k)
%%% performs several checks to avoid infinite or nan values
    x(x > 1) = 1;
    x(x < 0) = 0;
    x = CBIG_our_exp(x,k);
    x(isinf(x) & x > 0) = CBIG_our_exp(1-eps,k);
    x(~isreal(x) & real(x) > 0) = CBIG_our_exp(1-eps,k);

    x(isinf(x) & x < 0) = CBIG_our_exp(0+eps,k);
    x(~isreal(x) & real(x) < 0) = CBIG_our_exp(0+eps,k);
end

function y=CBIG_our_exp(p,k)
    y=exp(-k*p)-exp(-k);
end

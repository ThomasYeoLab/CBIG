function [lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation( ...
         project_dir, mesh, num_sess, num_clusters, subid, w, c, subject_set)
% [lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation(project_dir,mesh,num_sess,num_clusters,subid,w,c)
%
% This script will estimate the individual-level parcellation with 
% num_clusters networks for a single subject using num_sess sessions. The
% weight of group spatial prior should be specified by w, the weight of MRF
% smoothness prior should be specified by c.
%
% To generate the individual-level parcellation, we assume the group priors
% are already estimated by running CBIG_MSHBM_estimate_group_priors.m. The 
% estimated group priors should be saved in 
% project_dir/priors/Params_Final.mat.
%
% We also assume the functional connectivity profiles of test subjects with 
% num_sess sessions are already generated. The lists of functional 
% connectivity profiles should be saved in 
% project_dir/profile_list/<subject_set>/lh_sess<?>.txt
% project_dir/profile_list/<subject_set>/rh_sess<?>.txt
% for data in 'fsaverage4/fsaverage5/fsaverage6/fsaverage'
% or
% project_dir/profile_list/<subject_set>/sess<?>.txt
% for data in 'fs_LR_32k'.
%
% Input:
%   - project_dir:
%
%     The project directory.
%     1) project_dir/priors/Params_Final.mat.
%        contains the results of group priors estimated by running 
%        CBIG_MSHBM_estimate_group_priors.m on training set.
%        Params_Final.mat contain a struct variable Params with fields
%        corresponding to each group prior:
%        1) Params.epsil
%           Inter-subject functional connectivity variability. 
%        2) Params.mu
%           Group-level connectivity profiles for each network. 
%        3) Params.sigma
%           Intra-subject functional connectivity variability. 
%        4) Params.theta
%           Spatial prior which denotes the probability of each network 
%           occurring at each location. 
%                                                          
%     2) project_dir/profile_list/<subject_set>/lh_sess<?>.txt
%        project_dir/profile_list/<subject_set>/rh_sess<?>.txt
%        contain the functional connectivity profile lists of left
%        hemisphere and right hemisphere for each session. The functional
%        connectivity profiles are assumed to be pre-computed before run
%        the current script. These profile lists are assumed to be
%        pre-generated. For S test subjects and T sessions, there
%        should be T lh_sess<?>.txt and rh_sess<?>.txt lists for data in
%        'fsaverage4/fsaverage5/fsaverage6/fsaverage', or T sess<?>.txt
%        lists for data in 'fs_LR_32k'.
%        For example:
%        project_dir/profile_list/<subject_set>/lh_sess1.txt
%        project_dir/profile_list/<subject_set>/rh_sess1.txt
%        project_dir/profile_list/<subject_set>/lh_sess2.txt
%        project_dir/profile_list/<subject_set>/rh_sess2.txt
%        or
%        project_dir/profile_list/<subject_set>/sess1.txt
%        project_dir/profile_list/<subject_set>/sess2.txt
%        for 2 sessions. Each list should contain S rows, where each row 
%        is the full file path of the functional connectivity profile for 
%        each test subject.
%
%   - mesh: (string)
%     
%     The data surface space. 'fsaverage5/fsaverage6/fsaverage' or 'fs_LR_32k'. 
%
%   - num_sess: (string)
%
%     The number of sessions the user want to use to estimate the
%     individual-level parcellation. For example, '4'.
%
%   - num_clusters: (string)
%
%     The number of networks of the parcellations. For example, '17'.
%
%   - subid: (string)
%
%     The test subject number. For example, '4' indicates the 4-th
%     subject in the project_dir/profile_list/<subject_set>/?h_sess<?>.txt or 
%     project_dir/profile_list/<subject_set>/sess<?>.txt.
%
%   - w: (string)
%   
%     The weight of group spatial prior Params.theta. For example, '100'.
%     A large w indicates strong weight of Params.theta. The estimated
%     individual-level parcellation will be very similar to the group-level
%     parcellation with very large w.
%
%   - c: (string)
%
%     The weight of MRF smoothness prior. For example, '50'. A large c
%     indicates more penalty for neighboring vertices being assigned to
%     different networks.
%
%   - subject_set: (string)
%
%     'validation_set' or 'test_set'. If this input argument is not given,
%     the default value will be 'test_set'. This argument is used if the
%     user has a validation set to determine parameters w and c.
%
% Output:
%   
%   - lh_labels, rh_labels:
%
%     The labels of left and right hemisphere of individual parcellation. 
%
% Example:
% [lh_labels, rh_labels] = CBIG_MSHBM_generate_individual_parcellation(project_dir,'fsaverage5','5','17','4','100','50')
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));

if(nargin < 8)
    subject_set = 'test_set';
end

%% setting parameters
setting_params.num_sub = '1'; 
setting_params.num_clusters = str2double(num_clusters);
setting_params.epsilon = 1e-4;
setting_params.mesh = mesh;
setting_params.subid = str2double(subid);
setting_params.w = str2double(w); % weight of spatial prior theta
setting_params.c = str2double(c); % weight of MRF smoothness prior
% When generating individual parcellation, each time we will only generate the parcellation for one subject.
setting_params.num_sub = 1; 
setting_params.num_session = str2double(num_sess);
   
%% read group priors
group_prior_file = fullfile(project_dir, 'priors', 'Params_Final.mat');   
Ini = load(group_prior_file);

%% read data
data = fetch_data(project_dir, setting_params.num_session, setting_params.subid, ...
    setting_params.mesh,subject_set);
fprintf('read data...DONE!\n');
setting_params.dim = size(data.series,2) - 1;
setting_params.num_verts = size(data.series,1);

if(setting_params.dim < 1200) 
    setting_params.ini_concentration = 500; 
elseif(setting_params.dim >= 1200 && setting_params.dim < 1800 )
    setting_params.ini_concentration = 650;                                          
else
    error(['Data dimension is higher than 1800,' ...
          'please set the minimal concentration parameter value ini_concentration,' ...
          'where besseli((D/2-1),ini_concentration) > 0 and relatively small\n']);
end    

%% setting estimated group priors
% mu: DxL. The group-level functional connectivity profiles of networks.
% The group prior mu will be 
Params.mu = Ini.Params.mu;

% epsil: 1xL. The inter-subject concentration parameter, which represents
% inter-subject functional connectivity variability. A large epsil_l 
% indicats low inter-subject functional connectivity variability for 
% network l.
Params.epsil = Ini.Params.epsil;

% sigma: 1xL. The intra-subject concentration parameter, which represents
% intra-subject functional connectivity variability. A large sigma_l
% indicates low intra-subject functional connectivity variability for
% network l.
Params.sigma = Ini.Params.sigma;

%theta: NxL. The spatial prior denotes the probability of networks
%occurring at each spatial location.
Params.theta = Ini.Params.theta;

%% paramter initialization

% s_psi: DxLxS. The functional connectivity profiles of L networks for S
% subjects.
Params.s_psi = repmat(Ini.Params.mu,1,1,setting_params.num_sub);

% s_t_nu: DxLxTxS. The functional connectivity profiles of L networks for S
% subjects and each subject has T sessions.
Params.s_t_nu = repmat(Ini.Params.mu,1,1,setting_params.num_session,setting_params.num_sub);

% kappa: 1xL. The inter-region concentration parameter, which represents
% inter-region functional connectivity variability. A large kappa_l 
% indicates low inter-region functional variability for network l. However,
% please note in this script, we assume kappa to be the same across 
% networks.
Params.kappa = setting_params.ini_concentration*ones(1,setting_params.num_clusters);

% s_lambda: NxLxS. The posterior probability of the individual-specific
% parcellation of each subject. 
log_vmf = permute(Params.s_t_nu,[1,2,4,3]);
log_vmf = mtimesx(data.series,log_vmf);%NxLxSxT
log_vmf = bsxfun(@times,permute(log_vmf,[2,1,3,4]),transpose(Params.kappa));%LxNxSxT
log_vmf = bsxfun(@plus,Cdln(transpose(Params.kappa),setting_params.dim),log_vmf);%LxNxSxT
log_vmf = sum(log_vmf,4);
log_vmf = permute(log_vmf,[2,1,3]);
s_lambda = bsxfun(@minus,log_vmf,max(log_vmf,[],2));
mask = repmat((sum(s_lambda,2)==0),1,setting_params.num_clusters,1);
s_lambda = exp(s_lambda);
Params.s_lambda = bsxfun(@times,s_lambda,1./sum(s_lambda,2));
Params.s_lambda(mask) = 0;

%% introduce MRF smoothness prior
if(~isempty(strfind(mesh,'fsaverage')))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex');

    lmw = transpose(find(sum(Params.s_lambda(1:setting_params.num_verts/2,:),2)==0));
    rmw = transpose(find(sum(Params.s_lambda((setting_params.num_verts/2 + 1):end,:),2)==0));
    l1 = transpose(find(sum(Params.s_lambda(1:setting_params.num_verts/2,:),2)~=0));
    l2 = transpose(find(sum(Params.s_lambda((setting_params.num_verts/2 + 1):end,:),2)~=0));
elseif(~isempty(strfind(mesh,'fs_LR_32k')))
    lh_avg_mesh = CBIG_read_fslr_surface('lh','fs_LR_32k','inflated');
    rh_avg_mesh = CBIG_read_fslr_surface('rh','fs_LR_32k','inflated');
    
    lmw = transpose(find(sum(Params.s_lambda(1:setting_params.num_verts/2,:),2)==0));
    rmw = transpose(find(sum(Params.s_lambda(setting_params.num_verts/2+1:end,:),2)==0));
    l1 = transpose(find(sum(Params.s_lambda(1:setting_params.num_verts/2,:),2)~=0));
    l2 = transpose(find(sum(Params.s_lambda(setting_params.num_verts/2+1:end,:),2)~=0));
end
    
lh_neigh = lh_avg_mesh.vertexNbors;
rh_neigh = rh_avg_mesh.vertexNbors;
lh_neigh = lh_neigh(:,l1);
rh_neigh = rh_neigh(:,l2);


for t = lmw
    lh_neigh(lh_neigh == t) = 0;
end
unique(lh_neigh);
for l = 1:length(l1)
    lh_neigh(lh_neigh == l1(l)) = l;
end
unique(lh_neigh);
for t = rmw
    rh_neigh(rh_neigh == t) = 0;
end
unique(rh_neigh);
for l = 1:length(l2)
    rh_neigh(rh_neigh == l2(l)) = l;
end

rh_neigh = rh_neigh + size(lh_neigh,2);
rh_neigh(rh_neigh == size(lh_neigh,2)) = 0;
neighborhood = [lh_neigh rh_neigh];
setting_params.neighborhood = double(neighborhood);

setting_params.V_diff = ones(size(neighborhood));
setting_params.V_same = zeros(size(neighborhood));


%% EM
Params.iter_inter = 1;
%% Intra subject variability
cost = 0;
iter_intra_em = 0;
stop_intra_em = 0;
while(stop_intra_em == 0)
    iter_intra_em = iter_intra_em + 1;
    
    fprintf('Inter-region iteration %d ...\n', iter_intra_em);
    Params.kappa = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
    Params.s_t_nu = repmat(Ini.Params.mu, 1, 1, setting_params.num_session, setting_params.num_sub);
    Params = vmf_clustering_subject_session(Params, setting_params, data);
    
    fprintf('Intra-subject variability level...\n');
    Params = intra_subject_var(Params, setting_params);

    update_cost = bsxfun(@times, Params.s_psi, permute(Params.s_t_nu, [1,2,4,3]));
    update_cost = sum(bsxfun(@times, Params.sigma, update_cost), 1);
    update_cost = bsxfun(@plus, Cdln(Params.sigma, setting_params.dim), update_cost);
    update_cost = sum(sum(sum(update_cost, 2), 3), 4) ...
        + sum(sum(bsxfun(@plus, sum(bsxfun(@times, bsxfun(@times, Params.mu, Params.s_psi), ...
        Params.epsil), 1), Cdln(Params.epsil,setting_params.dim)),2),3);
    update_cost = update_cost + sum(Params.cost_em);
    Params.Record(iter_intra_em) = update_cost;
    if(abs(abs(update_cost-cost)./cost) <= setting_params.epsilon)
        stop_intra_em = 1;
        Params.cost_intra = update_cost;
    end
    if(iter_intra_em >= 50)
        stop_intra_em = 1;
        Params.cost_intra = update_cost;
    end
    cost = update_cost;
end

%% generate parcellation
labels = zeros(setting_params.num_verts, 1);
[~,labels(sum(Params.s_lambda,2)~=0)] = max(Params.s_lambda(sum(Params.s_lambda,2)~=0, :), [], 2);
lh_labels = labels(1:setting_params.num_verts/2);
rh_labels = labels((setting_params.num_verts/2 + 1):end);
out_dir = fullfile(project_dir, 'ind_parcellation', subject_set);

if(~exist(out_dir))
    mkdir(out_dir);
end
save(fullfile(out_dir, ...
    ['Ind_parcellation_MSHBM_sub',num2str(subid),'_w',num2str(setting_params.w),'_MRF',num2str(c),'.mat']), ...
    'lh_labels','rh_labels');

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));

end

%% sub-functions

function Params=intra_subject_var(Params,setting_params)

% Intra-subject variability level

iter_intra=0;
iter_intra=iter_intra+1;
fprintf('It is inter interation %d intra iteration %d..update s_psi and sigma..\n',Params.iter_inter,iter_intra);

% update s_psi
s_psi_update = sum(bsxfun(@times, Params.s_t_nu, repmat(Params.sigma, size(Params.s_t_nu, 1), ...
    1, size(Params.s_t_nu, 3), size(Params.s_t_nu, 4))), 3);
s_psi_update = reshape(s_psi_update, size(s_psi_update, 1), size(s_psi_update, 2), ...
    size(s_psi_update, 3)*size(s_psi_update, 4));
s_psi_update = bsxfun(@plus, s_psi_update, bsxfun(@times, Params.epsil, Params.mu));
s_psi_update = bsxfun(@times, s_psi_update, 1./sqrt(sum((s_psi_update).^2)));

Params.s_psi = s_psi_update;

end

function Params=vmf_clustering_subject_session(Params,setting_params,data)

stop_em = 0;
iter_em = 0;
cost = zeros(1, setting_params.num_sub);
while(stop_em == 0)
    iter_em = iter_em + 1;
    fprintf('It is EM iteration.. %d..\n', iter_em);
    %% Mstep
    flag_nu = zeros(setting_params.num_sub,setting_params.num_session);

    stop_m = 0;
    iter_m = 0;
    fprintf('M-step..\n');
    while(stop_m == 0)
        iter_m = iter_m + 1;

        s_lambda = Params.s_lambda; %NxLxS
        kappa_update = mtimesx(data.series, permute(Params.s_t_nu, [1,2,4,3]));
        kappa_update = bsxfun(@times,s_lambda,kappa_update);
        kappa_update = sum(sum(sum(sum(kappa_update, 1), 4), 3));
        kappa_update = kappa_update./sum((setting_params.num_session.*sum(sum(s_lambda,1),3)));
        kappa_update = invAd(setting_params.dim, kappa_update);
        kappa_update = repmat(kappa_update, 1, setting_params.num_clusters);
       
        for s = 1:setting_params.num_sub
            for t = 1:setting_params.num_session

                checknu=[];
                X = data.series(:,:,s,t);
                s_lambda = Params.s_lambda(:,:,s);
                lambda_X = bsxfun(@times,kappa_update,X'*s_lambda) ...
                    + bsxfun(@times,Params.sigma,Params.s_psi(:,:,s));
                s_t_nu_update = bsxfun(@times, lambda_X, 1./sqrt(sum((lambda_X).^2)));
                checknu = diag(s_t_nu_update'*Params.s_t_nu(:,:,t,s));
                checknu_flag = (sum(1-checknu < setting_params.epsilon) < setting_params.num_clusters);
                Params.s_t_nu(:,:,t,s) = s_t_nu_update;

                if(checknu_flag < 1)
                    flag_nu(s,t) = 1;
                end
                
            end
        end
        if((sum(sum(flag_nu))==setting_params.num_sub*setting_params.num_session) ...
          && (mean(abs(Params.kappa-kappa_update)./Params.kappa) < setting_params.epsilon))
            stop_m = 1;
        end
        Params.kappa = kappa_update;
    end
    %% Estep
    fprintf('Estep..\n');
    stop_lambda = 0;
    checklam = 0;
    lambda_iter = 0;
    while stop_lambda == 0
        lambda_iter = lambda_iter + 1;

        log_vmf = permute(Params.s_t_nu,[1,2,4,3]);
        log_vmf = mtimesx(data.series,log_vmf);%NxLxSxT
        log_vmf = bsxfun(@times,permute(log_vmf,[2,1,3,4]), transpose(Params.kappa)); %LxNxSxT
        log_vmf(:,sum(log_vmf==0,1)==0)=bsxfun(@plus, ...
                                        Cdln(transpose(Params.kappa),setting_params.dim), ...
                                        log_vmf(:,sum(log_vmf==0,1)==0));%LxNxSxT
        log_vmf = sum(log_vmf,4);%NxLxS
        idx = sum(log_vmf==0,1)~=0;

        tmp_lambda = Params.s_lambda(sum(Params.s_lambda,2)~=0,:);
        V_lambda = CBIG_MSHBM_V_lambda_Product(setting_params.neighborhood, ...
            setting_params.V_same, setting_params.V_diff, double(tmp_lambda));
        V_temp = zeros(size(Params.s_lambda));
        V_temp(sum(Params.s_lambda,2)~=0,:) = single(V_lambda);
        log_vmf = bsxfun(@plus, permute(log_vmf, [2,1,3]), ...
            setting_params.w * log(Params.theta) - 2 * setting_params.c * V_temp);
      
        s_lambda = bsxfun(@minus, log_vmf, max(log_vmf, [], 2));
        s_lambda = exp(s_lambda);
        update_s_lambda = bsxfun(@times, s_lambda, 1./sum(s_lambda, 2));
        
        update_s_lambda = permute(update_s_lambda, [2,1,3]);
        update_s_lambda(:,idx) = 0;
        update_s_lambda = permute(update_s_lambda, [2,1,3]);
        
        checklam_update = mean(mean(abs(update_s_lambda - Params.s_lambda)));
        if(abs(checklam_update - checklam) <= setting_params.epsilon)
            stop_lambda = 1;
        end
        checklam = checklam_update;
        if(lambda_iter > 100)
            stop_lambda = 1;
            warning('lambda can not converge');
        end
        Params.s_lambda = update_s_lambda;
    end
        
    %% em stop criteria
    for s = 1:setting_params.num_sub    
        for t = 1:setting_params.num_session
            X = data.series(:,:,s,t);
            setting_params.num_verts = size(X,1);
            log_vmf = vmf_probability(X, Params.s_t_nu(:,:,t,s), Params.kappa);
            if(t == 1)
                log_lambda_prop = log_vmf;
            else
                log_lambda_prop = log_lambda_prop + log_vmf;
            end
        end
        
        theta_cost = Params.theta;
        log_theta_cost = log(theta_cost);
        log_theta_cost(isinf(log_theta_cost)) = log(eps.^20);

        s_lambda_cost = Params.s_lambda(:,:,s);
        log_s_lambda_cost = log(s_lambda_cost);
        log_s_lambda_cost(isinf(log_s_lambda_cost)) = log(eps.^20);

        update_cost(:,s) = sum(sum(s_lambda_cost.*log_lambda_prop)) ...
                         + sum(sum(s_lambda_cost.*setting_params.w.*log_theta_cost)) ...
                         - sum(sum(s_lambda_cost.*log_s_lambda_cost)) ...
                         - sum(sum(s_lambda_cost.*setting_params.c.*V_temp));    
    
    end
    sub_set = find((abs(abs(update_cost-cost)./cost)>1e-4)==0);
    if(length(sub_set)==setting_params.num_sub)
        stop_em = 1;
        Params.cost_em = cost;
    end
    if(iter_em>100)
        stop_em = 1;
        Params.cost_em = update_cost;
        warning('vem can not converge');
    end
    cost = update_cost;
end
fprintf('EM..Done\n');
end




function log_vmf = vmf_probability(X,nu,kap)

% log of von Mises-Fisher distribution
% X: input data
% nu: mean direction
% kap: concentration parameter. kap is a 1xL vector

dim = size(X,2) - 1;
log_vmf = bsxfun(@plus,Cdln(kap,dim),bsxfun(@times,kap,X*nu));
end

function out = Ad(in,D)
out = besseli(D/2,in) ./ besseli(D/2-1,in);
end

function out = Cdln(k,d,k0)
k = double(k);

% Computes the logarithm of the partition function of vonMises-Fisher as
% a function of kappa

sizek = size(k);
k = k(:);

out = (d/2-1).*log(k)-log(besseli((d/2-1)*ones(size(k)),k));
if(d<1200)
    k0 = 500;
elseif(d>=1200 && d<1800)
    k0 = 650;
else
    error('dimension is too high, need  to specify k0');
end
fk0 = (d/2-1).*log(k0)-log(besseli(d/2-1,k0));
nGrids = 1000;

maskof = find(k>k0);
nkof = length(maskof);

% The kappa values higher than the overflow

if nkof > 0

    kof = k(maskof);

    ofintv = (kof - k0)/nGrids;
    tempcnt = (1:nGrids) - 0.5;
    ks = k0 + repmat(tempcnt,nkof,1).*repmat(ofintv,1,nGrids);
    adsum = sum( 1./((0.5*(d-1)./ks) + sqrt(1+(0.5*(d-1)./ks).^2)) ,2);

    out(maskof) =  fk0 - ofintv .* adsum;

end

out = single(reshape(out,sizek));
end

function [outu] = invAd(D,rbar)

rbar = double(rbar);

outu = (D-1).*rbar./(1-rbar.^2) + D/(D-1).*rbar;

[i] = besseli(D/2-1,outu);


if ((i == Inf)||(isnan(i)) || (i==0))
    out = outu - D/(D-1)*rbar/2;
    exitflag = Inf;
else
    [outNew, fval exitflag]  = fzero(@(argum) Ad(argum,D)-rbar,outu);
    if exitflag == 1
        out = outNew;
    else
        out = outu - D/(D-1)*rbar/2;
    end
end
end


function data = fetch_data(project_dir,num_session,subid,mesh,subject_set)

% read in input functional connectivity profiles
if(~isempty(strfind(mesh,'fs_LR_32k')))
    load(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', ...
        'lib', 'fs_LR_32k_medial_mask.mat'));
    for t = 1:num_session
        data_profile = fullfile(project_dir,'profile_list',subject_set,['sess' num2str(t) '.txt']);
        profile_name = table2cell(readtable(data_profile,'Delimiter',' ','ReadVariableNames',false));
        for i = 1
            fprintf('Session %d...It is subj %d...\n',t,subid);
            avg_file = profile_name{subid,1};

            [~, series, ~] = CBIG_MSHBM_read_fmri(avg_file);       
            series(~medial_mask, :) = 0;

            series = bsxfun(@minus,series,mean(series, 2));
            series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), ...
                                         sqrt(sum(series(all(series,2)~=0,:).^2,2)));
            data.series(:,:,i,t) = series;
        end
    end

elseif(~isempty(strfind(mesh,'fsaverage')))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex'); 
    for t = 1:num_session
        lh_data_profile = fullfile(project_dir,'profile_list',subject_set,...
            ['lh_sess' num2str(t) '.txt']);
        rh_data_profile = fullfile(project_dir,'profile_list',subject_set,...
            ['rh_sess' num2str(t) '.txt']);
        lh_profile_name = table2cell(readtable(lh_data_profile,'Delimiter',...
            ' ','ReadVariableNames',false));
        rh_profile_name = table2cell(readtable(rh_data_profile,'Delimiter',...
            ' ','ReadVariableNames',false));
        for i = 1
            fprintf('Session %d...It is subj %d...\n',t,subid);
            lh_avg_file = lh_profile_name{subid,1};
            rh_avg_file = rh_profile_name{subid,1};

            [~, lh_series, ~] = CBIG_MSHBM_read_fmri(lh_avg_file);
            [~, rh_series, ~] = CBIG_MSHBM_read_fmri(rh_avg_file);

            lh_series(lh_avg_mesh.MARS_label == 1,:) = 0;
            rh_series(rh_avg_mesh.MARS_label == 1,:) = 0;

            series = [lh_series;rh_series];

            series = bsxfun(@minus,series,mean(series, 2));
            series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), ...
                                         sqrt(sum(series(all(series,2)~=0,:).^2,2)));
            data.series(:,:,i,t) = series;
        end
    end
    
end 
    
end

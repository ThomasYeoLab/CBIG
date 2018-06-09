function Params = CBIG_MSHBM_estimate_group_priors(project_dir,mesh,num_sub,num_sess,num_clusters)

% Params = CBIG_MSHBM_generate_individual_parcellation_FS(project_dir,mesh,sub,sess,num_clusters)
%
% This script will estimate group priors for individual-level parcellation
% generation. The estimated group priors include:
% 1) inter-subject functional connectivity variability -- Params.epsil
% 2) group-level connectivity profiles for each network -- Params.mu
% 3) intra-subject functional connectivity variability -- Params.sigma
% 4) spatial prior which denotes the probability of each network occurring
%    at each location -- Params.theta
%
% To estimate the group priors, we assume the group-level
% parcellation algorithm (Yeo et al., 2011) was already done to obtain
% initialization parameters, which are saved in
% project_dir/group/group.mat. 
%
% The functional connectivity profiles of num_sub training subjects with 
% num_sess sessions are assumed already generated. The lists of functional 
% connectivity profiles are assumed in 
% project_dir/profile_list/training_set/lh_sess<?>.txt
% project_dir/profile_list/training_set/rh_sess<?>.txt.
%
% Input:
%   - project_dir:
%
%     The project directory.
%     1) project_dir/group/group.mat 
%        contains the results of group-level parcellation algorithm. This
%        file is assumed to be pre-computed before run the current script.
%        "group.mat" includes the von Mises-Fisher mean direction parameter
%        "clustered.mtc". If the dimension of functional connectivity profile
%        is NxD, then "clustered.mtc" should be num_clustersxD.
%
%     2) project_dir/profile_list/training_set/lh_sess<?>.txt
%        project_dir/profile_list/training_set/rh_sess<?>.txt
%        contain the functional connectivity profile lists of left
%        hemisphere and right hemisphere for each session. The functional
%        connectivity profiles are assumed to be pre-computed before run
%        the current script. These profile lists are assumed to be
%        pre-generated. For S training subjects and T sessions, there
%        should be T lh_sess<?>.txt and rh_sess<?>.txt lists for data in
%        'fsaverage4/fsaverage5/fsaverage6/fsaverage', or T sess<?>.txt
%        lists for data in 'fs_LR_32k'.
%        For example:
%        project_dir/profile_list/training_set/lh_sess1.txt
%        project_dir/profile_list/training_set/rh_sess1.txt
%        project_dir/profile_list/training_set/lh_sess2.txt
%        project_dir/profile_list/training_set/rh_sess2.txt
%        or
%        project_dir/profile_list/training_set/sess1.txt
%        project_dir/profile_list/training_set/sess2.txt
%        for 2 sessions. Each list should contain S rows, where each row 
%        is the full file path of the functional connectivity profile for 
%        each test subject.
%
%   - mesh: (string)
%     
%     The data surface space. 'fsaverage5/fsaverage6/fsaverage' or 
%     'fs_LR_32k'. 
%
%   - num_sub: (string)
%
%     The number of subjects the user want to use to estimate the group
%     priors. For example, '40'.
%
%   - num_sess: (string)
%
%     The number of sessions the user want to use to estimate the group
%     priors. For example, '4'.
%
%   - num_clusters: (string)
%
%     The number of networks of the parcellations. For example, '17'.
%
% Output:
%   
%   - Params: (struct)
%     
%     D: data dimension
%     N: #vertices
%     L: #networks
%     S: #subjects
%     T: #sessions
%
%     Params.mu: DxL. 
%     The group-level functional connectivity profiles of networks.
%
%     Params.epsil: 1xL. 
%     The inter-subject concentration parameter, which represents
%     inter-subject functional connectivity variability. A large epsil_l 
%     indicats low inter-subject functional connectivity variability for 
%     network l.
%    
%     Params.s_psi: DxLxS. 
%     The functional connectivity profiles of L networks for S subjects.
%
%     Params.sigma: 1xL. 
%     The intra-subject concentration parameter, which represents
%     intra-subject functional connectivity variability. A large sigma_l
%     indicates low intra-subject functional connectivity variability for
%     network l.
%
%     Params.s_t_nu: DxLxTxS. 
%     The functional connectivity profiles of L networks for S subjects 
%     and each subject has T sessions.
%
%     Params.kappa: 1xL. 
%     The inter-region concentration parameter, which represents
%     inter-region functional connectivity variability. A large kappa_l 
%     indicates low inter-region functional variability for network l. 
%     However, please note in this script, we assume kappa to be the same 
%     across networks.
%
%     Params.s_lambda: NxLxS. 
%     The posterior probability of the individual-specific parcellation of 
%     each subject. 
%
%     Params.theta: NxL. 
%     The spatial prior denotes the probability of networks occurring at 
%     each spatial location.
%
% Example:
%   Params = CBIG_MSHBM_estimate_group_priors(project_dir,'fsaverage5','37','2','17')
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath('../lib/');

%% setting parameters
setting_params.num_sub = str2double(num_sub); 
setting_params.num_session = str2double(num_sess);
setting_params.num_clusters = str2double(num_clusters);
setting_params.epsilon = 1e-4;
setting_params.mesh = mesh;

%% read in data
data = fetch_data(project_dir, setting_params.num_session, setting_params.num_sub, setting_params.mesh);
fprintf('read data...DONE!\n');
setting_params.dim = size(data.series,2) - 1;
if(setting_params.dim < 1200) 
    setting_params.ini_concentration = 500; 
elseif(setting_params.dim >= 1200 && setting_params.dim < 1800 )
    setting_params.ini_concentration = 650;                                          
else
    error('Data dimension is higher than 1800, please set the minimal concentration parameter value ini_concentration, where besseli((D/2-1),ini_concentration) > 0 and relatively small\n');
end    

%% paramter initialization
% load group parcellation and parameters
group = load(fullfile(project_dir, 'group', 'group.mat'));
setting_params.g_mu = transpose(group.clustered.mtc);

% mu: DxL. The group-level functional connectivity profiles of networks
Params.mu = setting_params.g_mu;

% epsil: 1xL. The inter-subject concentration parameter, which represents
% inter-subject functional connectivity variability. A large epsil_l 
% indicats low inter-subject functional connectivity variability for 
% network l.
Params.epsil = setting_params.ini_concentration*ones(1, setting_params.num_clusters);

% s_psi: DxLxS. The functional connectivity profiles of L networks for S
% subjects.
Params.s_psi = repmat(setting_params.g_mu, 1, 1, setting_params.num_sub);

% sigma: 1xL. The intra-subject concentration parameter, which represents
% intra-subject functional connectivity variability. A large sigma_l
% indicates low intra-subject functional connectivity variability for
% network l.
Params.sigma = setting_params.ini_concentration*ones(1, setting_params.num_clusters);

% s_t_nu: DxLxTxS. The functional connectivity profiles of L networks for S
% subjects and each subject has T sessions.
Params.s_t_nu = repmat(setting_params.g_mu, 1, 1, setting_params.num_session, setting_params.num_sub);

% kappa: 1xL. The inter-region concentration parameter, which represents
% inter-region functional connectivity variability. A large kappa_l 
% indicates low inter-region functional variability for network l. However,
% please note in this script, we assume kappa to be the same across 
% networks.
Params.kappa = setting_params.ini_concentration*ones(1, setting_params.num_clusters);

% s_lambda: NxLxS. The posterior probability of the individual-specific
% parcellation of each subject. 
log_vmf = permute(Params.s_t_nu, [1,2,4,3]);
log_vmf = mtimesx(data.series, log_vmf);%NxLxSxT
log_vmf = bsxfun(@times, permute(log_vmf,[2,1,3,4]), transpose(Params.kappa));%LxNxSxT
log_vmf = bsxfun(@plus,Cdln(transpose(Params.kappa), setting_params.dim), log_vmf);%LxNxSxT
log_vmf = sum(log_vmf, 4);
log_vmf = permute(log_vmf, [2,1,3]);
s_lambda = bsxfun(@minus, log_vmf, max(log_vmf,[],2));
mask = repmat((sum(s_lambda,2)==0), 1, setting_params.num_clusters, 1);
s_lambda = exp(s_lambda);
Params.s_lambda = bsxfun(@times, s_lambda, 1./sum(s_lambda,2));
Params.s_lambda(mask) = 0;

%theta: NxL. The spatial prior denotes the probability of networks
%occurring at each spatial location.
Params.theta = mean(Params.s_lambda, 3);
theta_num = sum(Params.s_lambda ~= 0, 3);
Params.theta(Params.theta ~= 0) = Params.theta(Params.theta ~= 0)./theta_num(Params.theta ~= 0);

%% EM
stop_inter = 0;
cost_inter = 0;
Params.iter_inter = 0;
while(stop_inter == 0)
    Params.iter_inter = Params.iter_inter + 1;
    %% Intra subject variability
    cost = 0;
    iter_intra_em = 0;
    stop_intra_em = 0;
    Params.sigma = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
    Params.s_psi = repmat(setting_params.g_mu, 1, 1, setting_params.num_sub);
    while(stop_intra_em == 0)
        iter_intra_em = iter_intra_em + 1;

        fprintf('Inter-region iteration %d ...\n', iter_intra_em);
        Params.kappa = setting_params.ini_concentration*ones(1, setting_params.num_clusters);
        Params.s_t_nu = repmat(setting_params.g_mu, 1, 1, setting_params.num_session, setting_params.num_sub);
        Params = vmf_clustering_subject_session(Params, setting_params, data);

        fprintf('Intra-subject variability level...\n');
        Params = intra_subject_var(Params, setting_params);

        update_cost = bsxfun(@times, Params.s_psi, permute(Params.s_t_nu, [1,2,4,3]));
        update_cost = sum(bsxfun(@times, Params.sigma, update_cost), 1);
        update_cost = bsxfun(@plus, Cdln(Params.sigma, setting_params.dim), update_cost);
        update_cost = sum(sum(sum(update_cost, 2), 3), 4) + sum(sum(bsxfun(@plus, sum(bsxfun(@times, bsxfun(@times, Params.mu, Params.s_psi),Params.epsil), 1), Cdln(Params.epsil, setting_params.dim)), 2), 3);
        update_cost = update_cost + sum(Params.cost_em);
        if(abs(abs(update_cost - cost)./cost) <= setting_params.epsilon)
            stop_intra_em = 1;
            Params.cost_intra = update_cost;
        end
        if(iter_intra_em >= 50)
            stop_intra_em = 1;
            Params.cost_intra = update_cost;
        end
        cost = update_cost;
    end

    %% Inter subject variability
    fprintf('Inter-subject variability level ...\n');
    Params = inter_subject_var(Params, setting_params);

    update_cost_inter = Params.cost_intra;    
    Params.Record(Params.iter_inter) = update_cost_inter;
    if(abs(abs(update_cost_inter - cost_inter)./cost_inter) <= 1e-6 || Params.iter_inter >= 50)
        stop_inter = 1;
        Params.cost_inter = update_cost_inter;

        % set s_lambda, s_psi, s_t_nu to be empty to save space and time
        Params.s_lambda = [];
        Params.s_psi = [];
        Params.s_t_nu = [];
        save(fullfile(project_dir, 'priors', 'Params_Final.mat'), 'Params');
    end
    cost_inter = update_cost_inter;
    mkdir(fullfile(project_dir, 'priors'));
    save(fullfile(project_dir, 'priors', ['Params_iteration',num2str(Params.iter_inter),'.mat']), 'Params');
end

rmpath('../lib/');

end


%% sub-functions

function Params=inter_subject_var(Params,setting_params)

% Inter-subject variability level

% update mu
mu_update = sum(Params.s_psi, 3);
mu_update = bsxfun(@times, mu_update, 1./sqrt(sum((mu_update).^2)));

% update epsil
epsil_update = bsxfun(@times, Params.s_psi, mu_update);
epsil_update = sum(sum(epsil_update, 1), 3);
epsil_update = epsil_update./setting_params.num_sub;
for i = 1:setting_params.num_clusters
    epsil_update(i) = invAd(setting_params.dim, epsil_update(i));
    if(epsil_update(i) < setting_params.ini_concentration)
        epsil_update(i) = setting_params.ini_concentration;
        fprintf('[WARNING] epsil of %d is less than the minimal value %d \n',i, setting_params.ini_concentration);
    end
    if(isinf(epsil_update(i)))
        epsil_update(i) = Params.epsil(i);
        fprintf('[WARNING] epsil of %d is Inf\n',i);
    end
end
Params.mu = mu_update;
Params.epsil = epsil_update;
end

function Params = intra_subject_var(Params,setting_params)

% Intra-subject variability level

flag_psi = zeros(setting_params.num_sub, 1);
stop_intra = 0;
iter_intra = 0;
while(stop_intra == 0)
    iter_intra = iter_intra + 1;
    fprintf('It is inter interation %d intra iteration %d..update s_psi and sigma..\n',Params.iter_inter,iter_intra);
    % update s_psi
    s_psi_update = sum(bsxfun(@times,Params.s_t_nu,repmat(Params.sigma,size(Params.s_t_nu,1),1,size(Params.s_t_nu,3),size(Params.s_t_nu,4))),3);
    s_psi_update = reshape(s_psi_update,size(s_psi_update,1),size(s_psi_update,2),size(s_psi_update,3)*size(s_psi_update,4));
    s_psi_update = bsxfun(@plus,s_psi_update,bsxfun(@times,Params.epsil,Params.mu));
    s_psi_update = bsxfun(@times,s_psi_update,1./sqrt(sum((s_psi_update).^2)));
   
    for s = 1:setting_params.num_sub
        checkpsi = diag(s_psi_update(:,:,s)'*Params.s_psi(:,:,s));
        checkpsi_flag = (sum(1-checkpsi < setting_params.epsilon) < setting_params.num_clusters);
        if(checkpsi_flag < 1)
            flag_psi(s,1) = 1;
        end
    end
    Params.s_psi = s_psi_update;
    
    % update sigma
    sigma_update = bsxfun(@times,Params.s_psi,permute(Params.s_t_nu,[1,2,4,3]));
    sigma_update = sum(sum(sum(sigma_update,1),3),4)./(setting_params.num_sub*setting_params.num_session);%1xLxSxT=>1xL
    for i = 1:setting_params.num_clusters
        sigma_update(i) = invAd(setting_params.dim,sigma_update(i));
    end
    
    if((sum(flag_psi) == setting_params.num_sub) && (mean(abs(Params.sigma-sigma_update)./Params.sigma) < setting_params.epsilon))
        stop_intra = 1;
    end
    Params.sigma = sigma_update;
end
end

function Params = vmf_clustering_subject_session(Params,setting_params,data)

% Inter-region level

stop_em = 0;
iter_em = 0;
cost = zeros(1, setting_params.num_sub);

while(stop_em == 0)
    iter_em = iter_em + 1;
    fprintf('It is EM iteration.. %d..\n',iter_em);
    %% Mstep
    flag_nu = zeros(setting_params.num_sub,setting_params.num_session);
    stop_m = 0;
    iter_m = 0;
    fprintf('M-step..\n');
    while(stop_m == 0)
        iter_m = iter_m + 1;
        
        % update kappa
        s_lambda = Params.s_lambda;%NxLxS
        kappa_update = mtimesx(data.series,permute(Params.s_t_nu,[1,2,4,3]));
        kappa_update = bsxfun(@times,s_lambda,kappa_update);
        kappa_update = sum(sum(sum(sum(kappa_update,1),4),3));
        kappa_update = kappa_update./sum((setting_params.num_session.*sum(sum(s_lambda,1),3)));
        kappa_update = invAd(setting_params.dim,kappa_update);
        kappa_update = repmat(kappa_update,1,setting_params.num_clusters);

        if (sum(kappa_update == Inf) ~= 0)
            fprintf('[WARNING] kappa is Inf !\n')
            kappa_update(kappa_update == Inf) = Params.kappa(kappa_update == Inf);
        end
        if (sum(kappa_update < setting_params.ini_concentration) ~= 0)
            kappa_update(kappa_update < setting_params.ini_concentration) = setting_params.ini_concentration;
            fprintf('[WARNING] kappa is less than the minimal value %d \n', setting_params.ini_concentration);
        end
       
        for s = 1:setting_params.num_sub
            for t = 1:setting_params.num_session
                
                % update s_t_nu
                checknu = [];
                X = data.series(:,:,s,t);
                s_lambda = Params.s_lambda(:,:,s);
                lambda_X = bsxfun(@times,kappa_update,X'*s_lambda) + bsxfun(@times,Params.sigma,Params.s_psi(:,:,s));
                s_t_nu_update = bsxfun(@times,lambda_X,1./sqrt(sum((lambda_X).^2)));
                checknu = diag(s_t_nu_update'*Params.s_t_nu(:,:,t,s));
                checknu_flag = (sum(1-checknu < setting_params.epsilon) < setting_params.num_clusters);
                Params.s_t_nu(:,:,t,s) = s_t_nu_update;

                if(checknu_flag < 1)
                    flag_nu(s,t) = 1;
                end
                
            end
        end
        if((sum(sum(flag_nu)) == setting_params.num_sub*setting_params.num_session)&&(mean(abs(Params.kappa-kappa_update)./Params.kappa)<setting_params.epsilon))
            stop_m=1;
        end
        Params.kappa = kappa_update;
    end
    %% Estep
    fprintf('Estep..\n');
    
    % estimate s_lambda
    log_vmf = permute(Params.s_t_nu,[1,2,4,3]);
    log_vmf = mtimesx(data.series,log_vmf);%NxLxSxT
    log_vmf = bsxfun(@times,permute(log_vmf,[2,1,3,4]),transpose(Params.kappa));%LxNxSxT
    log_vmf(:,sum(log_vmf==0,1)==0) = bsxfun(@plus,Cdln(transpose(Params.kappa),setting_params.dim),log_vmf(:,sum(log_vmf==0,1)==0));%LxNxSxT
    log_vmf = sum(log_vmf,4);%NxLxS
    idx = sum(log_vmf == 0,1) ~= 0;
    
    log_vmf = bsxfun(@plus,permute(log_vmf,[2,1,3]),log(Params.theta));
    s_lambda = bsxfun(@minus,log_vmf,max(log_vmf,[],2));
    s_lambda = exp(s_lambda);
    Params.s_lambda = bsxfun(@times,s_lambda,1./sum(s_lambda,2));
    Params.s_lambda = permute(Params.s_lambda,[2,1,3]);
    Params.s_lambda(:,idx) = 0;
    Params.s_lambda = permute(Params.s_lambda,[2,1,3]);
    
    % estimate theta
    Params.theta = mean(Params.s_lambda,3);

    %% em stop criteria
    for s = 1:setting_params.num_sub    
        for t = 1:setting_params.num_session
            X = data.series(:,:,s,t);
            setting_params.num_verts = size(X,1);
            log_vmf = vmf_probability(X,Params.s_t_nu(:,:,t,s), Params.kappa);
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

        update_cost(:,s) = sum(sum(s_lambda_cost.*log_lambda_prop))+sum(sum(s_lambda_cost.*log_theta_cost))-sum(sum(s_lambda_cost.*log_s_lambda_cost));    
    end
    sub_set=find((abs(abs(update_cost-cost)./cost) > setting_params.epsilon) == 0);
    if(length(sub_set) == setting_params.num_sub)
        stop_em = 1;
        Params.cost_em = cost;
    end
    if(iter_em > 100)
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

function data = fetch_data(project_dir,num_session,num_sub,mesh)

% read in input functional connectivity profiles
if(~isempty(strfind(mesh,'fs_LR_32k')))
    load('../lib/fs_LR_32k_medial_mask.mat');
    for t = 1:num_session
        data_profile = fullfile(project_dir,'profile_list','training_set',['sess' num2str(t) '.txt']);
        profile_name = table2cell(readtable(data_profile,'Delimiter',' ','ReadVariableNames',false));
        for i = 1:num_sub
            fprintf('Session %d...It is subj %d...\n',t,i);
            avg_file = profile_name{i,1};

            [~, series, ~] = CBIG_MSHBM_read_fmri(avg_file);       
            series(~medial_mask, :) = 0;

            series = bsxfun(@minus,series,mean(series, 2));
            series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), sqrt(sum(series(all(series,2)~=0,:).^2,2)));
            data.series(:,:,i,t) = series;
        end
    end

elseif(~isempty(strfind(mesh,'fsaverage')))
    lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh, 'inflated','cortex');
    rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh, 'inflated', 'cortex'); 
    for t = 1:num_session
        lh_data_profile = fullfile(project_dir,'profile_list','training_set',['lh_sess' num2str(t) '.txt']);
        rh_data_profile = fullfile(project_dir,'profile_list','training_set',['rh_sess' num2str(t) '.txt']);
        lh_profile_name = table2cell(readtable(lh_data_profile,'Delimiter',' ','ReadVariableNames',false));
        rh_profile_name = table2cell(readtable(rh_data_profile,'Delimiter',' ','ReadVariableNames',false));
        for i = 1:num_sub
            fprintf('Session %d...It is subj %d...\n',t,i);
            lh_avg_file = lh_profile_name{i,1};
            rh_avg_file = rh_profile_name{i,1};

            [~, lh_series, ~] = CBIG_MSHBM_read_fmri(lh_avg_file);
            [~, rh_series, ~] = CBIG_MSHBM_read_fmri(rh_avg_file);

            lh_series(lh_avg_mesh.MARS_label == 1,:) = 0;
            rh_series(rh_avg_mesh.MARS_label == 1,:) = 0;

            series = [lh_series;rh_series];

            series = bsxfun(@minus,series,mean(series, 2));
            series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), sqrt(sum(series(all(series,2)~=0,:).^2,2)));
            data.series(:,:,i,t) = series;
        end
    end
    
end 
    
end

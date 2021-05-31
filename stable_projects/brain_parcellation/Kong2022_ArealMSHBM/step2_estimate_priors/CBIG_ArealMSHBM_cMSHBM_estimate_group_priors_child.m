function CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_child(s, mesh, tmp_dir)

% CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_child(s, mesh, tmp_dir)
%
% This is the child script to estimate group priors for the strictly contiguous individual
% areal-level parcellation generation. To speed things up, we parallize the script
% as parent and child scripts with multiple CPUs. This training step requires multiple
% CPUs to run. For num_sub training subjects, user should submit 1 parent job and
% num_sub child jobs.
% Check CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_parent.m for more details.
%
% Input:
%
%   - s: (string)
%
%     The subject id. If num_sub subjects is specified in the parent script. User should
%     run num_sub child scripts for each subject. For example, '1','2','3'.
%
%   - mesh: (string)
%     
%     The data surface space. 'fsaverage5/fsaverage6/fsaverage' or 
%     'fs_LR_32k'. 
%
%   - tmp_dir: (string)
%
%     The directory to save temporary results. This should be defined as the same 
%     temporary directory as the parent script. 
%
% Example:
%   CBIG_ArealMSHBM_cMSHBM_estimate_group_priors_child('1', 'fs_LR_32k', tmp_dir)
%
% Written by Ru(by) Kong and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));
addpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2022_ArealMSHBM', 'lib'));

s = str2num(s);

% read in mesh & parcellation
if(~isempty(strfind(mesh, 'fsaverage')))
    lh_avg_mesh=CBIG_ReadNCAvgMesh('lh', mesh, 'sphere','cortex');
    rh_avg_mesh=CBIG_ReadNCAvgMesh('rh', mesh, 'sphere','cortex');
elseif(~isempty(strfind(mesh, 'fs_LR')))
    lh_avg_mesh = CBIG_read_fslr_surface('lh', mesh, 'sphere','medialwall.annot');
    rh_avg_mesh = CBIG_read_fslr_surface('rh', mesh, 'sphere','medialwall.annot');
end

spatial_data = transpose([lh_avg_mesh.vertices rh_avg_mesh.vertices]);

spatial_data = bsxfun(@rdivide, spatial_data, sqrt(sum(spatial_data.^2,2)));

while(~exist([tmp_dir '/all_sub_profile_list.mat']))
    pause(1);
end
load([tmp_dir '/all_sub_profile_list.mat']);

while(~exist([tmp_dir '/setting_params.mat']))
    pause(1);
end
load([tmp_dir '/setting_params.mat']);

%% read in profiles for each session
if(~isempty(strfind(mesh,'fs_LR_32k')))
    medial_mask = [lh_avg_mesh.MARS_label == 1;rh_avg_mesh.MARS_label == 1];
    for t = 1:setting_params.num_session
        avg_file = profile_list{s,t};

        [~, series, ~] = CBIG_MSHBM_read_fmri(avg_file);       
        series(medial_mask, :) = 0;

        series = bsxfun(@minus,series,mean(series, 2));
        series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), ...
                                        sqrt(sum(series(all(series,2)~=0,:).^2,2)));
        data_series(:,:,t) = series;
    end

elseif(~isempty(strfind(mesh,'fsaverage')))
    for t = 1:setting_params.num_session
        lh_avg_file = lh_profile_list{s,t};
        rh_avg_file = rh_profile_list{s,t};

        [~, lh_series, ~] = CBIG_MSHBM_read_fmri(lh_avg_file);
        [~, rh_series, ~] = CBIG_MSHBM_read_fmri(rh_avg_file);

        if(length(size(lh_series)) > 2)
            lh_series = reshape(lh_series, size(lh_avg_mesh.vertices,2), size(lh_series,4));
            rh_series = reshape(rh_series, size(rh_avg_mesh.vertices,2), size(rh_series,4));
        end
        
        lh_series(lh_avg_mesh.MARS_label == 1,:) = 0;
        rh_series(rh_avg_mesh.MARS_label == 1,:) = 0;

        series = [lh_series;rh_series];

        series = bsxfun(@minus,series,mean(series, 2));
        series(all(series,2)~=0,:) = bsxfun(@rdivide,series(all(series,2)~=0,:), ...
                                        sqrt(sum(series(all(series,2)~=0,:).^2,2)));
        data_series(:,:,t) = series;
    
    end
    
end 

%% parameter initialization
fprintf('Start initialization!\n')
group = load(setting_params.group);
if(setting_params.dim == 1482 )
    ini_val = 650;
else   
    ini_val = CBIG_ArealMSHBM_initialize_concentration(setting_params.dim);
end
ini_kappa=ini_val*ones(1,setting_params.num_clusters);
g_mu=group.mtc;

for t = 1:setting_params.num_session
    ini_log_vmf_tmp =g_mu;
    ini_log_vmf_tmp = data_series(:,:,t)*ini_log_vmf_tmp;
    ini_log_vmf_tmp = bsxfun(@times,ini_log_vmf_tmp',ini_kappa');
    ini_log_vmf_tmp = transpose(bsxfun(@plus,Cdln(transpose(ini_kappa),setting_params.dim),ini_log_vmf_tmp));%LxNxSxT
    if(t == 1)
        ini_log_vmf = ini_log_vmf_tmp;
    else
        ini_log_vmf = ini_log_vmf +ini_log_vmf_tmp;
    end
end
ini_s_lambda=bsxfun(@minus,ini_log_vmf,max(ini_log_vmf,[],2));
mask=repmat((sum(ini_s_lambda,2)==0),1,setting_params.num_clusters,1);
ini_s_lambda=exp(ini_s_lambda);

if(~exist([tmp_dir '/ini']))
    mkdir([tmp_dir '/ini']);
end
save([tmp_dir '/ini/sub' num2str(s) '.mat'],'ini_s_lambda','mask');
save_dump([tmp_dir '/ini/sub' num2str(s) '.mat']);

%% EM

iter_inter_em = 0;
while(1)
    fprintf('Start inter iteration!\n');

    iter_inter_em = iter_inter_em + 1;
    iter_intra_em = 0;
    while(1)
        fprintf('Start intra iteration!\n');

        iter_intra_em = iter_intra_em + 1;
        iter_em = 0;
        while(1)
            fprintf('Start EM iteration!\n');
            iter_em = iter_em + 1;
            iter_m = 0;
            while(1)
                fprintf('compute kappa step1 ...\n')
                iter_m = iter_m + 1;
                Mstep_params_tmp = load_until_exist([tmp_dir '/Mstep/setting_params/sub'...
                     num2str(s) '_' num2str(iter_m) '.mat']);

                %save out s_lambda.*(X*s_t_nu)
                kappa_update_tmp=mtimesx(data_series,Mstep_params_tmp.tmp_s_t_nu);
                kappa_update_tmp=bsxfun(@times,Mstep_params_tmp.tmp_s_lambda,kappa_update_tmp);
                kappa_update_tmp=sum(sum(sum(kappa_update_tmp)));

                if(~exist([tmp_dir '/Mstep/kappa_step1']))
                    mkdir([tmp_dir '/Mstep/kappa_step1']);
                end
                save([tmp_dir '/Mstep/kappa_step1/sub' num2str(s) '_' num2str(iter_m) '.mat'],'kappa_update_tmp');
                save_dump([tmp_dir '/Mstep/kappa_step1/sub' num2str(s) '_' num2str(iter_m) '.mat']);

                % load in kappa temporary files...
                kappa_tmp = load_until_exist([tmp_dir '/Mstep/kappa_step2/kappa_' num2str(iter_m) '.mat']); 

                for t=1:setting_params.num_session
                    checknu=[];
                    X=data_series(:,:,t);    
                    lambda_X=bsxfun(@times,kappa_tmp.kappa_update,X'*Mstep_params_tmp.tmp_s_lambda)+...
                        bsxfun(@times,Mstep_params_tmp.tmp_sigma,Mstep_params_tmp.tmp_psi);
                    s_t_nu_update=bsxfun(@rdivide,lambda_X,sqrt(sum((lambda_X).^2)));
                    checknu= diag(s_t_nu_update'*Mstep_params_tmp.tmp_s_t_nu(:,:,t));
                    checknu_flag= (sum(1-checknu < setting_params.epsilon) < setting_params.num_clusters);
                    tmp_s_t_nu(:,:,t)=s_t_nu_update;

                    if(checknu_flag <1)
                        flag_nu(t) = 1;
                    else
                        flag_nu(t) = 0;
                    end

                end

                if(~exist([tmp_dir '/Mstep/s_t_nu']))
                    mkdir([tmp_dir '/Mstep/s_t_nu']);
                end
                save([tmp_dir '/Mstep/s_t_nu/sub' num2str(s) '_' num2str(iter_m) '.mat'],'tmp_s_t_nu','flag_nu');
                save_dump([tmp_dir '/Mstep/s_t_nu/sub' num2str(s) '_' num2str(iter_m) '.mat']);

                Mstep_iter = load_until_exist([tmp_dir '/Mstep/converge_flag/iter_' num2str(iter_m) '.mat']);

                if(Mstep_iter.stop_m == 1)
                    break;
                end

            end

            fprintf('Mstep done ...\n')

            Mstep_converge_flag = load([tmp_dir '/Mstep/converge_flag/Mstep_converge.mat']);
            if(~isfield(Mstep_converge_flag,'stop_m'))
                error('stop_m does not exist')
            end
            if(Mstep_converge_flag.stop_m == 0)
                error('stop_m = 0')
            end

            fprintf('Spatial connectness prior ... \n');
            
            lambda_X_coonnect = spatial_data'*Mstep_params_tmp.tmp_s_lambda;
            tmp_s_muc = bsxfun(@rdivide,lambda_X_coonnect,sqrt(sum((lambda_X_coonnect).^2)));
            tmp_spatial_connect_vmf = spatial_data*tmp_s_muc; % 64984x3 * 3x400
            tmp_spatial_connect_vmf(1:size(tmp_spatial_connect_vmf,1)/2, setting_params.num_clusters/2+1:end) = 0;
            tmp_spatial_connect_vmf(size(tmp_spatial_connect_vmf,1)/2+1:end, 1:setting_params.num_clusters/2) = 0;
            
            tmp_gamma = nansum(sum(bsxfun(@times, Mstep_params_tmp.tmp_s_lambda, tmp_spatial_connect_vmf)));
            
            if(~exist([tmp_dir '/Mstep/s_muc']))
                mkdir([tmp_dir '/Mstep/s_muc']);
            end
            save([tmp_dir '/Mstep/s_muc/sub' num2str(s) '.mat'],'tmp_s_muc','tmp_gamma');
            save_dump([tmp_dir '/Mstep/s_muc/sub' num2str(s) '.mat']);

            if(~exist([tmp_dir '/Estep/converge_flag/Estep_cost.mat.dump']))
                fprintf('compute s_lambda step1 ...\n')

                Estep_params_tmp = load_until_exist([tmp_dir '/Estep/setting_params/sub' num2str(s) '.mat']);
                log_vmf = [];
                for t = 1:setting_params.num_session
                    setting_params.num_verts=size(data_series,1);

                    log_vmf_tmp = Estep_params_tmp.tmp_s_t_nu(:,:,t);
                    log_vmf_tmp = data_series(:,:,t)*log_vmf_tmp;
                    log_vmf_tmp = bsxfun(@times,log_vmf_tmp',Estep_params_tmp.tmp_kappa');%LxNxSxT
                    log_vmf_tmp(:,sum(log_vmf_tmp==0,1)==0)=bsxfun(@plus,Cdln(transpose(Estep_params_tmp.tmp_kappa),...
                        setting_params.dim),log_vmf_tmp(:,sum(log_vmf_tmp==0,1)==0));
                    log_vmf(:,:,t) = log_vmf_tmp;
                end


                log_vmf=transpose(sum(log_vmf,3));
                tmp_log_lambda_prop = log_vmf;
                tmp_idx = repmat(sum(log_vmf==0,2)~=0,1,setting_params.num_clusters);

                log_vmf=bsxfun(@plus,log_vmf,log(Estep_params_tmp.tmp_theta));
                
                %% add xyz connectness prior
                gamma_and_connect = load_until_exist([tmp_dir '/Mstep/gamma/gamma.mat']);
                
                log_connect_vmf = tmp_spatial_connect_vmf;    
                log_connect_vmf = log_connect_vmf.*gamma_and_connect.gamma_update;%NxL
                log_connect_vmf = log_connect_vmf +  Cdln(gamma_and_connect.gamma_update,3);
                
                log_vmf = log_vmf + setting_params.beta.*log_connect_vmf;
                
                
                tmp_s_lambda=bsxfun(@minus,log_vmf,max(log_vmf,[],2));
                tmp_s_lambda=exp(tmp_s_lambda);

                if(~exist([tmp_dir '/Estep/s_lambda_step1']))
                    mkdir([tmp_dir '/Estep/s_lambda_step1']);
                end
                save([tmp_dir '/Estep/s_lambda_step1/sub' num2str(s) '.mat'], 'tmp_s_lambda','tmp_idx');
                save_dump([tmp_dir '/Estep/s_lambda_step1/sub' num2str(s) '.mat']);

                s_lambda_step2 = load_until_exist([tmp_dir '/Estep/s_lambda_step2/sub' num2str(s) '.mat']);
                theta_estep = load_until_exist([tmp_dir '/Estep/theta/theta.mat']);

                theta_cost=theta_estep.tmp_theta;
                log_theta_cost = log(theta_cost);
                log_theta_cost(isinf(log_theta_cost)) = log(eps.^20);

                s_lambda_cost = s_lambda_step2.tmp_lambda;
                log_s_lambda_cost = log(s_lambda_cost);
                log_s_lambda_cost(isinf(log_s_lambda_cost)) = log(eps.^20);
                
                log_connect_vmf(isnan(log_connect_vmf)) = log(eps.^20);
                log_connect_vmf(isinf(log_connect_vmf)) = log(eps.^20);

                cost = sum(sum(s_lambda_cost.*tmp_log_lambda_prop))+sum(sum(s_lambda_cost.*log_theta_cost))-...
                    sum(sum(s_lambda_cost.*log_s_lambda_cost))+...
                    setting_params.beta.* sum(sum(s_lambda_cost.*log_connect_vmf));    

                if(~exist([tmp_dir '/Estep/cost']))
                    mkdir([tmp_dir '/Estep/cost']);
                end
                save([tmp_dir '/Estep/cost/sub' num2str(s) '.mat'],'cost');
                save_dump([tmp_dir '/Estep/cost/sub' num2str(s) '.mat']);

            end

            EMstep_iter = load_until_exist([tmp_dir '/inter_' num2str(iter_inter_em) '/intra_'...
                num2str(iter_intra_em) '/EMstep_iter' num2str(iter_em) '/converge_flag/iter_' num2str(iter_em) '.mat']);

            if(EMstep_iter.stop_em == 1)
                break;
            end
        end

        intra_iter = load_until_exist([tmp_dir '/inter_' num2str(iter_inter_em) '/intra_'...
             num2str(iter_intra_em) '/converge_flag/iter_' num2str(iter_intra_em) '.mat']);
        if(intra_iter.stop_intra_em == 1)
            break;
        end
    end

    inter_iter = load_until_exist([tmp_dir '/inter_' num2str(iter_inter_em) '/converge_flag/iter_'...
        num2str(iter_inter_em) '.mat']);

    if(inter_iter.stop_inter == 1)
        fprintf('Finished!\n')
        break;
    end
end

rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2019_MSHBM', 'lib'));
rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', 'Kong2022_ArealMSHBM', 'lib'));

end

function out = Ad(in,D)
out = besseli(D/2,in) ./ besseli(D/2-1,in);
end

function out = Cdln(k,d,k0)

% Computes the logarithm of the partition function of vonMises-Fisher as
% a function of kappa

k = double(k);

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

function out_var = load_until_exist(filename)

while(~exist([filename '.dump']))
    pause(1);
end
out_var = load(filename);

end
    
function save_dump(filename)

fclose(fopen([filename '.dump'],'w'));

end

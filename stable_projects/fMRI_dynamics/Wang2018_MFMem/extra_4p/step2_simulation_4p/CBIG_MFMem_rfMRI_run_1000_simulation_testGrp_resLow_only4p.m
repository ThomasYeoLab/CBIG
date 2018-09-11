function CBIG_MFMem_rfMRI_run_1000_simulation_testGrp_resLow_only4p(save_dir_input)

%--------------------------------------------------------------------------
% CBIG_mfm_rfMRI_run_1000_simulation_testGrp_resLow
%
% Input parameters:
%       save_dir_input: Output file directory
%
% This function will run 1000 simluations of FC with estimated parameter and
% SC of test group, simulation will run in low time resolution t = 0.01s
%
% Written by Peng Wang and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md
%--------------------------------------------------------------------------

%% setup directory for lib, data, save
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 
save_dir = save_dir_input;

cd('..')
cd('..')
high_dir = pwd;
lib_dir = fullfile(high_dir, 'lib');
addpath(lib_dir);




%% load model parameter
parameter_file_name = 'Estimated_Parameter_for_only_4p_2018821.mat';
load ([data_dir '/' parameter_file_name],'Para_E');


%% load FC, SC
FCSC_file_name = 'FCSC_Desikan68_Raphael_Wang.mat';
load([data_dir '/' FCSC_file_name]);

SC = SC_test;
FC = FC_test;


%scaling the SC
SC = SC./max(max(SC)).*0.2;

%prepare FC: use the entries above main diagonal
FC_mask = tril(ones(size(FC,1),size(FC,1)),0);
y = FC(~FC_mask);




%% begin simulation


funcP = @(Para_E,Nstate) CBIG_MFMem_rfMRI_nsolver_eul_sto_resLH_4p(Para_E,SC,y,FC_mask,Nstate,14.4,0.72,0);
    
numSimulation = 1000;

for i = 1:numSimulation

disp(['Sim:' num2str(i)]);

Nstate = rng;
[h_output(i,:), CC_check(i)] = funcP(Para_E,Nstate);

CC_check(i)

end

%% plot result
set(figure,'Position',[100 130 600 400],'Color','w')

plot([1:numSimulation],CC_check,'ko-','markerfacecolor','k')
xlabel('simulation number','FontSize',9)
ylabel('Similarity','FontSize',9)

%% save result
saved_date = fix(clock);

save( [save_dir '/1000simulation_testGrp_resLow_' num2str(saved_date(1)) num2str(saved_date(2)) num2str(saved_date(3))],'CC_check');
rmpath(lib_dir);

end

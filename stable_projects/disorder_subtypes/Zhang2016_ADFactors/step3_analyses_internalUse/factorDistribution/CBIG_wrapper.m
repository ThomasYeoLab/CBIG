clear;
clc;
close all;

addpath(genpath('../functions/'));
CBIG_plotSetup;

SUBINFO_FILE = '../../subInfo/rid_sex_mDoB_yDoB_mExam_yExam_ageExam_dx_icv_bl.csv';
GMICV_FILE = '../../vbm_810bl/gmVolAndICV/rid_gmVol_icv_bl.csv';

%!!!!!!
K = 3;

% AD
GRP = 'ad_188'; % possible values: ad_188 mci_394 a+_mci_147 converter_mci_198 a+_hc_43
CBIG_visualize_factorComp_mixedAmyloid(SUBINFO_FILE, GMICV_FILE, K, GRP);
title('188 AD Dementia Patients');
% a+ MCI
GRP = 'a+_mci_147'; % possible values: ad_188 mci_394 a+_mci_147 converter_mci_198 a+_hc_43
CBIG_visualize_factorComp(SUBINFO_FILE, GMICV_FILE, K, GRP);
title('147 Amyloid+ MCI Participants');
% a+ CN
GRP = 'a+_hc_43'; % possible values: ad_188 mci_394 a+_mci_147 converter_mci_198 a+_hc_43
CBIG_visualize_factorComp(SUBINFO_FILE, GMICV_FILE, K, GRP);
title('43 Amyloid+ CN Participants');

% Save all figures
h = get(0, 'children');
for idx = 1:numel(h)
    saveas(h(idx), ['k=' num2str(K) '_fig' num2str(idx) '.pdf']);
end
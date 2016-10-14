function rid = CBIG_select_predementedGrp(rid_dx, grpName)

% rid = CBIG_select_predementedGrp(rid_dx, grpName)
% Possible values: ad_188 mci_394 a+_mci_147 converter_mci_198 a+_hc_43 a+_hc&mci_190

%
% Written by CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

switch grpName
    case 'ad_188'
        rid = rid_dx(rid_dx(:, 2)==3, 1);
        fprintf('\nCohort: ad_188 (N = %i)\n', numel(rid));
    case 'mci_394'
        rid = rid_dx(rid_dx(:, 2)==2, 1);
        fprintf('\nCohort: mci_394 (N = %i)\n', numel(rid));
    case 'converter_mci_198'
        mci_rid = rid_dx(rid_dx(:, 2)==2, 1);
        rid_date_age = get_convDate_convAge(mci_rid);
        rid = cell2mat(rid_date_age(:, 1));
        fprintf('\nCohort: converter_mci_198 (N = %i)\n', numel(rid));
    case 'a+_mci_147'
        mci_rid = rid_dx(rid_dx(:, 2)==2, 1);
        amyloid = CBIG_get_amyloid(mci_rid);
        rid_ab = cell2mat(amyloid(:, [1 4]));
        rid = rid_ab(rid_ab(:, 2)<192, 1);
        fprintf('\nCohort: a+_mci_147 (N = %i)\n', numel(rid));
    case 'a+_hc_43'
        hc_rid = rid_dx(rid_dx(:, 2)==1, 1);
        amyloid = CBIG_get_amyloid(hc_rid);
        rid_ab = cell2mat(amyloid(:, [1 4]));
        rid = rid_ab(rid_ab(:, 2)<192, 1);
        fprintf('\nCohort: a+_hc_43 (N = %i)\n', numel(rid));
    case 'a+_hc&mci_190'
        hc_mci_rid = rid_dx(rid_dx(:, 2)<3, 1);
        amyloid = CBIG_get_amyloid(hc_mci_rid);
        rid_ab = cell2mat(amyloid(:, [1 4]));
        rid = rid_ab(rid_ab(:, 2)<192, 1);
        fprintf('\nCohort: a+_hc&mci_190 (N = %i)\n', numel(rid));
    otherwise
        error('Such group not configured!');
end

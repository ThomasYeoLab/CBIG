function CBIG_MMLDA_get_subinfo_wrapper(out_dir)
% Get subject information 
% Written by Nanbo Sun and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

if nargin < 1 || isempty(out_dir)
    out_dir = [getenv('CBIG_REPDATA_DIR') '/stable_projects/disorder_subtypes' ...
    '/Sun2019_ADJointFactors/step2_MMLDA/data'];
end

CBIG_CODE_DIR = getenv('CBIG_CODE_DIR');
addpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

% get the descript of each behavior score and it's the same for both ADNI1
% and ADNI2, so we just need to use either one
[ADAS, MMSE, NEUROBAT] = CBIG_MMLDA_behavior_choose('ADNI2');
ADNI_DOC_DIR = getenv('CBIG_MMLDA_ANDI_DOC_DIR');
ADNI_spreadsheet_path = [ADNI_DOC_DIR '/All/ADNI_161017/documentation'];
ADNI_VBM_path = [getenv('CBIG_MMLDA_ADNI_DIR') '/Sun2019_SPM_VBM'];

%%%
% ADNI2 bl
%%%
% get rid
rid_file = [ADNI_VBM_path '/preprocessing/lists/ADNI2_bl/ADNI2_bl_rid_rm_seg_fail.txt'];
rid = load(rid_file);
rid_cell = CBIG_MMLDA_matrix2cellstr(rid);
phase_cell = repmat({'ADNI2'}, length(rid_cell), 1);

% get age gender dx by combining 'sc' and 'bl'
DXSUM_file = [ADNI_spreadsheet_path '/Assessments/DXSUM_PDXCONV_ADNIALL.csv'];
PTDEMOG_file = [ADNI_spreadsheet_path '/Subject_Characteristics/PTDEMOG.csv'];

viscode_cell = repmat({'sc'}, length(rid_cell), 1);
[sc_age, sc_gender, sc_dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);
viscode_cell = repmat({'bl'}, length(rid_cell), 1);
[bl_age, bl_gender, bl_dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);
age = sc_age;
age(isnan(sc_age)) = bl_age(isnan(sc_age));
gender = sc_gender;
gender(isnan(sc_gender)) = bl_gender(isnan(sc_gender));
dx = sc_dx;
dx(isnan(sc_dx)) = bl_dx(isnan(sc_dx));

% get icv gmVol
gmVolAndICV_file = [ADNI_VBM_path '/preprocessing/output/ADNI2_bl/gmVolAndICV/id_gmVol_icv.mat']; 
load(gmVolAndICV_file)
icv = cell2mat(id_gmVol_icv(:, 3));
gmVol = cell2mat(id_gmVol_icv(:, 2));

% get behavioral scores by combing 'sc' and 'bl'
ADAS_ADNI1_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADASSCORES.csv'];
ADAS_ADNI2_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADAS_ADNIGO2.csv'];
MMSE_spreadsheet = [ADNI_spreadsheet_path '/Assessments/MMSE.csv'];
NEUROBAT_spreadsheet = [ADNI_spreadsheet_path '/Assessments/NEUROBAT.csv'];

viscode_cell = repmat({'bl'}, length(rid_cell), 1);
bl_behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, ... 
    NEUROBAT_spreadsheet, 'ADNI2', repmat({'ADNI2'}, length(rid), 1), rid_cell, viscode_cell);
viscode_cell = repmat({'sc'}, length(rid_cell), 1);
sc_behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, ...
    NEUROBAT_spreadsheet, 'ADNI2', repmat({'ADNI2'}, length(rid), 1), rid_cell, viscode_cell);
bl_sc_behavior_scores = bl_behavior_scores;
bl_sc_behavior_scores(isnan(bl_behavior_scores)) = sc_behavior_scores(isnan(bl_behavior_scores));

% create cellarray and write it into csv file
subinfo = [rid gender age dx icv gmVol bl_sc_behavior_scores];
subinfo_cell = CBIG_MMLDA_matrix2cellstr(subinfo);
var_names = [{'RID', 'GENDER', 'AGE', 'DX', 'ICV', 'gmVol'} ADAS.DESCRIPTION' MMSE.DESCRIPTION' NEUROBAT.DESCRIPTION'];
cell2csv([out_dir '/ADNI2_bl_subinfo.csv'], [var_names; subinfo_cell] , ',')

%%%
% ADNI2 m12
%%%
% get rid
rid_file = [ADNI_VBM_path '/preprocessing/lists/ADNI2_m12/ADNI2_m12_rid_rm_seg_fail.txt'];
rid = load(rid_file);
rid_cell = CBIG_MMLDA_matrix2cellstr(rid);
phase_cell = repmat({'ADNI2'}, length(rid_cell), 1);

% get age gender dx 
DXSUM_file = [ADNI_spreadsheet_path '/Assessments/DXSUM_PDXCONV_ADNIALL.csv'];
PTDEMOG_file = [ADNI_spreadsheet_path '/Subject_Characteristics/PTDEMOG.csv'];
viscode_cell = repmat({'m12'}, length(rid_cell), 1);
[age, gender, dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);

% get icv gmVol
gmVolAndICV_file = [ADNI_VBM_path '/preprocessing/output/ADNI2_m12/gmVolAndICV/id_gmVol_icv.mat']; 
load(gmVolAndICV_file)
icv = cell2mat(id_gmVol_icv(:, 3));
gmVol = cell2mat(id_gmVol_icv(:, 2));

% get behavioral scores 
ADAS_ADNI1_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADASSCORES.csv'];
ADAS_ADNI2_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADAS_ADNIGO2.csv'];
MMSE_spreadsheet = [ADNI_spreadsheet_path '/Assessments/MMSE.csv'];
NEUROBAT_spreadsheet = [ADNI_spreadsheet_path '/Assessments/NEUROBAT.csv'];

viscode_cell = repmat({'m12'}, length(rid_cell), 1);
behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, ... 
    NEUROBAT_spreadsheet, 'ADNI2', repmat({'ADNI2'}, length(rid), 1), rid_cell, viscode_cell);

% create cellarray and write it into csv file
subinfo = [rid gender age dx icv gmVol behavior_scores];
subinfo_cell = CBIG_MMLDA_matrix2cellstr(subinfo);
var_names = [{'RID', 'GENDER', 'AGE', 'DX', 'ICV', 'gmVol'} ADAS.DESCRIPTION' MMSE.DESCRIPTION' NEUROBAT.DESCRIPTION'];
cell2csv([out_dir '/ADNI2_m12_subinfo.csv'], [var_names; subinfo_cell] , ',')

%%%
% ADNI1 bl
%%%
% get rid
rid_file = [ADNI_VBM_path '/preprocessing/lists/ADNI1_bl/ADNI1_bl_rid_rm_seg_fail.txt'];
rid = load(rid_file);
rid_cell = CBIG_MMLDA_matrix2cellstr(rid);
phase_cell = repmat({'ADNI2'}, length(rid_cell), 1);

% get age gender dx by combining 'sc' and 'bl'
DXSUM_file = [ADNI_spreadsheet_path '/Assessments/DXSUM_PDXCONV_ADNIALL.csv'];
PTDEMOG_file = [ADNI_spreadsheet_path '/Subject_Characteristics/PTDEMOG.csv'];

viscode_cell = repmat({'sc'}, length(rid_cell), 1);
[sc_age, sc_gender, sc_dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);
viscode_cell = repmat({'bl'}, length(rid_cell), 1);
[bl_age, bl_gender, bl_dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);
age = sc_age;
age(isnan(sc_age)) = bl_age(isnan(sc_age));
gender = sc_gender;
gender(isnan(sc_gender)) = bl_gender(isnan(sc_gender));
dx = sc_dx;
dx(isnan(sc_dx)) = bl_dx(isnan(sc_dx));

% get icv gmVol
gmVolAndICV_file = [ADNI_VBM_path '/preprocessing/output/ADNI1_bl/gmVolAndICV/id_gmVol_icv.mat']; 
load(gmVolAndICV_file)
icv = cell2mat(id_gmVol_icv(:, 3));
gmVol = cell2mat(id_gmVol_icv(:, 2));

% get behavioral scores by combing 'sc' and 'bl'
ADAS_ADNI1_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADASSCORES.csv'];
ADAS_ADNI2_spreadsheet = [ADNI_spreadsheet_path '/Assessments/ADAS_ADNIGO2.csv'];
MMSE_spreadsheet = [ADNI_spreadsheet_path '/Assessments/MMSE.csv'];
NEUROBAT_spreadsheet = [ADNI_spreadsheet_path '/Assessments/NEUROBAT.csv'];

viscode_cell = repmat({'bl'}, length(rid_cell), 1);
bl_behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, ... 
    NEUROBAT_spreadsheet, 'ADNI1', repmat({'ADNI1'}, length(rid), 1), rid_cell, viscode_cell);
viscode_cell = repmat({'sc'}, length(rid_cell), 1);
sc_behavior_scores = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, ...
    NEUROBAT_spreadsheet, 'ADNI1', repmat({'ADNI1'}, length(rid), 1), rid_cell, viscode_cell);
bl_sc_behavior_scores = bl_behavior_scores;
bl_sc_behavior_scores(isnan(bl_behavior_scores)) = sc_behavior_scores(isnan(bl_behavior_scores));

% create cellarray and write it into csv file
subinfo = [rid gender age dx icv gmVol bl_sc_behavior_scores];
subinfo_cell = CBIG_MMLDA_matrix2cellstr(subinfo);
var_names = [{'RID', 'GENDER', 'AGE', 'DX', 'ICV', 'gmVol'} ADAS.DESCRIPTION' MMSE.DESCRIPTION' NEUROBAT.DESCRIPTION'];
cell2csv([out_dir '/ADNI1_bl_subinfo.csv'], [var_names; subinfo_cell] , ',')

%%%
% ADNI2 or ADNI3 PET
%%%
% get rid viscode
rid_viscode = importdata([ADNI_VBM_path '/preprocessing/lists/ADNI23_PETtau/RID_Viscode_269all.txt']);
tmp = regexp(rid_viscode, '_', 'split');
tmp_cell = vertcat(tmp{:});
rid_cell = tmp_cell(:, 1);
viscode_cell = tmp_cell(:, 2);
rid = CBIG_MMLDA_cellstr2matrix(rid_cell);
phase_cell = importdata([ADNI_VBM_path '/preprocessing/lists/ADNI23_PETtau/Phase_269all.txt']);

% get age gender dx
DXSUM_file = [ADNI_DOC_DIR '/All/ADNI_180413/documentation/DXSUM_PDXCONV_ADNIALL.csv'];
PTDEMOG_file = [ADNI_DOC_DIR '/All/ADNI_180413/documentation/PTDEMOG.csv'];
[age, gender, dx] = CBIG_MMLDA_get_age_gender_dx(DXSUM_file, PTDEMOG_file, phase_cell, rid_cell, viscode_cell);

% get icv gmVol 
load([getenv('CBIG_MMLDA_ADNI_DIR') 'Sun2019_SPMVBM/preprocessing/output/ADNI23_PETtau/gmVolAndICV/id_gmVol_icv.mat']);
icv = cell2mat(id_gmVol_icv(:, 3));
gmVol = cell2mat(id_gmVol_icv(:, 2));

% get behavioral scores
ADAS_ADNI1_spreadsheet = [ADNI_DOC_DIR '/All/ADNI_161017/documentation/ADASSCORES.csv'];
ADAS_ADNI2_spreadsheet = [ADNI_DOC_DIR '/All/ADNI_180413/documentation/ADAS_ADNIGO23.csv'];
MMSE_spreadsheet = [ADNI_DOC_DIR '/All/ADNI_180413/documentation/MMSE.csv'];
NEUROBAT_spreadsheet = [ADNI_DOC_DIR '/All/ADNI_180413/documentation/NEUROBAT.csv'];
behavior_scores_bl = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, NEUROBAT_spreadsheet, 'ADNI3', phase_cell, rid_cell, viscode_cell);
behavior_scores_sc = CBIG_MMLDA_get_behavior_scores(ADAS_ADNI1_spreadsheet, ...
    ADAS_ADNI2_spreadsheet, MMSE_spreadsheet, NEUROBAT_spreadsheet, 'ADNI3', ...
    phase_cell, rid_cell, strrep(viscode_cell, 'bl', 'sc'));
behavior_scores_bl_sc = behavior_scores_bl;
behavior_scores_bl_sc(isnan(behavior_scores_bl)) = behavior_scores_sc(isnan(behavior_scores_bl));


% create cellarray and write it into csv file
subinfo = [rid gender age dx icv gmVol behavior_scores_bl_sc];
subinfo_cell = CBIG_MMLDA_matrix2cellstr(subinfo);
var_names = [{'RID', 'GENDER', 'AGE', 'DX', 'ICV', 'gmVol'} ADAS.DESCRIPTION' MMSE.DESCRIPTION' NEUROBAT.DESCRIPTION'];
cell2csv([out_dir '/ADNI23_PETtau_subinfo.csv'], [var_names; subinfo_cell] , ',')

rmpath([CBIG_CODE_DIR '/stable_projects/disorder_subtypes/Sun2019_ADJointFactors/utilities'])

# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# UKBB
## FNN
### Age
python3 ../cbig/He2019/CBIG_ukbb_fnn.py --pred_item 1 --weight_decay 2.798810e-05 \
--lr 9.122903e-03 --dropout 0.308501 --scheduler_decrease 32 --n_l1 9 --n_layer 2 --gpu 0
# Best validation at index:  59
# Average validation corr: 0.608733191418575 , MAE: 4.87371812687222
# Average test corr 0.5944267356378252 , MAE 4.910241783172809
# Final ensemble test corr 0.598939930395491 , MAE 4.899412085222512
# file saved: log/fnn_pred_1_2019_05_20_22_39.npz
# time spent: 1947.9170

### Pairs matching
python3 ../cbig/He2019/CBIG_ukbb_fnn.py --pred_item 2 --weight_decay 1.140774e-06 \
--lr 1.456960e-03 --dropout 0.284967 --scheduler_decrease 196 --n_l1 415 --n_l2 437 --gpu 1
# Best validation at index:  41
# Average validation corr: 0.06095474783625908 , MAE: 0.5714195852144149
# Average test corr 0.03248601991646989 , MAE 0.5828878256063926
# Final ensemble test corr 0.04550521868479115 , MAE 0.5668671243996108
# file saved: log/fnn_pred_2_2019_05_20_22_38.npz
# time spent: 1874.2467

### Fluid intelligence
python3 ../cbig/He2019/CBIG_ukbb_fnn.py --pred_item 3 --weight_decay 1.424868e-04 \
--lr 1.690104e-03 --dropout 0.526189 --scheduler_decrease 113 --n_l1 318 --n_l2 357 --gpu 2
# Best validation at index:  146
# Average validation corr: 0.2388634849080966 , MAE: 1.6198834025237037
# Average test corr 0.23012958759505997 , MAE 1.6155740883261185
# Final ensemble test corr 0.23859183696516087 , MAE 1.613521670429741
# file saved: log/fnn_pred_3_2019_05_20_22_37.npz
# time spent: 1844.6149

### Sex
python3 ../cbig/He2019/CBIG_ukbb_fnn_sex.py --weight_decay 2.662384e-04 \
--lr 1.132743e-04 --dropout 0.002753 --scheduler_decrease 150 --n_l1 3 --n_layer 2 --gpu 3
# Best validation at index:  150
# Average validation aucc: 0.8828000426292419
# Average test aucc: 0.8988000273704528
# Final averaged test aucc 0.916
# file saved: log/fnn_sex_2019_05_20_22_38.npz
# time spent: 1861.2247

## BrainNetCNN
### Age
python3 ../cbig/He2019/CBIG_ukbb_brainnetcnn.py --pred_item 1 --weight_decay 4.023178e-04 \
--lr 1.237070e-04 --dropout 0.573231 --scheduler_decrease 146 --e2e 22 --e2n 79 \
--n2g 91 --gpu 0
# Best validation at index:  121
# Average validation corr: 0.6080955422351598 , MAE: 4.785217372621038
# Average test corr 0.5912458592449898 , MAE 4.850027225787653
# Final ensemble test corr 0.5975887468775711 , MAE 4.8235957190804575
# file saved: log/brainnetcnn_pred_1_2019_05_21_00_32.npz
# time spent: 8694.7439

### Pairs matching
python3 ../cbig/He2019/CBIG_ukbb_brainnetcnn.py --pred_item 2 --weight_decay 1.286714e-05 \
--lr 8.659643e-03 --dropout 0.264064 --scheduler_decrease 57 --e2e 27 --e2n 29 \
--n2g 54 --gpu 1
# Best validation at index:  34
# Average validation corr: 0.06915346293643472 , MAE: 0.5470395188680041
# Average test corr 0.06351325822564248 , MAE 0.5568478520835881
# Final ensemble test corr 0.06715597865068382 , MAE 0.5533606540988921
# file saved: log/brainnetcnn_pred_2_2019_05_21_00_26.npz
# time spent: 8307.3032

### Fluid intelligence
python3 ../cbig/He2019/CBIG_ukbb_brainnetcnn.py --pred_item 3 --weight_decay 1.010687e-04 \
--lr 1.152218e-03 --dropout 0.776375 --scheduler_decrease 126 --e2e 40 --e2n 60 \
--n2g 41 --gpu 2
# Best validation at index:  87
# Average validation corr: 0.25079064221920744 , MAE: 1.60885481216976
# Average test corr 0.23131304675585582 , MAE 1.61393375143398
# Final ensemble test corr 0.23538333496825428 , MAE 1.6104855664108566
# file saved: log/brainnetcnn_pred_3_2019_05_21_00_15.npz
# time spent: 7686.2845

### Sex
python3 ../cbig/He2019/CBIG_ukbb_brainnetcnn_sex.py --weight_decay 3.237138e-07 \
--lr 1.389995e-04 --dropout 0.463221 --scheduler_decrease 200 --e2e 38 --e2n 58 \
--n2g 7 --gpu 3
# Best validation at index:  126
# Average validation aucc: 0.8910000443458557
# Average test aucc: 0.9110000491142273
# Final averaged test aucc 0.914
# file saved: log/brainnetcnn_sex_2019_05_21_00_25.npz
# time spent: 8186.3888

## GCNN
### Age
python3 ../cbig/He2019/CBIG_ukbb_gcnn.py --pred_item 1 --graph_setup 8868Subject_corr_option_3_param_5 \
--lr 2.160718e-04 --l2_regularizer 9.180659e-07 --dropout 0.316318 --n_l1 10 --optimizer Adam --gpu 0
# Best validation at index:  273
# Average validation corr: 0.5861164431652895 , MAE: 4.953960851136457
# Average test corr 0.5826725455814895 , MAE 4.943359771762527
# Final ensemble test corr 0.5932505151246215 , MAE 4.895236611094313
# file saved: log/gcnn_pred_1_2019_05_21_03_35.npz
# time spent: 4087.4489

### Pairs matching
python3 ../cbig/He2019/CBIG_ukbb_gcnn.py --pred_item 2 --graph_setup 8868Subject_corr_option_1_param_5 \
--lr 9.898400e-03 --l2_regularizer 4.716467e-07 --dropout 0.308161 --n_l1 3 --gpu 1
# Best validation at index:  943
# Average validation corr: 0.10567290624954688 , MAE: 0.5015095236984087
# Average test corr 0.0054346627916686385 , MAE 0.5098499840006613
# Final ensemble test corr 0.007950924654701286 , MAE 0.49656696061944455
# file saved: log/gcnn_pred_2_2019_05_21_03_35.npz
# time spent: 4059.5661

### Fluid intelligence
python3 ../cbig/He2019/CBIG_ukbb_gcnn.py --pred_item 3 --graph_setup 8868Subject_corr_option_3_param_5 \
--lr 3.767978e-04 --l2_regularizer 7.182952e-04 --dropout 0.555404 --n_l1 72 --optimizer Adam --gpu 2
# Best validation at index:  176
# Average validation corr: 0.22933381754799192 , MAE: 1.6298617692250748
# Average test corr 0.2211120455674774 , MAE 1.6336887923018935
# Final ensemble test corr 0.2323602114278351 , MAE 1.612857167789439
# file saved: log/gcnn_pred_3_2019_05_21_03_36.npz
# time spent: 4078.6902

### Sex
python3 ../cbig/He2019/CBIG_ukbb_gcnn_sex.py --graph_setup 8868Subject_corr_option_3_param_5 \
--lr 1.048051e-03 --l2_regularizer 3.343659e-04 --dropout 0.015034 --n_l1 71 --optimizer Adam --gpu 3
# Best validation at index:  378
# Average validation aucc: 0.8869999999999999
# Average test aucc: 0.9071999999999999
# Final averaged test aucc 0.916
# file saved: log/gcnn_sex_2019_05_21_03_35.npz
# time spent: 4053.2327

# HCP
## FNN
python3 ../cbig/He2019/CBIG_hcp_fnn.py --dropout 0.6 --n_l1 223 --n_l2 128 --n_l3 192 \
--lr 0.01 --lr_decay 1e-4 --l2_regularizer 0.02 --gpu 0
# Optimal index for each fold at: [207 133 227 222 221 136 111 150 195 143 124 153 150 183 172 177 149 186
#  225 183]
# Optimal result for each fold: [0.0831874  0.06920238 0.02598082 0.06790322 0.12293432 0.14864182
#  0.09959732 0.09883037 0.17443786 0.11601304 0.13400943 0.15618891
#  0.05696233 0.08602551 0.11133959 0.25771613 0.0507679  0.28369494
#  0.15156227 0.12618322]
# Final test result: 0.12105893812349562
# log saved at: HCP_fnn_2019_05_22_01_17.npz
# time spent: 49290.9864

## BrainNetCNN
python ../cbig/He2019/CBIG_hcp_brainnetcnn.py --lr_decay 2.2e-07 --lr 6.66e-03 \
--dropout 0.4 --leaky_alpha 0.3 --e2e 18 --e2n 19 --n2g 84 --gpu 1
# Optimal index for each fold at: [79 88 79 89 79 80 92 98 89 82 84 85 77 96 86 81 83 76 81 79]
# Optimal result for each fold: [0.10034946 0.08520263 0.02917019 0.07020571 0.09436217 0.08464824
#  0.12025285 0.05335622 0.19127117 0.12422879 0.13161608 0.13879307
#  0.06719132 0.09390611 0.10750996 0.18411366 0.07586971 0.22671835
#  0.13921217 0.16183644]
# Final test result: 0.11399071384942239
# log saved at: HCP_brainnetcnn_2019_05_27_23_32.npz
# time spent: 272279.6233

## GCNN
python3 ../cbig/He2019/CBIG_hcp_gcnn.py --graph_setup 953Subject_corr_option_3_param_5 \
--dropout 0.3 --n_l1 256 --lr 0.01 --lr_decay 1e-5 --l2_regularizer 8e-4 --gpu 2
# Optimal index for each fold at: [188 173 175 166 207 164 176 130 163 196 181 136 144 167 147 188 141 177
#  156 153]
# Optimal result for each fold: [ 0.01547325  0.04896941 -0.00139652  0.05447468  0.08495046  0.13294162
#   0.07033918  0.04742601  0.10807401  0.10366147  0.04739561  0.09902692
#   0.0245089   0.00255759  0.07704467  0.15410892  0.09480941  0.15498725
#   0.04905781  0.07503035]
# Final test result: 0.07217205028059999
# log saved at: HCP_gcnn_2019_06_01_23_46.npz
# time spent: 207844.4465

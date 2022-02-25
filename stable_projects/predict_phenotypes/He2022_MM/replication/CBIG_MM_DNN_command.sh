#!/bin/sh
# Written by Tong He and CBIG under MIT license: https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md

# Experiment 1 new env
# The result shown here is little bit different from result in manuscript
# The reason is that we updated our conda environment due to experiment 2
# Here is the result you will get with current enviroment
# We kept "Experiment 1 old env" section where results match manuscript for reference

# Base model training
python ../cbig/He2022/CBIG_ukbb_dnn.py --gpu 0 --seed 1 --epochs 1000 --metric cod \
--weight_decay 8.447320e-04 --lr 3.645653e-03 --dropout 0.241964 \
--scheduler_decrease 312 --n_l1 87 --n_l2 386 --n_l3 313 --n_l4 349 --n_layer 4
# Best validation at index:  108
# file saved: /home/the/deepGround/code/2002/He2022_MM/replication/output_dnn/dnn_base.npz
# Average validation corr: 0.32973410861647945 , COD: 0.15799483422497265 , MAE: 0.28677470754758266
# time spent: 5293.3907

# Meta-matching with DNN
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm.py
# rng 99 at 40.751468658447266s: cor 0.14537, cod 0.00691

# Meta-matching with DNN stacking
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py
# rng 99 at 739.8069829940796s: cor 0.15651, cod 0.04084

# We added the option for restricting the alpha of KRR (stacking)
# This option is added during our development of experiment 2
# The result of experiment 1 in manuscript does not use restricted-alpha
# Here we include both for your reference 
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py --not-restricted-alpha
# rng 99 at 3325.576665878296s: cor 0.15463, cod -0.02479

# Meta-matching with DNN finetune
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_finetune.py --gpu 0
# rng 99 at 518072.8960635662s: cor 0.14537, tuned 0.14951, cod 0.00691, tuned 0.00551


# Experiment 1 old env 
# (for reference only)
# (conda envoriment in config folder will get result above "Experiment 1 new env")

# Base model training
python ../cbig/He2022/CBIG_ukbb_dnn.py --gpu 0 --seed 1 --epochs 1000 --metric cod \
--weight_decay 8.447320e-04 --lr 3.645653e-03 --dropout 0.241964 \
--scheduler_decrease 312 --n_l1 87 --n_l2 386 --n_l3 313 --n_l4 349 --n_layer 4
# Best validation at index:  118
# file saved: /home/user/deepGround/code/1906/metaLearning/He2022_MM/replication/output_dnn/dnn_base.npz
# Average validation corr: 0.3341948884783289 , COD: 0.160241500064145 , MAE: 0.2867968463373428
# time spent: 5107.5932

# Meta-matching with DNN
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm.py
#rng 99 at 1683.1999824047089s: cor 0.14556, cod 0.00473

# Meta-matching with DNN stacking
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py
#rng 99 at 5008.264993429184s: cor 0.15607, cod -0.01366

# Meta-matching with DNN finetune
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_finetune.py --gpu 0
# rng 99 at 238891.5737681389s: cor 0.14556, tuned 0.14968, cod 0.00473, tuned 0.00472


# Experiment 2
# Base model training
python ../cbig/He2022/CBIG_ukbb_dnn.py --across_dataset True --gpu 0 --seed 1 --epochs 1000 \
--metric cod --weight_decay 1.062206e-07 --lr 6.876846e-03 --dropout 0.288360 \
--scheduler_decrease 64 --n_l1 118 --n_l2 445 --n_l3 353 --n_l4 212 --n_layer 4
# Best validation at index:  167
# file saved: /home/the/deepGround/code/2002/He2022_MM/replication/output_dnn/dnn_base.npz
# Average validation corr: 0.31945954016739897 , COD: 0.14196395888981858 , MAE: 158.01032618455525
# time spent: 48654.9432

# Meta-matching with DNN
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm.py --across-dataset
# rng 99 at 48.09036111831665s: cor 0.14691, cod -0.01775

# Meta-matching with DNN stacking
python3 ../cbig/He2022/CBIG_ukbb_dnn_mm_stacking.py --across-dataset
# rng 99 at 735.3229515552521s: cor 0.15609, cod 0.02795
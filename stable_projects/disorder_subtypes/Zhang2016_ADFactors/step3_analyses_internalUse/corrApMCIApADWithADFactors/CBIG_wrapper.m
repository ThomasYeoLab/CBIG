close all;
clear;
clc;

%% Factors estiamted with 147 amyloid+ MCI

resultDir = '../../lda/postprocess/147blApMCI/k3/';
r = dir([resultDir 'r*']);
beta_apMCI = exp(load([resultDir r.name '/final.beta']));

%% Factors estiamted with 188 AD

resultDir = '../../lda/postprocess/188blAD/k3/';
r = dir([resultDir 'r*']);
beta_AD = exp(load([resultDir r.name '/final.beta']));

%% Factors estiamted with 91 amyloid+ AD

resultDir = '../../lda/postprocess/91blApAD/k3/';
r = dir([resultDir 'r*']);
beta_apAD = exp(load([resultDir r.name '/final.beta']));

%% Hungarian matching and computing average correlations

% 91 with 188
avgCorr_91with188 = CBIG_hunMatch(beta_AD, beta_apAD)

% 147 with 188
avgCorr_147with188 = CBIG_hunMatch(beta_AD, beta_apMCI)

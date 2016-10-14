function [mu, sigma2] = CBIG_compute_mu_sigma2(stats, idx1, idx2)

try % LME
    beta1 = stats.Bhat(idx1);
    beta2 = stats.Bhat(idx2);
    beta1_var = stats.CovBhat(idx1, idx1);
    beta2_var = stats.CovBhat(idx2, idx2);
    beta12_cov = stats.CovBhat(idx1, idx2);
catch % glmfit
    beta1 = stats.beta(idx1);
    beta2 = stats.beta(idx2);
    beta1_var = stats.covb(idx1, idx1);
    beta2_var = stats.covb(idx2, idx2);
    beta12_cov = stats.covb(idx1, idx2);
end

mu = beta1-beta2;

sigma2 = beta1^2+beta1_var+beta2^2+beta2_var-2*(beta12_cov+beta1*beta2)-(beta1-beta2)^2;

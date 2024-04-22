function [sigma, beta, gamma] = parameter_personalityF(w, pb_w, b, personality, alpha, Fs_mean, Fb_mean, Fbs_mean, F_juv_mean)

% create the parameter arrays
sigma= zeros(w,b); beta = sigma; gamma = sigma;

% Survival, sigma. States order in the array: PB, SB, FB, PSB, PFB, NB
surv_PB = F_juv_mean.mean_phi_Pre'; % this goes up to 16, but I want it to go to 32
sigma(:,1) = [surv_PB, repelem(surv_PB(end),w-pb_w)];

sigma(:,2) = invlogit(logit(Fs_mean.mean_Fphi_SB(1:w)) + alpha(1)*personality);
sigma(:,3) = invlogit(logit(Fs_mean.mean_Fphi_FB(1:w)) + alpha(1)*personality);
sigma(:,4) = invlogit(logit(Fs_mean.mean_Fphi_PSB(1:w)) + alpha(1)*personality);
sigma(:,5) = invlogit(logit(Fs_mean.mean_Fphi_PFB(1:w)) + alpha(1)*personality);
sigma(:,6) = invlogit(logit(Fs_mean.mean_Fphi_NB(1:w)) + alpha(1)*personality);
sigma(isnan(sigma))=0;

% Breeding, beta. States order: PB, SB, FB, PSB, PFB, NB
breed_PB = F_juv_mean.mean_rho_Pre'; % this goes up to 16, but I want it to go to 32
beta(:,1) = [breed_PB, repelem(breed_PB(end),w-pb_w)];

beta(:,2) = invlogit(logit(Fs_mean.mean_Frho_SB(1:w)) + alpha(2)*personality);
beta(:,3) = invlogit(logit(Fs_mean.mean_Frho_FB(1:w)) + alpha(2)*personality);
beta(:,4) = invlogit(logit(Fs_mean.mean_Frho_PSB(1:w)) + alpha(2)*personality);
beta(:,5) = invlogit(logit(Fs_mean.mean_Frho_PFB(1:w)) + alpha(2)*personality);
beta(:,6) = invlogit(logit(Fs_mean.mean_Frho_NB(1:w)) + alpha(2)*personality);
beta(isnan(beta))=0;

% Survival, sigma. States order: PB, SB, FB, PSB, PFB, NB
succ_PB = F_juv_mean.mean_pi_Pre'; % this goes up to 16, but I want it to go to 32
gamma(:,1) = [succ_PB, repelem(succ_PB(end),w-pb_w)];

gamma(:,2) = invlogit(logit(Fs_mean.mean_Fpi_SB(1:w)) + alpha(3)*personality);
gamma(:,3) = invlogit(logit(Fs_mean.mean_Fpi_FB(1:w)) + alpha(3)*personality);
gamma(:,4) = invlogit(logit(Fs_mean.mean_Fpi_PSB(1:w)) + alpha(3)*personality);
gamma(:,5) = invlogit(logit(Fs_mean.mean_Fpi_PFB(1:w)) + alpha(3)*personality);
gamma(:,6) = invlogit(logit(Fs_mean.mean_Fpi_NB(1:w)) + alpha(3)*personality);
gamma(isnan(gamma))=0;

end


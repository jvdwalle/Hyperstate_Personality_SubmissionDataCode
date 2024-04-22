function [sigma, beta, gamma] = parameter_personalityM(w, pb_w, b, personality, alpha, Ms_mean, Mb_mean, Mbs_mean, M_juv_mean)

% create the parameter arrays
sigma= zeros(w,b); beta = sigma; gamma = sigma;

% Survival, sigma. States order in the array: PB, SB, FB, PSB, PFB, NB
surv_PB = M_juv_mean.mean_phi_Pre'; % this goes up to 16, but I want it to go to 32
sigma(:,1) = [surv_PB, repelem(surv_PB(end),w-pb_w)];

sigma(:,2) = invlogit(logit(Ms_mean.mean_Mphi_SB(1:w)) + alpha(1)*personality);
sigma(:,3) = invlogit(logit(Ms_mean.mean_Mphi_FB(1:w)) + alpha(1)*personality);
sigma(:,4) = invlogit(logit(Ms_mean.mean_Mphi_PSB(1:w)) + alpha(1)*personality);
sigma(:,5) = invlogit(logit(Ms_mean.mean_Mphi_PFB(1:w)) + alpha(1)*personality);
sigma(:,6) = invlogit(logit(Ms_mean.mean_Mphi_NB(1:w)) + alpha(1)*personality);
sigma(isnan(sigma))=0;

% Breeding, beta. States order: PB, SB, FB, PSB, PFB, NB
breed_PB = M_juv_mean.mean_rho_Pre'; % this goes up to 16, but I want it to go to 32
beta(:,1) = [breed_PB, repelem(breed_PB(end),w-pb_w)];

beta(:,2) = invlogit(logit(Ms_mean.mean_Mrho_SB(1:w)) + alpha(2)*personality);
beta(:,3) = invlogit(logit(Ms_mean.mean_Mrho_FB(1:w)) + alpha(2)*personality);
beta(:,4) = invlogit(logit(Ms_mean.mean_Mrho_PSB(1:w)) + alpha(2)*personality);
beta(:,5) = invlogit(logit(Ms_mean.mean_Mrho_PFB(1:w)) + alpha(2)*personality);
beta(:,6) = invlogit(logit(Ms_mean.mean_Mrho_NB(1:w)) + alpha(2)*personality);
beta(isnan(beta))=0;

% Survival, sigma. States order: PB, SB, FB, PSB, PFB, NB
succ_PB = M_juv_mean.mean_pi_Pre'; % this goes up to 16, but I want it to go to 32
gamma(:,1) = [succ_PB, repelem(succ_PB(end),w-pb_w)];

gamma(:,2) = invlogit(logit(Ms_mean.mean_Mpi_SB(1:w)) + alpha(3)*personality);
gamma(:,3) = invlogit(logit(Ms_mean.mean_Mpi_FB(1:w)) + alpha(3)*personality);
gamma(:,4) = invlogit(logit(Ms_mean.mean_Mpi_PSB(1:w)) + alpha(3)*personality);
gamma(:,5) = invlogit(logit(Ms_mean.mean_Mpi_PFB(1:w)) + alpha(3)*personality);
gamma(:,6) = invlogit(logit(Ms_mean.mean_Mpi_NB(1:w)) + alpha(3)*personality);
gamma(isnan(gamma))=0;

end


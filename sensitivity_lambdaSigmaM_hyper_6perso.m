% Sensitivity analyses  (run model first)
% Authors: Joanie Van de Walle, Silke van Daalen

%% SENSITIVITY OF LAMBDA TO SURVIVAL, SIGMA

%% Preliminary steps - useful vectors and matrices

%  Get important Identity and permutation matrices
[Is1, Iss1, Is1square, Kss1] = Imat(s1, s);
%[Is2, Iss2, Is2square, Kss2] = Imat(s2, s);
%[Is3, Iss3, Is3square, Kss3] = Imat(s3, s);


% Get useful matrices
% Get the Ujks matrices
Ujks_m = {}; % initiate empty cell array

for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
        Ujk_m = a_m{1}(:,:,j,k);
        Ujks_m{end + 1} = Ujk_m;
    end
end

% Get useful vectors
% Get the sigmajk vectors
sigmajks_m = {};
for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
        sigmajk_m = THETA_mean{4}(:,j,k); % The first cell of THETA is sigma for males
        sigmajks_m{end + 1} = sigmajk_m;
    end
end

% Get the betajk vectors
betajks_m = {};
for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
        betajk_m = THETA_mean{5}(:,j,k); % The second cell of THETA is beta for females
        betajks_m{end + 1} = betajk_m;
    end
end

% Get the gammajk vectors
gammajks_m = {};
for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
        gammajk_m = THETA_mean{6}(:,j,k); % The third cell of THETA is gamma for females
        gammajks_m{end + 1} = gammajk_m;
    end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two pathways from Lambda to sigma (do one at a time)
% 1) lambda - Atilde - Utilde - Ublock - Ujk - sigma
% 2) lambda - Atilde - Ftilde - Rblock - Rjk - sigma

%%%%%%%%%%%%%%%%%%
% Equations
% Sentitivity of lambda to sigma (dlambda_dsigma)
% 1 dlambda_dsigma = dlambda_dvecAtilde*dvecAtilde_dsigma;
% 2 dlambda_dvecAtilde = kron(w_vec', v_vec');
% 3 dvecAtilde_dsigma = dvecUtilde_dsigma + dvecFtilde_dsigma

% Because this is computationnally extensive, do one path at a time
% path 1 - Utilde
% 4 dvecUtilde_dsigma = dvecUtilde_dvecUblock*dvecUblock_dsigma
% 5 dvecUtilde_dvecUblock = kron(sparse(right_Ublock), sparse(left_Ublock))
% 6 dvecUblock_dsigma = dvecUblock_dvecUjk*dvecUjk_dsigma
% 7 dvecUblock_dvecUjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 8 dvecUjk_dsigma = DvecYud(kron(Iw' , Onew)

% path 2 - Ftilde
% 9 dvecFtilde_dsigma = dvecFtilde_dvecRblock*dvecRblock_dsigma
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock), sparse(left_Rblock))
% 11 dvecRblock_dsigma = dvecRblock_dvecRjk*dvecRjk_dsigma
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dsigma =

%%
% 1 dlambda_dsigma = dlambda_dvecAtilde*dvecAtilde_dsigma;
%   2 dlambda_dvecAtilde = kron(w_vec', v_vec'); (eq. 90 in Caswell et al. 2018)

% % Get eigenvectors
% [wmat_m, dmat, vmat]=eig(A_tilde);
% lambda=diag(dmat);
% imax=find(lambda==max(lambda));
% lambda1=lambda(imax); % this is the dominant eigenvalue, i.e. population growth rate
% %clear lambda


[w_m,lam_m,v_m]=eig(A_tilde_m);
lambda_m = diag(lam_m);
[I_m,~]=find(lam_m==max(diag(lam_m)));
lambda1_m=lambda_m(I_m);
w_m=w_m(:,I_m);
v_m=v_m(:,I_m);
dlambda_dvecAtilde_m=kron(w_m',v_m')/(v_m'*w_m);

% % stable age-stage distribution, w = right eigenvector
% w_vec=wmat(:,imax); % right eigenvector corresponding to lambda_1
% w_vec=w_vec./(sum(w_vec)); % normalize w to sum to 1
% 
% v_vec = vmat(:, imax); % left eigenvector corresponding to lambda_1
% v_vec=v_vec./(sum(v_vec)); % normalize v to sum to 1

% dlambda_dvecAtilde = kron(w_vec', v_vec')%/(v_vec'.w_vec);
%clear w_vec v_vec vmat wmat dmat

%%  3 dvecAtilde_dsigma = dvecUtilde_dsigma + dvecFtilde_dsigma


%% PATH 1 - U TILDE
%%
%  4 dvecUtilde_dsigma = dvecUtilde_dvecUblock*dvecUblock_dsigma
%  5 dvecUtilde_dvecUblock = kron(sparse(right_Ublock), sparse(left_Ublock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Ublock_m = U_m{1}; % aging process
Bblock_m = U_m{2}; % reproductive proccess (change in reproductive states)
Pblock_m = U_m{3}; % changes in personality

right_Ublock_m = Iblock';
left_Ublock_m = (K3*K2)'*Pblock_m*K3*Bblock_m*K2;

dvecUtilde_dvecUblock_m = kron(sparse(right_Ublock_m), sparse(left_Ublock_m));
%clear right_Ublock left_Ublock K2 K3 Iblock Ublock Bblock Pblock


%   6 dvecUblock_dsigma = dvecUblock_dvecUjk*dvecUjk_dsigma
%   7 dvecUblock_dvecUjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
%   8 dvecUjk_dsigma = DvecYud(kron(Iw' , Onew)

% Since we are in the first dimension (Ujk), k = 1
%%% Get the kronecker prod Iss1, Kss1, Is1
part_1_u =  kron(sparse(Kss1), sparse(Is1));
bigthing = kron(sparse(Iss1), part_1_u);   % this is constant and will be applied to all the ssk matrices.
%clear part_1_u Iss1 Kss1 Is1

% Get the Yu matrix and its associated vectors
Yu = zeros(w,w);
Yu = diag(ones(1, w-1),-1);
Yu(w,w) = 1;
vecYu = Yu(:);
DvecYu = diag(vecYu);
%clear Yu vecYu

% Get the Onew and Iw matrices
Onew = ones(w, 1);
Iw = eye(w);
Iw1w = kron(Iw', Onew);
%clear Onew Iw

% Initialize output
totalsum_u_m = zeros((w*b*g)^2,w); %
sens_60_u_m = {}; % Sensitivity of all matrix entries to sigma for each of the 60 matrices
dlambda_dsigmajk_u_m= {}; % This will give the sensitivity of lambda to sigma for ages 1-31 for each of the combination of reproductive state and personality
elambda_esigmajk_u_m= {}; % This will give the elasticities of lambda to sigma for ages 1-31 for each of the combination of reproductive state and personality
% Loop over all the s/s1 (60) matrices
for i=1:s/s1
    E =  Emat(i,i,s/s1);
    % Then apply the vec transformation
    vecE = E(:);
    %clear E
    
    % Then kronecker product with the Is1square matrix
    vecEI = kron(sparse(vecE), sparse(Is1square));
    %clear vecE
    
    % Then multiply them
    dvecUblock_dvecUjk_m = bigthing*vecEI;
    %clear vecEI
    
    dvecUjk_dsigma_m = DvecYu*Iw1w;
    dvecUblock_dsigma_m =  dvecUblock_dvecUjk_m*dvecUjk_dsigma_m;
    %clear dvecUblock_dvecUjk dvecUjk_dsigma
    
    sens_60_u_m{i} = dvecUblock_dsigma_m;
    
    totalsum_u_m = totalsum_u_m + dvecUblock_dsigma_m;
    dlambda_dsigmajk_u_m{i} = dlambda_dvecAtilde_m*dvecUtilde_dvecUblock_m*dvecUblock_dsigma_m;
    elambda_esigmajk_u_m{i} = (dlambda_dsigmajk_u_m{i}.*sigmajks_m{i}')/lambda1_m;
    %clear dvecUblock_dsigma
end


%%%% Finalise calculation for the path between sigma and U
dvecUblock_dsigma_m = totalsum_u_m;
dvecUtilde_dsigma_m = dvecUtilde_dvecUblock_m*dvecUblock_dsigma_m;


% Save the output
save('dvecUtilde_dsigma_m_6perso.mat', "dvecUtilde_dsigma_m");
save("dlambda_dsigmajk_u_m_6perso.mat", "dlambda_dsigmajk_u_m")
save("elambda_esigmajk_u_m_6perso.mat", "elambda_esigmajk_u_m")





%% PATH 2 - F TILDE
%%
% 9 dvecFtilde_dsigma = dvecFtilde_dvecRblock*dvecRblock_dsigma
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock), sparse(left_Rblock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Rblock_m = F_m{1};
Fblock_m = F_m{2};
Hblock_m = F_m{3};

right_Rblock_m = Iblock';
left_Rblock_m = (K3*K2)'*Hblock_m*K3*Fblock_m*K2;

dvecFtilde_dvecRblock_m = kron(sparse(right_Rblock_m), sparse(left_Rblock_m));


% 11 dvecRblock_dsigma = dvecRblock_dvecRjk*dvecRjk_dsigma
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dsigma =


% Since we are in the first dimension (Rjk), k = 1
%%% Get the kronecker prod Iss1, Kss1, Is1
part_1_r =  kron(sparse(Kss1), sparse(Is1));
bigthing = kron(sparse(Iss1), part_1_r);   % this is constant and will be applied to all the ssk matrices.
%clear part_1_r Iss1 Kss1 Is1

% Get the Yr matrix (will be used in the loop)
Onew = ones(w, 1);
Iw = eye(w);
Iw1w = kron(Iw', Onew);
Yr = zeros(w,w);
Yr(1,:) = 1;

% Initiate the summation
totalsum_r_m = zeros((w*b*g)^2,31);
sens_60_r_m = {};
dlambda_dsigmajk_r_m= {}; % This will give the sensitivity of lambda to sigma for ages 1-31 for each of the combination of reproductive state and personality
elambda_esigmajk_r_m = {}; % This will give the elasticity of lambda to sigma for for ages 1-31 for each of the combination of reproductive state and personality
% Loop over all the s/s1 (60) matrices
for i=1:s/s1
    E =  Emat(i,i,s/s1);
    % Then apply the vec transformation
    vecE = E(:);
    % Then kronecker product with the Is1square matrix
    vecEI = kron(sparse(vecE), sparse(Is1square));
    
    % Then multiply them
    dvecRblock_dvecRjk_m = bigthing*vecEI; % this is equation 12
    
    indiag = (rho*Yr).*(Onew*betajks_m{i}').*(Onew*gammajks_m{i}');
    dvec_indiag = diag(indiag(:));
    
    dvecRjk_dsigma_m = dvec_indiag*Iw1w; % this is equation 13
    dvecRblock_dsigma_m = dvecRblock_dvecRjk_m*dvecRjk_dsigma_m;
    
    sens_60_r_m{i} = dvecRblock_dsigma_m;
    
    totalsum_r_m = totalsum_r_m + dvecRblock_dsigma_m; % this is equation 11
    dlambda_dsigmajk_r_tmp_m = dlambda_dvecAtilde_m*dvecFtilde_dvecRblock_m*dvecRblock_dsigma_m;
    dlambda_dsigmajk_r_m{i} = dlambda_dsigmajk_r_tmp_m;
    elambda_esigmajk_r_m{i} = (dlambda_dsigmajk_r_tmp_m.*sigmajks_m{i}')/lambda1_m;
end

%%%% Finalise calculation for the path between sigma and F
dvecRblock_dsigma_m = totalsum_r_m;
dvecFtilde_dsigma_m = dvecFtilde_dvecRblock_m*dvecRblock_dsigma_m;

% Save the output
save('dvecFtilde_dsigma_m_6perso.mat', "dvecFtilde_dsigma_m");
save("dlambda_dsigmajk_r_m_6perso.mat", "dlambda_dsigmajk_r_m");
save("elambda_esigmajk_r_m_6perso.mat", "elambda_esigmajk_r_m");



%%
% Combine all the information from the two paths to get total sensitivities

% Upload the sensitivities
load("dvecUtilde_dsigma_m_6perso.mat");
load("dvecFtilde_dsigma_m_6perso.mat");

load("dlambda_dsigmajk_u_m_6perso.mat")%, "dlambda_dsigmajk_u");
load("dlambda_dsigmajk_r_m_6perso.mat")%, "dlambda_dsigmajk_r");

load("elambda_esigmajk_u_m_6perso.mat")%, "elambda_esigmajk_u");
load("elambda_esigmajk_r_m_6perso.mat")%, "elambda_esigmajk_r");

% Combine them along the two paths
dlambda_dsigma_u_m = dlambda_dvecAtilde_m*(dvecUtilde_dsigma_m);
dlambda_dsigma_m_m = dlambda_dvecAtilde_m*(dvecFtilde_dsigma_m);
dlambda_dsigma_m = dlambda_dvecAtilde_m*(dvecUtilde_dsigma_m + dvecFtilde_dsigma_m);


plot(dlambda_dsigma_m)
hold on
plot(dlambda_dsigma_u_m)
hold on
plot(dlambda_dsigma_m_m)
hold off



% Get the vector
Sens_lambda_sigma = [];
Elas_lambda_sigma = [];
for jk = 1:b*g
    sens_tmp_u = dlambda_dsigmajk_u_m{jk};
    sens_tmp_r = dlambda_dsigmajk_r_m{jk};
    sens_tmp = sens_tmp_u + sens_tmp_r;
    Sens_lambda_sigma = [Sens_lambda_sigma, sens_tmp];
    
    elas_tmp_u = elambda_esigmajk_u_m{jk};
    elas_tmp_r = elambda_esigmajk_r_m{jk};
    elas_tmp = elas_tmp_u + elas_tmp_r;
    Elas_lambda_sigma = [Elas_lambda_sigma, elas_tmp];
end

Sens_lambda_sigma = Sens_lambda_sigma';
Elas_lambda_sigma = Elas_lambda_sigma';

save('Sens_lambda_sigma_m_6perso.mat', "Sens_lambda_sigma");
save('Elas_lambda_sigma_m_6perso.mat', "Elas_lambda_sigma");

% 
% 
% % 1) Get overall sensitivities and elasticities across ages
% sens_sigma_age_m = zeros(w, 1);
% elas_sigma_age_m = zeros(w,1);
% for jk = 1:b*g
%     tmp1_u = dlambda_dsigmajk_u_m{jk};
%     tmp1_r = dlambda_dsigmajk_r_m{jk};
%     tmp1 = tmp1_u + tmp1_r;
%     tmp1 = tmp1';
%     sens_sigma_age_m = sens_sigma_age_m + tmp1;
%     
%     tmp2_u = elambda_esigmajk_u_m{jk};
%     tmp2_r = elambda_esigmajk_r_m{jk};
%     tmp2 = tmp2_u + tmp2_r;
%     tmp2 = tmp2';
%     elas_sigma_age_m = elas_sigma_age_m + tmp2;
% end
% 
% % make nice plot
% font=14;
% hfont=12;
% cco=[0.8500 0.3250 0.0980;
%     0.9290 0.6940 0.1250;
%     0.6350 0.0780 0.1840;
%     0.4660 0.6740 0.1880;
%     0.3010 0.7450 0.9330;
%     0.6350 0.0780 0.1840];
% newcolors=brewermap(6,'Spectral');
% 
% plot(1:w, sens_sigma_age_m)
% plot(1:w, elas_sigma_age_m)
% 
% figure;
% plot(1:w, sens_sigma_age_m, 'Linewidth',2);
% ylabel("Sensitivity of \lambda to \sigma",'fontsize', font)
% xlabel("Age",'fontsize', font)
% set(gca,'fontsize',font) 
% saveas(gcf, "Sensitivity_Lambda_Sigma_Age_m.png")
% 
% 
% figure;
% plot(1:w, elas_sigma_age_m, 'Linewidth',2);
% ylabel("Elasticity of \lambda to \sigma",'fontsize', font)
% xlabel("Age",'fontsize', font)
% set(gca,'fontsize',font) 
% saveas(gcf, "Elasticity_Lambda_Sigma_Age_m.png")
% 
% 
% 
% % 2) Get overall sensitivities and elasticities across stages
% sens_age_m = zeros(b,w);
% elas_age_m = zeros(b,w);
% 
% sens_state_perso_m = zeros(b,g);
% elas_state_perso_m = zeros(b,g);
% 
% PB_mat =   [1 7   13 19 25 31 37 43 49 55];
% SB_mat =   [2 8   14 20 26 32 38 44 50 56];
% FB_mat =   [3 9   15 21 27 33 39 45 51 57];
% PSB_mat = [4 10 16 22 28 34 40 46 52 58];
% PFB_mat = [5 11 17 23 29 35 41 47 53 59];
% NB_mat =   [6 12 18 24 30 36 42 48 54 60];
% 
% State_mat = [PB_mat
%     SB_mat
%     FB_mat
%     PSB_mat
%     PFB_mat
%     NB_mat];
% 
% for j=1:b
%     mats = State_mat(j,:);
%     for jk=1:g
%         tmp1_u = dlambda_dsigmajk_u_m{mats(jk)};
%         tmp1_r  = dlambda_dsigmajk_r_m{mats(jk)};
%         tmp1 = tmp1_u + tmp1_r;
%         sens_state_perso_m(j,jk) = sum(tmp1); % this is the sum of sens across ages for each state and pero combi
%                 
%         tmp2_u = elambda_esigmajk_u_m{mats(jk)};
%         tmp2_r =  elambda_esigmajk_r_m{mats(jk)};
%         tmp2 = tmp2_u + tmp2_r;
%         elas_state_perso_m(j,jk) = sum(tmp2); % this is the sum of sens across ages for each state and pero combi
%         
%         sens_age_m(j,:) = tmp1 + sens_age_m(j,:);
%         elas_age_m(j,:) = tmp2 + elas_age_m(j,:);
%     end
% end
% 
% figure
% for i=1:b
% plot(1:w, sens_age_m(i,:), 'Linewidth',2)
% set(gca,'fontsize',font)
% hold on
% end
% colororder(newcolors)
% ylabel("Sensitivity of \lambda to \sigma",'fontsize', font)
% xlabel("Age",'fontsize', font)
% h=legend('PB', 'SB', 'FB','PSB','PFB','NB');
% set(gca,'fontsize',font)
% hold off
% saveas(gcf, "Sensitivity_Lambda_Sigma_AgeState_m.png")
% 
% 
% 
% figure
% for i=1:b
% plot(1:w, elas_age_m(i,:), 'Linewidth',2)
% set(gca,'fontsize',font)
% hold on
% end
% colororder(newcolors)
% ylabel("Elasticity of \lambda to \sigma",'fontsize', font)
% xlabel("Age",'fontsize', font)
% h=legend('PB', 'SB', 'FB','PSB','PFB','NB');
% set(gca,'fontsize',font)
% hold off
% saveas(gcf, "Elasticity_Lambda_Sigma_AgeState_m.png")
% 
% 
% 
% % 2) Get overall sensitivities and elasticities across age for breeding state and perso
% % This is integrated in the loop above
% bar(sens_state_perso_m);
% colororder(newcolors)
% ylabel("Sensitivity of \lambda to \sigma",'fontsize', font)
% xlabel("Breeding state",'fontsize', font)
% h=legend('Perso1', 'Perso2', 'Perso3','Perso4','Perso5','Perso6','Perso7','Perso8','Perso9','Perso10');
% set(gca, 'XTickLabel', {'PB' 'SB' 'FB' 'PSB' 'PFB' 'NB'})
% saveas(gcf, "Sensitivity_Lambda_Sigma_StatePerso_m.png")
% 
% bar(elas_state_perso_m);
% colororder(newcolors)
% ylabel("Elasticity of \lambda to \sigma",'fontsize', font)
% xlabel("Breeding state",'fontsize', font)
% h=legend('Perso1', 'Perso2', 'Perso3','Perso4','Perso5','Perso6','Perso7','Perso8','Perso9','Perso10');
% set(gca, 'XTickLabel', {'PB' 'SB' 'FB' 'PSB' 'PFB' 'NB'})
% saveas(gcf, "Elasticity_Lambda_Sigma_StatePerso_m.png")
% 
% 
% 
% 
% age=1:w;
% sens_by_perso_m=zeros(w,g);
% elas_by_perso_m = zeros(w,g);
% for i=1:g
%     j=i*6-5;
%     sens_by_perso_tmp_u_m = dlambda_dsigmajk_u_m{j}+dlambda_dsigmajk_u_m{j+1}+dlambda_dsigmajk_u_m{j+2}+dlambda_dsigmajk_u_m{j+3}+dlambda_dsigmajk_u_m{j+4}+dlambda_dsigmajk_u_m{j+5};
%     sens_by_perso_tmp_r_m = dlambda_dsigmajk_r_m{j}+dlambda_dsigmajk_r_m{j+1}+dlambda_dsigmajk_r_m{j+2}+dlambda_dsigmajk_r_m{j+3}+dlambda_dsigmajk_r_m{j+4}+dlambda_dsigmajk_r_m{j+5};
%     sens_by_perso_m(:,i)=sens_by_perso_tmp_u_m + sens_by_perso_tmp_r_m;
%     
%     elas_by_perso_tmp_u_m = elambda_esigmajk_u_m{j}+elambda_esigmajk_u_m{j+1}+elambda_esigmajk_u_m{j+2}+elambda_esigmajk_u_m{j+3}+elambda_esigmajk_u_m{j+4}+elambda_esigmajk_u_m{j+5};
%     elas_by_perso_tmp_r_m = elambda_esigmajk_r_m{j}+elambda_esigmajk_r_m{j+1}+elambda_esigmajk_r_m{j+2}+elambda_esigmajk_r_m{j+3}+elambda_esigmajk_r_m{j+4}+elambda_esigmajk_r_m{j+5};
%     elas_by_perso_m(:,i)=elas_by_perso_tmp_u_m + elas_by_perso_tmp_r_m;
%     
% end
% 
% 
% % figure code
% c_brew=brewermap(10,'Spectral');
% colormap(c_brew)
% font=16;
% hfont=12.5;
% scrsz = get(0,'ScreenSize'); % what resolution is your laptop screen?
% scrsz_laptop=[1 1 1280 720]; % mine is this
% scrsz=scrsz_laptop; % you could (optionally) redefine your screen size
% ax=axes;
% ax.LineStyleOrder = {'-','--'};
% 
% figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(3)/2]);
% p=plot(age,sens_by_perso_m,'LineWidth',1.5);
% xlabel('Survival at age','fontsize', font)
% ylabel('Sensitivity of \lambda to','fontsize', font)
% set(gca,'fontsize',font)
% h=legend(p,'Perso 1','Perso 2','Perso 3','Perso 4','Perso 5','Perso 6','Perso 7','Perso 8','Perso 9','Perso 10')
% h.FontSize=hfont;
% saveas(gcf, "Sensitivity_Lambda_Sigma_AgePerso_m.png")
% 
% figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(3)/2]);
% p=plot(age,elas_by_perso_m,'LineWidth',1.5);
% xlabel('Survival at age','fontsize', font)
% ylabel('Elasticity of \lambda to','fontsize', font)
% set(gca,'fontsize',font)
% h=legend(p,'Perso 1','Perso 2','Perso 3','Perso 4','Perso 5','Perso 6','Perso 7','Perso 8','Perso 9','Perso 10')
% h.FontSize=hfont;
% saveas(gcf, "Elasticity_Lambda_Sigma_AgePerso_m.png")
% 
% % mixing and marginalizing
% % e.g., sum across both stage and perso
% Iw=eye(w);
% one_bg=ones(b*g,1);
% lambda_beta_w_sumbg=kron(one_bg',Iw)*Sens_lambda_beta;
% 
% font=14;
% hfont=12;
% cco=[0.8500 0.3250 0.0980;
%     0.9290 0.6940 0.1250;
%     0.6350 0.0780 0.1840;
%     0.4660 0.6740 0.1880;
%     0.3010 0.7450 0.9330;
%     0.6350 0.0780 0.1840];
% newcolors=brewermap(6,'Spectral');
% 
% figure;
% plot(lambda_beta_w_sumbg, 'Linewidth',2);
% ylabel("Sensitivity of \lambda to \beta",'fontsize', font)
% xlabel("Age",'fontsize', font)
% set(gca,'fontsize',font)
% 


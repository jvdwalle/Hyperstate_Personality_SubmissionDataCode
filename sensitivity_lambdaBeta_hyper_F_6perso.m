% Sensitivity analyses  (run model first)
% Authors: Joanie Van de Walle, Silke van Daalen

%% SENSITIVITY OF LAMBDA TO BREEDING PROBABILITY, BETA

% Get important matrices
[Is1, Iss1, Is1square, Kss1] = Imat(s1, s);
[Is2, Iss2, Is2square, Kss2] = Imat(s2, s);


% Get useful matrices
 % Get the beta ik vectors
 betaiks_f = {};
   for k = 1:g % for each personality class
    for i = 1:w % for each breeding state
     betaik_f = THETA_mean{2}(i,:,g); % The second cell of THETA is beta for females
     betaiks_f{end + 1} = betaik_f;
    end  
   end
   
 % Get the gamma ik vectors
  gammaiks_f = {};
   for k = 1:g % for each personality class
    for i = 1:w % for each breeding state
     gammaik_f = THETA_mean{3}(i,:,g); % The third cell of THETA is gamma for females
     gammaiks_f{end + 1} = gammaik_f;
    end  
   end
   
 % Get the gammajk vectors
 gammajks_f = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     gammajk_f = THETA_mean{3}(:,j,k); % The third cell of THETA is gamma for females
     gammajks_f{end + 1} = gammajk_f;
    end  
   end
   
   
 % Get the betajk vectors
 betajks_f = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     betajk_f = THETA_mean{2}(:,j,k); % The second cell of THETA is beta for females
     betajks_f{end + 1} = betajk_f;
    end  
   end
   
 % Get the sigmajk vectors
 sigmajks_f = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     sigmajk_f = THETA_mean{1}(:,j,k); % The first cell of THETA is sigma for females
     sigmajks_f{end + 1} = sigmajk_f;
    end  
   end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two pathways from Lambda to beta (do one at a time)
% 1) lambda - Atilde - Utilde - Bblock - Bik - beta
% 2) lambda - Atilde - Ftilde - Rblock - Rjk - beta

%%%%%%%%%%%%%%%%%%
% Equations
% Sentitivity of lambda to beta (dlambda_dbeta)
% 1 dlambda_dbeta = dlambda_dvecAtilde*dvecAtilde_dbeta;
% 2 dlambda_dvecAtilde = kron(w_vec', v_vec')/v'w';
% 3 dvecAtilde_dbeta = dvecUtilde_dbeta + dvecFtilde_dbeta

% Because this is computationnally extensive, do one path at a time
% path 1 - Utilde
% 4 dvecUtilde_dbeta = dvecUtilde_dvecBblock*dvecBblock_dbeta
% 5 dvecUtilde_dvecBblock = kron(sparse(right_Bblock)', sparse(left_Bblock))
% 6 dvecBblock_dbeta = dvecBblock_dvecBik*dvecBik_dbeta
% 7 dvecBblock_dvecBik = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 8 dvecBik_dbeta = Dvec(Tplus + Tminus)(kron(Iw' , Onew) - Dvec(Yb)(kron(Iw', Onew)

% path 2 - Ftilde
% 9 dvecFtilde_dbeta = dvecFtilde_dvecRblock*dvecRblock_dbeta
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock)', sparse(left_Rblock))
% 11 dvecRblock_dbeta = dvecRblock_dvecRjk*dvecRjk_dbeta
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dbeta =




%%
% 1 dlambda_dbeta = dlambda_dvecAtilde*dvecAtilde_dbeta;
%   2 dlambda_dvecAtilde = kron(w_vec', v_vec')/v'w'; 
% 
% [w_f,lam_f,v_f]=eig(A_tilde_f);
% lambda_f = diag(lam_f);
% [I_f,~]=find(lam_f==max(diag(lam_f)));
% lambda1_f=lambda_f(I_f);
% w_f=w_f(:,I_f);
% v_f=v_f(:,I_f);

 % stable age-stage distribution, w = right eigenvector
 [wmat_mean, dmat_mean, vmat_mean]=eig(A_tilde_f);
    lambda_mean=diag(dmat_mean);
    imax_mean=find(lambda_mean==max(lambda_mean));
    lambda1_f=lambda_mean(imax_mean); % this is the dominant eigenvalue, i.e. population growth rate

 w_vec_mean=wmat_mean(:,imax_mean); % right eigenvector corresponding to lambda_1
w_f=w_vec_mean./(sum(w_vec_mean)); % normalize w to sum to 1

v_vec_mean = vmat_mean(:, imax_mean); % left eigenvector corresponding to lambda_1
v_f=v_vec_mean./(sum(v_vec_mean)); % normalize v to sum to 1
    


dlambda_dvecAtilde_f=kron(w_f',v_f')/(v_f'*w_f);

%% PATH 1 - U TILDE
%%
%  4 dvecUtilde_dbeta = dvecUtilde_dvecBblock*dvecBblock_dbeta
%  5 dvecUtilde_dvecBblock = kron(sparse(right_Bblock)', sparse(left_Bblock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Ublock_f = U_f{1}; % aging process
Bblock_f = U_f{2}; % reproductive proccess (change in reproductive states)
Pblock_f = U_f{3}; % changes in personality

right_Bblock_f = (K2*Ublock_f)';
left_Bblock_f = (K3*K2)'*Pblock_f*K3; 
 
dvecUtilde_dvecBblock_f = kron(sparse(right_Bblock_f), sparse(left_Bblock_f));

%   6 dvecUblock_dsigma = dvecUblock_dvecUjk*dvecUjk_dsigma
%   7 dvecUblock_dvecUjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
%   8 dvecUjk_dsigma = DvecYud(kron(Iw' , Onew)

% Since we are in the second dimension (Bik), k = 2
%%% Get the kronecker prod Iss2, Kss2, Is2
part_1_b =  kron(sparse(Kss2), sparse(Is2));
bigthing = kron(sparse(Iss2), part_1_b);   % this is constant and will be applied to all the ssk matrices.

%%% Decompose Bik in its constituents for those that can be calculated prior to looping
% Get the Y matrices and their associated vectors
Yb = zeros(b,b); % Negative entries of beta
Yb(1,1) = 1; Yb(4,2) = 1; Yb(5,3) = 1; Yb(6,4:6) = 1;
Eb = ones(b,b); % matrix of ones of dimension b by b
Y_plus = zeros(b,b); % Positive entries of gamma
Y_plus(2,:) = 1;
Y_minus = zeros(b,b); % Negative entries of gamma
Y_minus(3,:) = 1;

% Get the Oneb and Ib matrices
Oneb = ones(b, 1);
Ib = eye(b);
Ib1b = kron(Ib', Oneb);



% Initialize output
totalsum_b_f = zeros((w*b*g)^2,b);
sens_310_b_f = {};
dlambda_dbetaik_b_f= {}; % This will give the sensitivity of lambda to beta for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
elambda_ebetaik_b_f = {}; % This will give the elasticity of lambda to beta for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
% Loop over all the s/s2 (310) matrices
for i=1:s/s2
     E =  Emat(i,i,s/s2);
     % Then apply the vec transformation 
     vecE = E(:);
     
     % Then kronecker product with the Is1square matrix
     vecEI = kron(sparse(vecE), sparse(Is2square));
     
     % Then multiply them
     %bigthing = dvecBblock_dvecBik*vecEI;
     dvecBblock_dvecBik_f = bigthing*vecEI;
         
     %decompose Bik
     betaik_f = betaiks_f{i};
     gammaik_f = gammaiks_f{i};
     
     GAMMA_plus = Y_plus.*(Oneb*gammaik_f);
     GAMMA_minus = Y_minus.*(Eb-(Oneb*gammaik_f));
    
     GAMMAS = GAMMA_plus + GAMMA_minus;
     dvecGAMMAS = diag(GAMMAS(:));
     parta = dvecGAMMAS*Ib1b;
     
     dvecY = diag(Yb(:));
     partb = dvecY*Ib1b;
     
     dvecBik_dbeta_f = parta - partb;
     dvecBblock_dbeta_f =  dvecBblock_dvecBik_f*dvecBik_dbeta_f;
     
     sens_310_b_f{i} = dvecBblock_dbeta_f;

     totalsum_b_f = totalsum_b_f + dvecBblock_dbeta_f;
    
     dlambda_dbetaik_b_f{i} =   dlambda_dvecAtilde_f*dvecUtilde_dvecBblock_f*dvecBblock_dbeta_f;
    elambda_ebetaik_b_f{i} = (dlambda_dbetaik_b_f{i}.*betaiks_f{i})/lambda1_f;
     
    end

    
%%%% Finalise calculation for the path between sigma and U
dvecBblock_dbeta_f = totalsum_b_f;
dvecUtilde_dbeta_f = dvecUtilde_dvecBblock_f*dvecBblock_dbeta_f;


% Save the output
save('dvecUtilde_dbeta_f_6perso.mat', "dvecUtilde_dbeta_f");
save("dlambda_dbetaik_b_f_6perso.mat", "dlambda_dbetaik_b_f")
save("elambda_ebetaik_b_f_6perso.mat", "elambda_ebetaik_b_f")



%% PATH 2 - F TILDE
%%
% 9 dvecFtilde_dbeta = dvecFtilde_dvecRblock*dvecRblock_dbeta
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock'), sparse(left_Rblock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Rblock_f = F_f{1};
Fblock_f = F_f{2};
Hblock_f = F_f{3};

right_Rblock_f = Iblock';
left_Rblock_f = (K3*K2)'*Hblock_f*K3*Fblock_f*K2; 

dvecFtilde_dvecRblock_f = kron(sparse(right_Rblock_f), sparse(left_Rblock_f));

% 11 dvecRblock_dbeta = dvecRblock_dvecRjk*dvecRjk_dbeta
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dbeta =

% Since we are in the first dimension (Rjk), k = 1
%%% Get the kronecker prod Iss1, Kss1, Is1
part_1_r =  kron(sparse(Kss1), sparse(Is1));
bigthing = kron(sparse(Iss1), part_1_r);   % this is constant and will be applied to all the ssk matrices.

% Get the Yr matrix (will be used in the loop)
Onew = ones(w, 1);
Iw = eye(w);
Iw1w = kron(Iw', Onew);
Yr = zeros(w,w);
Yr(1,:) = 1;


% Initiate the summation
totalsum_r_f = zeros((w*b*g)^2,31);
sens_60_r_f = {};
dlambda_dbetajk_r_f= {}; % This will give the sensitivity of lambda to beta for ages 1-31 for each of the combination of reproductive state and personality
elambda_ebetajk_r_f = {}; % This will give the elasticity of lambda to beta for for ages 1-31 for each of the combination of reproductive state and personality

% Loop over all the s/s1 (60) matrices
for i=1:s/s1
   E =  Emat(i,i,s/s1);
    % Then apply the vec transformation
    vecE = E(:);
    % Then kronecker product with the Is1square matrix
    vecEI = kron(sparse(vecE), sparse(Is1square));
    
    % Then multiply them
    dvecRblock_dvecRjk_f = bigthing*vecEI; % this is equation 12
    
% Then multiply them
%bigthing = dvecRblock_dvecRjk*vecEI;

indiag = (rho*Yr).*(Onew*gammajks_f{i}').*(Onew*sigmajks_f{i}');
dvec_indiag = diag(indiag(:));


 dvecRjk_dbeta_f = dvec_indiag*Iw1w; % this is equation 13
    dvecRblock_dbeta_f = dvecRblock_dvecRjk_f*dvecRjk_dbeta_f;


sens_60_r_f{i} = dvecRblock_dbeta_f;

   totalsum_r_f = totalsum_r_f + dvecRblock_dbeta_f; % this is equation 11
    dlambda_dbetajk_r_tmp_f = dlambda_dvecAtilde_f*dvecFtilde_dvecRblock_f*dvecRblock_dbeta_f;
    dlambda_dbetajk_r_f{i} = dlambda_dbetajk_r_tmp_f;
    elambda_ebetajk_r_f{i} = (dlambda_dbetajk_r_tmp_f.*betajks_f{i}')/lambda1_f;
   end


%%%% Finalise calculation for the path between sigma and F
dvecRblock_dbeta_f = totalsum_r_f;
dvecFtilde_dbeta_f = dvecFtilde_dvecRblock_f*dvecRblock_dbeta_f;

% Save the output
save('dvecFtilde_dbeta_f_6perso.mat', "dvecFtilde_dbeta_f");
save("dlambda_dbetajk_r_f_6perso.mat", "dlambda_dbetajk_r_f")
save("elambda_ebetajk_r_f_6perso.mat", "elambda_ebetajk_r_f")


%%
% Combine all the information from the two paths to get total sensitivities

% Upload the sensitivities
load("dvecUtilde_dbeta_f_6perso.mat");
load("dvecFtilde_dbeta_f_6perso.mat");

load("dlambda_dbetaik_b_f_6perso.mat")%, "dlambda_dsigmajk_u_f");
load("dlambda_dbetajk_r_f_6perso.mat")%, "dlambda_dsigmajk_r_f");

load("elambda_ebetaik_b_f_6perso.mat")%, "elambda_esigmajk_u_f");
load("elambda_ebetajk_r_f_6perso.mat")%, "elambda_esigmajk_r_f");



% Combine them along the two paths ( This is where it doesn't work - not the same dimensions)
% Need to change arrangement.
% dvecUtilde_dbeta.mat is arranged stage within age within personality (ik)
% dvecFtilde_dbeta.mat is arranged age within stage within personality (jk)
% We want to rearrange dvecUtilde_dbeta.mat so it becomes age within stage within personality (jk)
% Transform into a table

% transform dlambda_dbetajk_r in a vector
dlambda_dbetajk_f = dlambda_dbetajk_r_f{1}';
elambda_ebetajk_f = elambda_ebetajk_r_f{1}';
for mat = 2:s/s1
    tmps= dlambda_dbetajk_r_f{mat}';
    tmpe= elambda_ebetajk_r_f{mat}';
    dlambda_dbetajk_f(end+1: end+w) = tmps;
    elambda_ebetajk_f(end+1: end+w) = tmpe;
end

% transform dlambda_dbetaik_b in a vector
dlambda_dbetaik_f = dlambda_dbetaik_b_f{1}';
elambda_ebetaik_f = elambda_ebetaik_b_f{1}';
for mat = 2:s/s2
    tmps = dlambda_dbetaik_b_f{mat}';
    tmpe = elambda_ebetaik_b_f{mat}';
    dlambda_dbetaik_f(end+1: end+b) = tmps;
    elambda_ebetaik_f(end+1: end+b) = tmpe;
end

% Ok, so now put this in a table
% Arrangement jk (age within stage within personality)
age_classes = [1:w];
stage_classes = [1:b];

Personality_jk = repelem(1:g, b*w)';
Stage_jk = repmat(repelem(stage_classes, w),1,g)';
Age_jk = repmat(age_classes, 1,b*g)';
Sens_jk = dlambda_dbetajk_f;
Elas_jk = elambda_ebetajk_f;

Tjk = table(Personality_jk, Stage_jk, Age_jk, Sens_jk, Elas_jk);

% Arrangement ik (stage within age within personality)
Personality_ik = repelem(1:g, b*w)';
Stage_ik = repmat(stage_classes, 1, w*g)';
Age_ik = repmat(repelem(age_classes, b)',g,1);
Sens_ik = dlambda_dbetaik_f;
Elas_ik = elambda_ebetaik_f;

Tik = table(Personality_ik, Stage_ik, Age_ik, Sens_ik, Elas_ik);


% Sort Tik to follow arrangement of Tjk
Tik_jk = sortrows(Tik, {'Personality_ik','Stage_ik','Age_ik'});

% Sum across dimensions
Sens_lambda_beta_Utilde = table2array(Tik_jk(:,4));
Sens_lambda_beta_Ftilde = table2array(Tjk(:,4));
Sens_lambda_beta = Sens_lambda_beta_Utilde + Sens_lambda_beta_Ftilde;
save('Sens_lambda_beta_f_6perso.mat', "Sens_lambda_beta");

Elas_lambda_beta_Utilde = table2array(Tik_jk(:,5));
Elas_lambda_beta_Ftilde = table2array(Tjk(:,5));
Elas_lambda_beta = Elas_lambda_beta_Utilde + Elas_lambda_beta_Ftilde;
save('Elas_lambda_beta_f_6perso.mat', "Elas_lambda_beta");


% plot(Sens_lambda_beta);
% hold on
% plot(Elas_lambda_beta);
% hold off
% 
% 
% % mixing and marginalizing
% % e.g., sum across both stage and perso
% Iw=eye(w);
% one_bg=ones(b*g,1);
% senslambda_beta_w_sumbg=kron(one_bg',Iw)*Sens_lambda_beta;
% elaslambda_beta_w_sumbg=kron(one_bg',Iw)*Elas_lambda_beta;
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
% plot(senslambda_beta_w_sumbg, 'Linewidth',2);
% ylabel("Sensitivity of \lambda to \beta",'fontsize', font)
% xlabel("Age",'fontsize', font)
% set(gca,'fontsize',font) 
% 
% figure;
% plot(elaslambda_beta_w_sumbg, 'Linewidth',2);
% ylabel("Elasticity of \lambda to \beta",'fontsize', font)
% xlabel("Age",'fontsize', font)
% set(gca,'fontsize',font) 
% 
% 
% 
% % e.g., sum across age, keep stage and perso intact
% one_w=ones(w,1);
% I_bg=eye(b*g);
% senslambda_beta_sumw_bg=kron(I_bg',one_w')*Sens_lambda_beta;
% senssumw_bg=reshape(senslambda_beta_sumw_bg,[b,g]);
% 
% elaslambda_beta_sumw_bg=kron(I_bg',one_w')*Elas_lambda_beta;
% elassumw_bg=reshape(elaslambda_beta_sumw_bg,[b,g]);
% 
% figure;
% bar(senssumw_bg')
% colororder(newcolors)
% ylabel("Sensitivity of \lambda to \beta",'fontsize', font)
% xlabel("Personality",'fontsize', font)
% h=legend('PB', 'SB', 'FB','PS','PF','NB');
% set(gca,'fontsize',font)
% h.FontSize=hfont;
% 
% 
% 
% bs=1:10;
% x_val=[bs' bs'+10 bs'+20 bs'+30 bs'+40 bs'+50]';
% figure;
% s=scatter(x_val(1,:),senssumw_bg(1,:),'filled');
% ylabel("Sensitivity of \lambda to \beta",'fontsize', font)
% xlabel("Personality within breeding state",'fontsize', font)
% s(1).MarkerFaceColor=newcolors(1,:);
% hold on
% for i=2:g
%     s=scatter(x_val(i,:),senssumw_bg(i,:),'filled');
% 
% end
% s(i).MarkerFaceColor=newcolors(i,:);
% xticks([5.5,15.5,25.5,35.5,45.5,55.5])
% xticklabels({'PB', 'SB', 'FB','PS','PF','NB'});
% set(gca,'fontsize',font)
% 
% bs=1:10;
% x_val=[bs' bs'+10 bs'+20 bs'+30 bs'+40 bs'+50]';
% figure;
% s=scatter(x_val(1,:),elassumw_bg(1,:),'filled');
% ylabel("Elasticity of \lambda to \beta",'fontsize', font)
% xlabel("Personality within breeding state",'fontsize', font)
% s(1).MarkerFaceColor=newcolors(1,:);
% hold on
% for i=2:g
%     s=scatter(x_val(i,:),elassumw_bg(i,:),'filled');
% 
% end
% s(i).MarkerFaceColor=newcolors(i,:);
% xticks([5.5,15.5,25.5,35.5,45.5,55.5])
% xticklabels({'PB', 'SB', 'FB','PS','PF','NB'});
% set(gca,'fontsize',font)
% 
% 

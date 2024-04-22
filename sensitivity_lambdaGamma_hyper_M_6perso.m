% Sensitivity analyses  (run model first)
% Authors: Joanie Van de Walle, Silke van Daalen

%% SENSITIVITY OF LAMBDA TO BREEDING SUCCESS PROBABILITY, GAMMA
%%

% Get important matrices
[Is1, Iss1, Is1square, Kss1] = Imat(s1, s);
[Is2, Iss2, Is2square, Kss2] = Imat(s2, s);

% Get useful matrices
 % Get the beta ik vectors
 betaiks_m = {};
   for k = 1:g % for each personality class
    for i = 1:w % for each breeding state
     betaik_m = THETA_mean{5}(i,:,g); % The second cell of THETA is beta for males
     betaiks_m{end + 1} = betaik_m;
    end  
   end
   
 % Get the gamma ik vectors
  gammaiks_m = {};
   for k = 1:g % for each personality class
    for i = 1:w % for each breeding state
     gammaik_m = THETA_mean{6}(i,:,g); % The third cell of THETA is gamma for males
     gammaiks_m{end + 1} = gammaik_m;
    end  
   end
   
 % Get the gammajk vectors
 gammajks_m = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     gammajk_m = THETA_mean{6}(:,j,k); % The third cell of THETA is gamma for males
     gammajks_m{end + 1} = gammajk_m;
    end  
   end
   
   
 % Get the betajk vectors
 betajks_m = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     betajk_m = THETA_mean{5}(:,j,k); % The second cell of THETA is beta for males
     betajks_m{end + 1} = betajk_m;
    end  
   end
   
 % Get the sigmajk vectors
 sigmajks_m = {};
   for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     sigmajk_m = THETA_mean{4}(:,j,k); % The first cell of THETA is sigma for males
     sigmajks_m{end + 1} = sigmajk_m;
    end  
   end
   
   
   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two pathways from Lambda to gamma (do one at a time)
% 1) lambda - Atilde - Utilde - Bblock - Bik - gamma
% 2) lambda - Atilde - Ftilde - Rblock - Rjk - gamma

%%%%%%%%%%%%%%%%%%
% Equations
% Sentitivity of lambda to gamma (dlambda_dgamma)
% 1 dlambda_dgamma = dlambda_dvecAtilde*dvecAtilde_dgamma;
% 2 dlambda_dvecAtilde = kron(w_vec', v_vec')/v'w';
% 3 dvecAtilde_dgamma = dvecUtilde_dgamma + dvecFtilde_dgamma

% Because this is computationnally extensive, do one path at a time
% path 1 - Utilde
% 4 dvecUtilde_dgamma = dvecUtilde_dvecBblock*dvecBblock_dgamma
% 5 dvecUtilde_dvecBblock = kron(sparse(right_Bblock)', sparse(left_Bblock))
% 6 dvecBblock_dgamma = dvecBblock_dvecBik*dvecBik_dgamma
% 7 dvecBblock_dvecBik = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 8 dvecBik_dgamma = Dvec(Aplus + Aminus)(kron(Iw' , Onew) - Dvec(Yg)(kron(Iw', Onew)

% path 2 - Ftilde
% 9 dvecFtilde_dgamma = dvecFtilde_dvecRblock*dvecRblock_dgamma
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock)', sparse(left_Rblock))
% 11 dvecRblock_dgamma = dvecRblock_dvecRjk*dvecRjk_dgamma
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dgamma = Dvec((rho*Yr).*(Onew*betajk').*(Onew*sigmajk'))*kron(Iw,Onew)



%%
% 1 dlambda_dgamma = dlambda_dvecAtilde*dvecAtilde_dgamma;
%   2 dlambda_dvecAtilde = kron(w_vec', v_vec')/v'w'; 

[w_m,lam_m,v_m]=eig(A_tilde_m);
lambda_m = diag(lam_m);
[I_m,~]=find(lam_m==max(diag(lam_m)));
lambda1_m=lambda_m(I_m);
w_m=w_m(:,I_m);
v_m=v_m(:,I_m);
dlambda_dvecAtilde_m=kron(w_m',v_m')/(v_m'*w_m);
   
%% PATH 1 - U TILDE
%%
%  4 dvecUtilde_dgamma = dvecUtilde_dvecBblock*dvecBblock_dgamma
%  5 dvecUtilde_dvecBblock = kron(sparse(right_Bblock)', sparse(left_Bblock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Ublock_m = U_m{1}; % aging process
Bblock_m = U_m{2}; % reproductive proccess (change in reproductive states)
Pblock_m = U_m{3}; % changes in personality

right_Bblock_m = (K2*Ublock_m)';
left_Bblock_m = (K3*K2)'*Pblock_m*K3; 
 
dvecUtilde_dvecBblock_m = kron(sparse(right_Bblock_m), sparse(left_Bblock_m));

%   6 dvecUblock_dgamma = dvecUblock_dvecUjk*dvecUjk_dgamma
%   7 dvecUblock_dvecUjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
%   8 dvecUjk_dgamma = DvecYud(kron(Iw' , Onew)


% Since we are in the second dimension (Bik), k = 2
%%% Get the kronecker prod Iss2, Kss2, Is2
part_1_b =  kron(sparse(Kss2), sparse(Is2));
bigthing = kron(sparse(Iss2), part_1_b);   % this is constant and will be applied to all the ssk matrices.


 %%% Decompose Bik in its constituents for those that can be calculated prior to looping
Yg = zeros(b,b); % Negative entries of gamma
Yg(3,:) = 1;
Eb = ones(b,b); % matrix of ones of dimension b by b
B_plus = zeros(b,b); % Positive entries of beta
B_plus(2:3,:) = 1;
B_minus = zeros(b,b); % Negative entries of beta
B_minus(1,1) = 1; B_minus(4,2) = 1; B_minus(5,3) = 1; B_minus(6,4:6) = 1; 

% Get the Oneb and Ib matrices
Oneb = ones(b, 1);
Ib = eye(b);
Ib1b = kron(Ib', Oneb);

% Initialize output
totalsum_b_m = zeros((w*b*g)^2,b);
sens_310_b_m = {};
dlambda_dgammaik_b_m= {}; % This will give the sensitivity of lambda to gamma for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
elambda_egammaik_b_m = {}; % This will give the elasticity of lambda to  gamma for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
% Loop over all the s/s2 (310) matrices
for i=1:s/s2
     E =  Emat(i,i,s/s2);
     % Then apply the vec transformation 
     vecE = E(:);
     
     % Then kronecker product with the Is1square matrix
     vecEI = kron(sparse(vecE), sparse(Is2square));
     
     
     % Then multiply them
     %bigthing = dvecBblock_dvecBik*vecEI;
     dvecBblock_dvecBik_m = bigthing*vecEI;
     
      %decompose Bik
     betaik = betaiks_m{i};
     gammaik = gammaiks_m{i};
     
     BETA_plus = B_plus.*(Oneb*betaik);
     BETA_minus = B_minus.*(Eb-(Oneb*betaik));
    
     BETAS = BETA_plus + BETA_minus;
     dvecBETAS = diag(BETAS(:));
     parta = dvecBETAS*Ib1b;
     
     dvecYg = diag(Yg(:));
     partb = dvecYg*Ib1b;
     
     %dvecB_dvecgammaik = parta - partb;
     dvecBik_dgamma_m = parta - partb;
     %insummation = bigthing*dvecB_dvecgammaik;
     dvecBblock_dgamma_m = dvecBblock_dvecBik_m*dvecBik_dgamma_m;
     
     sens_310_b_m{i} = dvecBblock_dgamma_m;
     
     totalsum_b_m = totalsum_b_m + dvecBblock_dgamma_m;
     
     dlambda_dgammaik_b_m{i} = dlambda_dvecAtilde_m*dvecUtilde_dvecBblock_m*dvecBblock_dgamma_m;
     elambda_egammaik_b_m{i} = (dlambda_dgammaik_b_m{i}.*gammaiks_m{i})/lambda1_m;
    
end

    
%%%% Finalise calculation for the path between sigma and U
dvecBblock_dgamma_m = totalsum_b_m;
dvecUtilde_dgamma_m = dvecUtilde_dvecBblock_m*dvecBblock_dgamma_m;


% Save the output
save('dvecUtilde_dgamma_m_6perso.mat', "dvecUtilde_dgamma_m");
save("dlambda_dgammaik_b_m_6perso.mat", "dlambda_dgammaik_b_m")
save("elambda_egammaik_b_m_6perso.mat", "elambda_egammaik_b_m")



%% PATH 2 - F TILDE
%%
% 9 dvecFtilde_dgamma = dvecFtilde_dvecRblock*dvecRblock_dgamma
% 10 dvecFtilde_dvecRblock = kron(sparse(right_Rblock'), sparse(left_Rblock))

K2 = vecperm_hyp(2,siz);
K3 = vecperm_hyp(3,siz);
Iblock = eye(w*b*g);

Rblock_m = F_m{1};
Fblock_m = F_m{2};
Hblock_m = F_m{3};

right_Rblock_m = Iblock';
left_Rblock_m = (K3*K2)'*Hblock_m*K3*Fblock_m*K2; 

dvecFtilde_dvecRblock_m = kron(sparse(right_Rblock_m), sparse(left_Rblock_m));


% 11 dvecRblock_dgamma = dvecRblock_dvecRjk*dvecRjk_dgamma
% 12 dvecRblock_dvecRjk = sum(kron(Issk, Kssk, Isk)*(kron(vecE, Isk2))
% 13 dvecRjk_dgamma = Dvec((rho*Yr).*(Onew*betajk').*(Onew*sigmajk'))*kron(Iw,Onew)


% Since we are in the first dimension (Rjk), k = 1
%%% Get the kronecker prod Iss1, Kss1, Is1
part_1_r =  kron(sparse(Kss1), sparse(Is1));
bigthing = kron(sparse(Iss1), part_1_r);   % this is constant and will be applied to all the ssk matrices.

 % Get the Onew and Iw matrices
 Onew = ones(w, 1);
 Iw = eye(w);
 Iw1w = kron(Iw', Onew);
 
 % Get the Yr matrix and its associated vectors
 Yr = zeros(w,w);
 Yr(1,:) = 1;
 
 % Get the Onew and Iw matrices
 Onew = ones(w, 1);
 Iw = eye(w);
 Iw1w = kron(Iw', Onew);
 
% Initiate the summation
totalsum_r_m = zeros((w*b*g)^2,31);
sens_60_r_m = {};
dlambda_dgammajk_r_m= {}; % This will give the sensitivity of lambda to gamma for ages 1-31 for each of the combination of reproductive state and personality
elambda_egammajk_r_m = {}; % This will give the elasticity of lambda to gamma for for ages 1-31 for each of the combination of reproductive state and personality

% Loop over all the s/s1 (60) matrices
for i=1:s/s1
 E =  Emat(i,i,s/s1);
 % Then apply the vec transformation
 vecE = E(:);
 
 % Then kronecker product with the Is1square matrix
 vecEI = kron(sparse(vecE), sparse(Is1square));
 
% Then multiply them
%bigthing = dvecRblock_dvecRjk*vecEI;
dvecRblock_dvecRjk_m = bigthing*vecEI; % this is equation 12

indiag = (rho*Yr).*(Onew*betajks_m{i}').*(Onew*sigmajks_m{i}');
dvec_indiag = diag(indiag(:));

 dvecRjk_dgamma_m = dvec_indiag*Iw1w; % this is equation 13
 dvecRblock_dgamma_m = dvecRblock_dvecRjk_m*dvecRjk_dgamma_m;

sens_60_r_m{i} = dvecRblock_dgamma_m;

totalsum_r_m = totalsum_r_m + dvecRblock_dgamma_m;

    dlambda_dgammajk_r_tmp_m = dlambda_dvecAtilde_m*dvecFtilde_dvecRblock_m*dvecRblock_dgamma_m;
    dlambda_dgammajk_r_m{i} = dlambda_dgammajk_r_tmp_m;
    elambda_egammajk_r_m{i} = (dlambda_dgammajk_r_tmp_m.*gammajks_m{i}')/lambda1_m;
end


%%%% Finalise calculation for the path between sigma and F
dvecRblock_dgamma_m = totalsum_r_m;
dvecFtilde_dgamma_m = dvecFtilde_dvecRblock_m*dvecRblock_dgamma_m;

% Save the output
save('dvecFtilde_dgamma_m_6perso.mat', "dvecFtilde_dgamma_m");
save("dlambda_dgammajk_r_m_6perso.mat", "dlambda_dgammajk_r_m")
save("elambda_egammajk_r_m_6perso.mat", "elambda_egammajk_r_m")


%%
% Combine all the information from the two paths to get total sensitivities

% Upload the sensitivities
load("dvecUtilde_dgamma_m_6perso.mat");
load("dvecFtilde_dgamma_m_6perso.mat");

load("dlambda_dgammaik_b_m_6perso.mat")%, "dlambda_dsigmajk_u_m");
load("dlambda_dgammajk_r_m_6perso.mat")%, "dlambda_dsigmajk_r_m");

load("elambda_egammaik_b_m_6perso.mat")%, "elambda_esigmajk_u_m");
load("elambda_egammajk_r_m_6perso.mat")%, "elambda_esigmajk_r_m");



% Combine them along the two paths ( This is where it doesn't work - not the same dimensions)
% Need to change arrangement.
% dvecUtilde_dgamma.mat is arranged stage within age within personality (ik)
% dvecFtilde_dgamma.mat is arranged age within stage within personality (jk)
% We want to rearrange dvecUtilde_dgamma.mat so it becomes age within stage within personality (jk)
% Transform into a table

% transform dlambda_dgammajk_r in a vector
dlambda_dgammajk_m = dlambda_dgammajk_r_m{1}';
elambda_egammajk_m = elambda_egammajk_r_m{1}';
for mat = 2:s/s1
    tmps= dlambda_dgammajk_r_m{mat}';
    tmpe= elambda_egammajk_r_m{mat}';
    dlambda_dgammajk_m(end+1: end+w) = tmps;
    elambda_egammajk_m(end+1: end+w) = tmpe;
end

% transform dlambda_dbetaik_b in a vector
dlambda_dgammaik_m = dlambda_dgammaik_b_m{1}';
elambda_egammaik_m = elambda_egammaik_b_m{1}';
for mat = 2:s/s2
    tmps = dlambda_dgammaik_b_m{mat}';
    tmpe = elambda_egammaik_b_m{mat}';
    dlambda_dgammaik_m(end+1: end+b) = tmps;
    elambda_egammaik_m(end+1: end+b) = tmpe;
end

% Ok, so now put this in a table
% Arrangement jk (age within stage within personality)
age_classes = [1:w];
stage_classes = [1:b];

Personality_jk = repelem(1:g, b*w)';
Stage_jk = repmat(repelem(stage_classes, w),1,g)';
Age_jk = repmat(age_classes, 1,b*g)';
Sens_jk = dlambda_dgammajk_m;
Elas_jk = elambda_egammajk_m;

Tjk = table(Personality_jk, Stage_jk, Age_jk, Sens_jk, Elas_jk);

% Arrangement ik (stage within age within personality)
Personality_ik = repelem(1:g, b*w)';
Stage_ik = repmat(stage_classes, 1, w*g)';
Age_ik = repmat(repelem(age_classes, b)',g,1);
Sens_ik = dlambda_dgammaik_m;
Elas_ik = elambda_egammaik_m;

Tik = table(Personality_ik, Stage_ik, Age_ik, Sens_ik, Elas_ik);


% Sort Tik to follow arrangement of Tjk
Tik_jk = sortrows(Tik, {'Personality_ik','Stage_ik','Age_ik'});

% Sum across dimensions
Sens_lambda_gamma_Utilde = table2array(Tik_jk(:,4));
Sens_lambda_gamma_mtilde = table2array(Tjk(:,4));
Sens_lambda_gamma = Sens_lambda_gamma_Utilde + Sens_lambda_gamma_mtilde;
save('Sens_lambda_gamma_m_6perso.mat', "Sens_lambda_gamma");

Elas_lambda_gamma_Utilde = table2array(Tik_jk(:,5));
Elas_lambda_gamma_mtilde = table2array(Tjk(:,5));
Elas_lambda_gamma = Elas_lambda_gamma_Utilde + Elas_lambda_gamma_mtilde;
save('Elas_lambda_gamma_m_6perso.mat', "Elas_lambda_gamma");




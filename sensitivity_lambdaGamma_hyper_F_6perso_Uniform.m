% Sensitivity analyses  (run model first)
% Authors: Joanie Van de Walle, Silke van Daalen

%% SENSITIVITY OF LAMBDA TO BREEDING SUCCESS PROBABILITY, GAMMA
%%

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

[w_f,lam_f,v_f]=eig(A_tilde_f);
lambda_f = diag(lam_f);
[I_f,~]=find(lam_f==max(diag(lam_f)));
lambda1_f=lambda_f(I_f);
v_f=v_f(:,I_f);

% Change w here
w_f = new_w_vec;
%w_f=w_f(:,I_f);
dlambda_dvecAtilde_f=kron(w_f',v_f')/(v_f'*w_f);

%% PATH 1 - U TILDE
%%
%  4 dvecUtilde_dgamma = dvecUtilde_dvecBblock*dvecBblock_dgamma
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
totalsum_b_f = zeros((w*b*g)^2,b);
sens_310_b_f = {};
dlambda_dgammaik_b_f= {}; % This will give the sensitivity of lambda to gamma for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
elambda_egammaik_b_f = {}; % This will give the elasticity of lambda to  gamma for reproductive states 1-6 for each of the combination of age and personality (310 matrices)
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
     betaik = betaiks_f{i};
     gammaik = gammaiks_f{i};
     
     BETA_plus = B_plus.*(Oneb*betaik);
     BETA_minus = B_minus.*(Eb-(Oneb*betaik));
    
     BETAS = BETA_plus + BETA_minus;
     dvecBETAS = diag(BETAS(:));
     parta = dvecBETAS*Ib1b;
     
     dvecYg = diag(Yg(:));
     partb = dvecYg*Ib1b;
     
     %dvecB_dvecgammaik = parta - partb;
     dvecBik_dgamma_f = parta - partb;
     %insummation = bigthing*dvecB_dvecgammaik;
     dvecBblock_dgamma_f = dvecBblock_dvecBik_f*dvecBik_dgamma_f;
     
     sens_310_b_f{i} = dvecBblock_dgamma_f;
     
     totalsum_b_f = totalsum_b_f + dvecBblock_dgamma_f;
     
     dlambda_dgammaik_b_f{i} = dlambda_dvecAtilde_f*dvecUtilde_dvecBblock_f*dvecBblock_dgamma_f;
     elambda_egammaik_b_f{i} = (dlambda_dgammaik_b_f{i}.*gammaiks_f{i})/lambda1_f;
    
end

    
%%%% Finalise calculation for the path between sigma and U
dvecBblock_dgamma_f = totalsum_b_f;
dvecUtilde_dgamma_f = dvecUtilde_dvecBblock_f*dvecBblock_dgamma_f;


% Save the output
save('dvecUtilde_dgamma_f_6perso_Uniform.mat', "dvecUtilde_dgamma_f");
save("dlambda_dgammaik_b_f_6perso_Uniform.mat", "dlambda_dgammaik_b_f")
save("elambda_egammaik_b_f_6perso_Uniform.mat", "elambda_egammaik_b_f")



%% PATH 2 - F TILDE
%%
% 9 dvecFtilde_dgamma = dvecFtilde_dvecRblock*dvecRblock_dgamma
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
totalsum_r_f = zeros((w*b*g)^2,31);
sens_60_r_f = {};
dlambda_dgammajk_r_f= {}; % This will give the sensitivity of lambda to gamma for ages 1-31 for each of the combination of reproductive state and personality
elambda_egammajk_r_f = {}; % This will give the elasticity of lambda to gamma for for ages 1-31 for each of the combination of reproductive state and personality

% Loop over all the s/s1 (60) matrices
for i=1:s/s1
 E =  Emat(i,i,s/s1);
 % Then apply the vec transformation
 vecE = E(:);
 
 % Then kronecker product with the Is1square matrix
 vecEI = kron(sparse(vecE), sparse(Is1square));
 
% Then multiply them
%bigthing = dvecRblock_dvecRjk*vecEI;
dvecRblock_dvecRjk_f = bigthing*vecEI; % this is equation 12

indiag = (rho*Yr).*(Onew*betajks_f{i}').*(Onew*sigmajks_f{i}');
dvec_indiag = diag(indiag(:));

 dvecRjk_dgamma_f = dvec_indiag*Iw1w; % this is equation 13
 dvecRblock_dgamma_f = dvecRblock_dvecRjk_f*dvecRjk_dgamma_f;

sens_60_r_f{i} = dvecRblock_dgamma_f;

totalsum_r_f = totalsum_r_f + dvecRblock_dgamma_f;

    dlambda_dgammajk_r_tmp_f = dlambda_dvecAtilde_f*dvecFtilde_dvecRblock_f*dvecRblock_dgamma_f;
    dlambda_dgammajk_r_f{i} = dlambda_dgammajk_r_tmp_f;
    elambda_egammajk_r_f{i} = (dlambda_dgammajk_r_tmp_f.*gammajks_f{i}')/lambda1_f;
end


%%%% Finalise calculation for the path between sigma and F
dvecRblock_dgamma_f = totalsum_r_f;
dvecFtilde_dgamma_f = dvecFtilde_dvecRblock_f*dvecRblock_dgamma_f;

% Save the output
save('dvecFtilde_dgamma_f_6perso_Uniform.mat', "dvecFtilde_dgamma_f");
save("dlambda_dgammajk_r_f_6perso_Uniform.mat", "dlambda_dgammajk_r_f")
save("elambda_egammajk_r_f_6perso_Uniform.mat", "elambda_egammajk_r_f")


%%
% Combine all the information from the two paths to get total sensitivities

% Upload the sensitivities
load("dvecUtilde_dgamma_f_6perso_Uniform.mat");
load("dvecFtilde_dgamma_f_6perso_Uniform.mat");

load("dlambda_dgammaik_b_f_6perso_Uniform.mat")%, "dlambda_dsigmajk_u_f");
load("dlambda_dgammajk_r_f_6perso_Uniform.mat")%, "dlambda_dsigmajk_r_f");

load("elambda_egammaik_b_f_6perso_Uniform.mat")%, "elambda_esigmajk_u_f");
load("elambda_egammajk_r_f_6perso_Uniform.mat")%, "elambda_esigmajk_r_f");



% Combine them along the two paths ( This is where it doesn't work - not the same dimensions)
% Need to change arrangement.
% dvecUtilde_dgamma.mat is arranged stage within age within personality (ik)
% dvecFtilde_dgamma.mat is arranged age within stage within personality (jk)
% We want to rearrange dvecUtilde_dgamma.mat so it becomes age within stage within personality (jk)
% Transform into a table

% transform dlambda_dgammajk_r in a vector
dlambda_dgammajk_f = dlambda_dgammajk_r_f{1}';
elambda_egammajk_f = elambda_egammajk_r_f{1}';
for mat = 2:s/s1
    tmps= dlambda_dgammajk_r_f{mat}';
    tmpe= elambda_egammajk_r_f{mat}';
    dlambda_dgammajk_f(end+1: end+w) = tmps;
    elambda_egammajk_f(end+1: end+w) = tmpe;
end

% transform dlambda_dbetaik_b in a vector
dlambda_dgammaik_f = dlambda_dgammaik_b_f{1}';
elambda_egammaik_f = elambda_egammaik_b_f{1}';
for mat = 2:s/s2
    tmps = dlambda_dgammaik_b_f{mat}';
    tmpe = elambda_egammaik_b_f{mat}';
    dlambda_dgammaik_f(end+1: end+b) = tmps;
    elambda_egammaik_f(end+1: end+b) = tmpe;
end

% Ok, so now put this in a table
% Arrangement jk (age within stage within personality)
age_classes = [1:w];
stage_classes = [1:b];

Personality_jk = repelem(1:g, b*w)';
Stage_jk = repmat(repelem(stage_classes, w),1,g)';
Age_jk = repmat(age_classes, 1,b*g)';
Sens_jk = dlambda_dgammajk_f;
Elas_jk = elambda_egammajk_f;

Tjk = table(Personality_jk, Stage_jk, Age_jk, Sens_jk, Elas_jk);

% Arrangement ik (stage within age within personality)
Personality_ik = repelem(1:g, b*w)';
Stage_ik = repmat(stage_classes, 1, w*g)';
Age_ik = repmat(repelem(age_classes, b)',g,1);
Sens_ik = dlambda_dgammaik_f;
Elas_ik = elambda_egammaik_f;

Tik = table(Personality_ik, Stage_ik, Age_ik, Sens_ik, Elas_ik);


% Sort Tik to follow arrangement of Tjk
Tik_jk = sortrows(Tik, {'Personality_ik','Stage_ik','Age_ik'});

% Sum across dimensions
Sens_lambda_gamma_Utilde = table2array(Tik_jk(:,4));
Sens_lambda_gamma_Ftilde = table2array(Tjk(:,4));
Sens_lambda_gamma = Sens_lambda_gamma_Utilde + Sens_lambda_gamma_Ftilde;
save('Sens_lambda_gamma_f_6perso_Uniform.mat', "Sens_lambda_gamma");

Elas_lambda_gamma_Utilde = table2array(Tik_jk(:,5));
Elas_lambda_gamma_Ftilde = table2array(Tjk(:,5));
Elas_lambda_gamma = Elas_lambda_gamma_Utilde + Elas_lambda_gamma_Ftilde;
save('Elas_lambda_gamma_f_6perso_Uniform.mat', "Elas_lambda_gamma");




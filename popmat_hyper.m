function [a, b2, siz, s, A_tilde] = popmat_hyper(m, w, b, g, sigmas, betas, ...
                                                                          gammas, rho, tau, edges_e, ...
                                                                           midpoint_e, Vp, intercept, PERSONALITY, SSD_init)
% Follow steps in Roth and Caswell 2006

%%%%%%%%
% 1.1, specify the number of dimensions
%%%%%%%%
m=m;

%%%%%%%%
% 1.2, list the k-stages in the model
%%%%%%%%
% dimension 1 represents age classes
% dimension 2 represents breeding states
% dimension 3 reprents personality classes

 w = w; %number of age classes
 b = b;  %number of breeding states
 g = g; %number of personality classes

%%%%%%%%
% 1.3, create the vector siz
%%%%%%%%
siz = [w,b,g];
s = w*b*g;

%useful matrices
Is=speye(length(siz));
es=ones(length(siz),1);
e1=eye([w,1]);

%%%%%%%%
% 1.4, create the k-stage-transition matrix for each dimension.
%%%%%%%%

% Define the index
% i = 1:w ; i is for age
% j = 1:b ; j is for breeding state
% k = 1:g ; k is for personality class

% There are 3 sets of matrices, one for each process: 1) Ageing (U), 2) Breeding (B), 3) Change in Personality, (P)

%_____________________
% Build the U matrices
%_____________________
% This is the aging process. U matrices are w X w. There is one U matrix per breeding state and personality combination. So U can be indexed Ujk
Ujk = zeros(w,w,b,g); % initiate empty cell array 

% Loop. This will create 1 matrix for each of the 6 breeding state of personality 1, and then repeat for personality 2, and so on
% For a total of 6X10 = 60 matrices
for k = 1:g % for each personality class
    for j = 1:b % for each breeding state
     surv=sigmas(:,j,k)';
     U_temp = diag(surv(1:(w-1)),-1);
     U_temp(w,w) = surv(w);
     Ujk(:,:,j,k)=U_temp;
    end   
end

%_____________________
% Build the B matrices
%_____________________

% This is the breeding process. B matrices are b X b. There is one B matrix per age and personality combination. So B can be indexed Bik
Bik = zeros(b,b,w,g); % initiate empty cell array 

% Loop. This will create 1 matrix for each of the 32 age classes of personality 1, and then repeat for personality 2, and so on
% For a total of 31X10 = 310 matrices
for k = 1:g
    for i = 1:w
        B_temp = zeros(b,b);
        %Pre-breeders
        B_temp(1,1) = (1-betas(i,1,k));
        B_temp(2,1) = betas(i,1,k)*gammas(i,1,k);
        B_temp(3,1) = betas(i,1,k)*(1-gammas(i,1,k));
        % Successful breeders
        B_temp(2,2) = betas(i,2,k)*gammas(i,2,k);
        B_temp(3,2) = betas(i,2,k)*(1-gammas(i,2,k));
        B_temp(4,2) = (1-betas(i,2,k));
        % Failed breeders
        B_temp(2,3) = betas(i,3,k)*gammas(i,3,k);
        B_temp(3,3) = betas(i,3,k)*(1-gammas(i,3,k));
        B_temp(5,3) = (1-betas(i,3,k));
        % Post-successful breeders
        B_temp(2,4) = betas(i,4,k)*gammas(i,4,k);
        B_temp(3,4) = betas(i,4,k)*(1-gammas(i,4,k));
        B_temp(6,4) = (1-betas(i,4,k));
        % Post-failed breeders
        B_temp(2,5) = betas(i,5,k)*gammas(i,5,k);
        B_temp(3,5) = betas(i,5,k)*(1-gammas(i,5,k));
        B_temp(6,5) = (1-betas(i,5,k));
        % Non-breeders
        B_temp(2,6) = betas(i,6,k)*gammas(i,6,k);
        B_temp(3,6) = betas(i,6,k)*(1-gammas(i,6,k));
        B_temp(6,6) = (1-betas(i,6,k));
        
        Bik(:,:,i,k) = B_temp;
    end
end

%_____________________
% Build the P matrices
%_____________________

% This is the change in personality process. P matrices are g X g. There is one P matrix per age and state combination. So P can be indexed Pij
Pij = zeros(g,g,w,b); % initiate empty cell array 

% Loop. This will create 1 matrix for each of the 32 age classes of breeding state 1, and then repeat for breeding state 2, and so on
% For a total of 31X6 = 186 matrices
for j = 1:b
    for i = 1:w
        P_temp = eye(g); % here we assume individuals do not change personality score between years
        Pij(:,:,i,j) = P_temp;
    end
end


a{1} = Ujk;
a{2} = Bik;
a{3} = Pij;



%_____________________
% Build the R matrices
%_____________________

% Here I add a new parameter, rho, which is the sex-ratio at birth
rho=rho;

% This is age-specific reproduction process. R matrices are w x w. There is one R matrix per breeding state and personality score combination. So R can be indexed Rjk
Rjk = zeros(w,w,b,g); % initiate the cell array

% Loop. This will create 1 matrix for each of the 6 breeding states of personality score 1, and then repeat for personality score 2, and so on
% For a total of 6X10 = 70 matrices
for k = 1:g
    for j = 1:b
        R_temp = zeros(w,w);
        R_temp(1,:) = sigmas(:,j,k).*betas(:,j,k).*gammas(:,j,k)*rho;
        Rjk(:,:,j,k) = R_temp;
    end
end


%_____________________
% Build the F matrices
%_____________________


% This is breeding state specific fertilities. F matrices are b x b. There is one F matrix per age class and personality score combination. So F can be indexed Fik
Fik = zeros(b,b,w,g);

% Loop. This will create 1 matrix for each of the 32 age classes of personality score 1, and then repeat for personality score 2, and so on
% For a total of 31X10 = 310 matrices
for k = 1:g
    for i =1:w
       F_temp = zeros(b,b);
       F_temp(1,:) = repelem(1,b);
       Fik(:,:,i,k) = F_temp;
    end
end


%_____________________
% Build the H matrices
%_____________________
slope=tau/2;

Hij = zeros(g,g,w,b);
% for j = 1:b
%     for i = 1:w
%         H_temp = zeros(g,g);
%             
%         % Here recalculate the proportion of the population when removing one personality class at a time (because the diagonal has to be fixed)
%         H1 = eye(g);
%         %tryout= rand(10,1);
%         %tryout = tryout/sum(tryout);
%         %H2 = repmat(tryout, 1, g);
%         H2 = repmat(SSD_init(:,3), 1, g);
%         H3 = H2-diag(diag(H2));
%         H4 = H3./sum(H3);
%         H_temp = tau*H1 + (1-tau)*H4;
%         
%         Hij(:,:,i,j) = H_temp;
% end
    
%Create the H mat. This is the transmission between mother and offspring phenotype. This is based on the mother daughter 
%hmat = h_matrix(tau, g, edges_e, midpoint_e, Vp, intercept);
binwidth_e = midpoint_e(2)-midpoint_e(1);
hmat = exp(-(midpoint_e'-(intercept+slope*midpoint_e)).^2./(2*Vp)).*binwidth_e./sum(exp(-(midpoint_e'-(intercept + slope*midpoint_e)).^2./(2*Vp)).*binwidth_e,1);

% [BVi,BVj] = meshgrid(PERSONALITY,PERSONALITY);
%  hmat = exp(-(BVi-(intercept+slope*BVj)).^2./(2*Vp))./sum(exp(-(BVi-(intercept+slope*BVj)).^2./(2*Vp)),1);
%   
% hmat = (1/(sqrt(Vp)*sqrt(2*pi)))*exp(-(midpoint_e'-midpoint_e).^2./(2*Vp))./sum((1/(sqrt(Vp)*sqrt(2*pi)))*exp(-(midpoint_e'-midpoint_e).^2./(2*Vp)),1);

 %imagesc(hmat)

% Loop. This will create 1 matrix for each of the 32 age classes of breeding state 1, and then repeat for breeding state 2, and so on
% For a total of 31X6 = 186 matrices
for j = 1:b
    for i = 1:w
        H_temp = hmat;        
        Hij(:,:,i,j) = H_temp;
    end
end

b2{1} = Rjk;
b2{2} = Fik;
b2{3} = Hij;

% Get the block transition matrices
for i=1:m
    U{i}=BD_proj_mat(a{i});
    F{i}=BD_proj_mat(b2{i});
end

% Get to projection matrix
tildeU=hyper_state_matrix(siz,U);
tildeF=hyper_state_matrix(siz,F);

A_tilde = tildeU + tildeF;

end



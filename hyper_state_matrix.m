;%hyper_state_matrix function
%
%Given
%   - the size vector siz
%   - the cell array A containing the block-diagonal projection matrices
%
%Returns
%   - out.tildA         the total projection matrix tildA 
%
%08/02/16
function hsm=hyper_state_matrix(siz,A)
m=length(siz);          %number of dimensions
s=prod(siz);

K=eye(s);   
tildA=A{1}; 

%computation of the product in equation (23)
for k=2:m
    tildA=A{k}*vecperm_hyp(k,siz)*tildA;
    K=vecperm_hyp(k,siz)*K;                 %the matrix K keeps track of 
                                            %the product of the 
                                           %vec-permutation matrices 
end
hsm=K'*tildA;
end
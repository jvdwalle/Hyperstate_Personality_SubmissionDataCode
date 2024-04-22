function [Isk, Issk, Isk2, Kperm] = Imat(sk, s)

% This function gets the I and K2 (Kperm) matrices in equation 34 from Roth and Caswell (2016) for a specified dimension

% Ratio s/sk
ssk=s/sk;

% Calculate the various I mat
Isk = eye(sk);
Issk = eye(ssk);
Isk2 = eye(sk*sk);

% vec permutation matrix, on two dimensions, for dimension 1, and with dimensions sk, s/sk
% vec permutation matrix
 p = zeros(sk*ssk);
 a = zeros(sk,ssk);
 for i = 1:sk
     for j = 1:ssk
         e = a;
         e(i,j) = 1;
         p = p + kron(e,e');
     end
 end
Kperm = p;


end


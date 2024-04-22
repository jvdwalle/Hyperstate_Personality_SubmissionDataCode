% vecperm_hyp function
%
% Given:
%   - sizes vector s=(s1, ...,sm)
%   - index k
%
% Returns 
%   - the kth vec-permutation matrix associated with the sizes s1, ..., sm
%
% 24/7/15

function p = vecperm_hyp(k,s)
[~,m]=size(s);
if k==1
    p=eye(prod(s));
elseif k == m
    p = Qmat(prod(s(1:m-1)),s(m))*kron(eye(prod(s(m))),Qmat(prod(s(1:m-2)),s(m-1)))';
else
    p=kron(eye(prod(s(k+1:m))),Qmat(prod(s(1:k-1)),s(k)))*kron(eye(prod(s(k:m))),Qmat(prod(s(1:k-2)),s(k-1)))';
end


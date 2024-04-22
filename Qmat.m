%Qmat function
%
%Given:
%   - sizes u and v
%
%Returns:
%   - the matrix Q of size uv x uv (defined in eq (24) of the 
%     hyper-state model's paper)
%
%24/09/14

function k= Qmat(u,v)

k=zeros(u*v);
a=zeros(u,v);
for i=1:u
    for j = 1:v
        e=a;
        e(i,j)=1;
        k = k+ kron(e,e');
    end
end

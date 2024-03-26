% First order difference
% For piecewise constant
function mat = mat1(psix,T,eps)

m = size(psix,2);
[~,index] = sort(T);
spsix= psix(index,:);
R = (diff(spsix))'*diff(spsix);
R = R+eps*eye(size(m));
mat = R;

end
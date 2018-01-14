function F_es = fundamental_matrix(x1_test, x2_test)

[x1_new, T1] = normalise2dpts(x1_test);
[x2_new, T2] = normalise2dpts(x2_test);

u1=x1_new(1,:)';
v1=x1_new(2,:)';

u2=x2_new(1,:)';
v2=x2_new(2,:)';

W=[u1.*u2 v1.*u2 u2 u1.*v2 v1.*v2 v2 u1 v1 ones(size(u1,1),1)];

[U,D,V] = svd(W);
F = reshape(V(:,9),3,3)';

[Uf,Df,Vf] = svd(F);
Df(3,3)=0;
F_es=Uf*Df*Vf';
F_es=T2'*F_es*T1;

end

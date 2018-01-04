function [H]=homography2d(x1_points,x2_points)
    [x1norm, T1]=normalization(x1_points);
    [x2norm, T2]=normalization(x2_points);
    x2 = x2norm(1,:);
    y2 = x2norm(2,:);
    w2 = x2norm(3,:);
    A=zeros(size(x1norm,2)*2,9);
    for i=1:size(x1norm,2)
        A(i*2-1,:)=[zeros(1,3)     -w2(i)*x1norm(:,i)'   y2(i)*x1norm(:,i)'];
        A(i*2,:)=[w2(i)*x1norm(:,i)'   zeros(3,1)'     -x2(i)*x1norm(:,i)'];
    end
    [U,D,V] = svd(A);
    H = reshape(V(:,9),3,3)';
    H = inv(T2) * H * T1;

end

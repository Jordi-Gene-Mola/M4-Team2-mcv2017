function [ X ] = triangulate( x1, x2, P1, P2, imsize )


    x1=homog(x1);
    x2=homog(x2);
        
    [x1, x2, P1, P2]=normalization(x1,x2,P1,P2,imsize);

    x1 = euclid(x1);
    x2 = euclid(x2);
    
    A=[ x1(1)*P1(3,:)-P1(1,:);
        x1(2)*P1(3,:)-P1(2,:);
        x2(1)*P2(3,:)-P2(1,:);
        x2(2)*P2(3,:)-P2(2,:)];
    
    
    [U,D,V] = svd(A);
    
    X=V(:,4)./V(4,4);

end


function [x1n, x2n, P1n, P2n]=normalization(x1,x2,P1,P2,imsize)
    H = [2/imsize(1), 0, -1;
        0, 2/imsize(2), -1;
        0, 0, 1;];
    x1n=H*x1;
    x2n=H*x2;
    P1n=H*P1;
    P2n=H*P2;
end

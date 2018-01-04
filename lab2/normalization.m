function [new_x, T] = normalization( x )
    pts(1,:)=x(1,:)./x(3,:);
    pts(2,:)=x(2,:)./x(3,:);
    pts(3,:)=1;
        
    c = mean(pts(1:2,:),2);            % Centroid of finite points
    newp(1,:) = pts(1,:)-c(1); % Shift origin to centroid.
    newp(2,:) = pts(2,:)-c(2);
    
    dist = sqrt(newp(1,:).^2 + newp(2,:).^2);
    meandist = mean(dist(:));  % Ensure dist is a column vector for Octave 3.0.1
    
    scale = sqrt(2)/meandist;
    
    T = [scale   0   -scale*c(1)
         0     scale -scale*c(2)
         0       0      1      ];
    
    new_x = T*pts;
end